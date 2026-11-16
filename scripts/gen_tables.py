"""Generate a crate's src/tables.rs from opencore-amr's C table files (Apache-2.0).

Deterministic and mechanical: every `const WordNN name[...] = { ... };` (also `int16`, 2-D
arrays flattened row-major) becomes `pub(crate) static NAME: [iNN; N]`. Fails loudly on
anything it cannot parse. Run via `just gen-tables`; never edit tables.rs by hand.
Usage: gen_tables.py nb|wb <gsm_amr/amr_nb or amr_wb> <out tables.rs>
       gen_tables.py wbenc <vo-amrwbenc-0.1.3/amrwbenc> <out enc_tables.rs> <decoder tables.rs>
(wbenc: vo-amrwbenc's file-level tables, also non-`const` ones; a table the decoder's tables.rs
already has with the same values is not repeated; ambiguous names get the file stem.)
"""

import pathlib
import re
import sys

CODEC = sys.argv[1]
ROOT = pathlib.Path(sys.argv[2])  # .../gsm_amr/amr_nb or amr_wb
OUT = pathlib.Path(sys.argv[3])
if CODEC == "nb":
    FILES = sorted(
        p
        for d in ("common/src", "dec/src", "enc/src")
        for p in (ROOT / d).glob("*.cpp")
        if re.search(r"(_tab|_tbl)\.cpp$", p.name)
    )
    # File-level tables outside the *_tab files (VAD2, its FFT, the MR475/515 track table);
    # listed by name, since these files also hold function-local arrays.
    EXTRA = {
        "vad2.cpp": {"ch_tbl_sh", "vm_tbl", "hangover_table", "burstcount_table", "vm_threshold_table"},
        "r_fft.cpp": {"phs_tbl", "ii_table"},
        "c2_9pf.cpp": {"trackTable"},
    }
    FILES += [ROOT / "common/src/vad2.cpp", ROOT / "common/src/r_fft.cpp", ROOT / "enc/src/c2_9pf.cpp"]
elif CODEC == "wb":
    # The decoder keeps its tables at file level, also outside the *_tab files (no
    # function-local initialised arrays: the generator would pick them up too).
    FILES = sorted((ROOT / "dec/src").glob("*.cpp")) + [ROOT / "dec/include/pvamrwbdecoder_api.h"]
    EXTRA = {}
elif CODEC == "wbenc":
    FILES = sorted((ROOT / "inc").glob("*.tab")) + [ROOT / "inc/log2_tab.h"] + sorted((ROOT / "src").glob("*.c"))
    EXTRA = {}
    # Names that repeat across files or say nothing on their own get the file stem.
    RENAME = {("hp6k.c", "fir_6k_7k"): "hp6k_fir_6k_7k", ("hp50.c", "a"): "hp50_a", ("hp50.c", "b"): "hp50_b", ("hp400.c", "a"): "hp400_a",
              ("hp400.c", "b"): "hp400_b", ("hp_wsp.c", "a"): "hp_wsp_a", ("hp_wsp.c", "b"): "hp_wsp_b",
              ("log2_tab.h", "table"): "log2_table", ("q_gain2.c", "pred"): "q_gain2_pred",
              ("lp_dec2.c", "h_fir"): "lp_dec2_h_fir"}
    KNOWN = {m.group(1): m.group(2) for m in re.finditer(
        r"pub\(crate\) static (\w+): \[\w+; \d+\] = \[(.*?)\];", pathlib.Path(sys.argv[4]).read_text(), re.S)}
else:
    sys.exit("usage: gen_tables.py nb|wb|wbenc <dir> <out tables.rs> [decoder tables.rs]")


def tables(f):
    """The (type, name, body) table declarations of `f` to generate."""
    if CODEC == "wbenc":
        # file level: column 0 in .c files (function-local tables are indented), any in .tab
        decl = DECL_TAB if f.suffix in (".tab", ".h") else DECL_C
        return [(t[1], RENAME.get((f.name, t[2]), t[2]), t[3])
                for t in decl.findall(strip_comments(f.read_text(errors="replace")))]
    found = DECL.findall(strip_comments(f.read_text(errors="replace")))
    return [t[1:] for t in found if f.name not in EXTRA or t[2] in EXTRA[f.name]]


def is_static(f, name):
    """A `static` table has no linker symbol, so the C oracle check can't read it."""
    return any(t[0] and t[2] == name for t in DECL.findall(strip_comments(f.read_text(errors="replace"))))


TYPES = {"Word16": "i16", "Word32": "i32", "UWord8": "u8", "Word8": "i8", "int16": "i16", "int32": "i32", "uint8": "u8"}
CTYPES = "Word16|Word32|UWord8|Word8|int16|int32|uint8"
DECL = re.compile(rf"(static\s+)?\bconst\s+({CTYPES})\s+(\w+)\s*(?:\[[^\]]*\]\s*)+=\s*\{{(.*?)\}}\s*;", re.S)
DECL_C = re.compile(rf"^((?:static\s+)?(?:const\s+)?)({CTYPES})\s+(\w+)\s*(?:\[[^\]]*\]\s*)+=\s*\{{(.*?)\}}\s*;", re.S | re.M)
DECL_TAB = re.compile(rf"((?:static\s+)?(?:const\s+)?)\b({CTYPES})\s+(\w+)\s*(?:\[[^\]]*\]\s*)+=\s*\{{(.*?)\}}\s*;", re.S)
NUM = re.compile(rf"^\(?\s*(?:\((?:{CTYPES})\)\s*)?(-?\s*(?:0[xX][0-9a-fA-F]+|\d+))[lLuU]*\s*\)?$")


def strip_comments(s: str) -> str:
    s = re.sub(r"/\*.*?\*/", "", s, flags=re.S)
    return re.sub(r"//[^\n]*", "", s)


DEFINES = {}
for h in sorted(ROOT.rglob("*.h")) + sorted(ROOT.rglob("*.cpp")):
    for name, val in re.findall(r"^\s*#define\s+(\w+)\s+([^\n]+)", strip_comments(h.read_text(errors="replace")), re.M):
        DEFINES.setdefault(name, val.strip())


def eval_int(expr: str, where: str, depth: int = 0) -> int:
    """Evaluate an integer constant expression, expanding #defines; arithmetic only."""
    if depth > 20:
        sys.exit(f"gen_tables: #define recursion in {where}")
    expr = re.sub(rf"\((?:{CTYPES})\)", "", expr)
    expr = re.sub(r"\b(0[xX][0-9a-fA-F]+|\d+)[lLuU]+\b", r"\1", expr)
    for ident in set(re.findall(r"\b[A-Za-z_]\w*\b", expr)):
        if ident not in DEFINES:
            sys.exit(f"gen_tables: unknown identifier {ident!r} in {where}")
        expr = re.sub(rf"\b{ident}\b", f"({eval_int(DEFINES[ident], where, depth + 1)})", expr)
    if not re.fullmatch(r"[\s0-9xXa-fA-F+\-*/()<>]*", expr):
        sys.exit(f"gen_tables: non-arithmetic expression {expr!r} in {where}")
    return int(eval(expr.replace("/", "//"), {"__builtins__": {}}))  # noqa: S307 (arithmetic only)


def parse_value(tok: str, ty: str, where: str) -> str:
    m = NUM.match(tok.strip())
    v = int(m.group(1).replace(" ", ""), 0) if m else eval_int(tok, where)
    if ty == "i32" and v > 0x7FFFFFFF:
        v -= 1 << 32
    if ty == "i16" and v > 0x7FFF:
        v -= 1 << 16
    return str(v)


SOURCE = "vo-amrwbenc 0.1.3" if CODEC == "wbenc" else "opencore-amr 0.1.6"
out = [
    f"// @generated by scripts/gen_tables.py ({CODEC}) from {SOURCE} (Apache-2.0).",
    "// Do not edit: run `just gen-tables`. Names are the C names in upper case.",
    "#![allow(dead_code)]" if CODEC != "wbenc" else "#![allow(dead_code, unused_imports)]",
    "",
]
seen = {}
for f in FILES:
    for ty_c, name, body in tables(f):
        ty = TYPES[ty_c]
        where = f"{f.name}:{name}"
        flat = body.replace("\n", " ").replace("{", " ").replace("}", " ")
        vals = [parse_value(t, ty, where) for t in flat.split(",") if t.strip()]
        rust = name.upper()
        if CODEC == "wbenc" and rust in KNOWN:
            known = [v.strip() for v in KNOWN[rust].replace("\n", " ").split(",") if v.strip()]
            if known != vals:
                sys.exit(f"gen_tables: {where} differs from the decoder's {rust}: add it to RENAME")
            seen[rust] = f.name
            out.append(f"pub(crate) use crate::tables::{rust}; // `{name}` from `{f.name}`, as the decoder's")
            out.append("")
            continue
        if rust in seen:
            sys.exit(f"gen_tables: duplicate table name {rust} ({seen[rust]} and {f.name})")
        seen[rust] = f.name
        out.append(f"/// `{name}` from `{f.name}`.")
        out.append(f"pub(crate) static {rust}: [{ty}; {len(vals)}] = [")
        for i in range(0, len(vals), 12):
            out.append("    " + ", ".join(vals[i : i + 12]) + ",")
        out.append("];")
        out.append("")

# Oracle check: every generated table equals the C symbol in the linked opencore library.
decls, asserts = [], []
for f in FILES:
    for _, name, _ in tables(f):
        if is_static(f, name) or f.suffix == ".h" or CODEC == "wbenc":
            continue
        rust = name.upper()
        m = re.search(rf"pub\(crate\) static {rust}: \[(\w+); (\d+)\]", "\n".join(out))
        ty, n = m.groups()
        decls.append(f'        #[link_name = "{name}"]\n        static C_{rust}: [{ty}; {n}];')
        asserts.append(f'        assert_eq!(super::{rust}, unsafe {{ C_{rust} }}, "{name}");')
out += [
    "#[cfg(test)]",
    "#[allow(unsafe_code)]",
    "mod vs_c {",
    "    use opencore_sys as _; // links the C oracle",
    "    unsafe extern \"C\" {",
    *decls,
    "    }",
    "",
    "    #[test]",
    "    fn all_tables_match_c() {",
    *asserts,
    "    }",
    "}",
    "",
]
OUT.write_text("\n".join(out))
print(f"gen_tables: {len(seen)} tables from {len(FILES)} files -> {OUT}")
