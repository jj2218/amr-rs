//! Compiles the 3GPP TS 26.173 reference C code (AMR-WB fixed point) from
//! `.cache/3gpp/ref26173/c-code` (run `just vectors`), the source list being its makefile's
//! (`DOBJ`, without the `decoder.c` main). Then every global symbol the library defines is
//! renamed `refwb_<name>` (`nm` + `objcopy --redefine-syms`), so it links next to opencore and
//! the TS 26.073 reference (`ref_`), which define some of the same names. 3GPP code is "all
//! rights reserved": an oracle only.
use std::path::PathBuf;
use std::process::Command;

/// The gitignored `.cache/` of this checkout, or of the main checkout when building in a
/// `git worktree` (as `opencore-sys`).
fn cache_dir() -> PathBuf {
    let ws = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
    if ws.join(".cache/3gpp/ref26173").exists() {
        return ws.join(".cache");
    }
    let out = Command::new("git")
        .args(["rev-parse", "--path-format=absolute", "--git-common-dir"])
        .current_dir(&ws)
        .output()
        .ok()
        .filter(|o| o.status.success());
    match out {
        Some(o) => PathBuf::from(String::from_utf8_lossy(&o.stdout).trim())
            .parent()
            .map_or_else(|| ws.join(".cache"), |main| main.join(".cache")),
        None => ws.join(".cache"),
    }
}

/// `tool` of the C compiler's binutils: its cross prefix when cross-compiling (the host's
/// `objcopy` leaves an aarch64 archive unchanged, without an error), else plain `tool`.
fn binutil(b: &cc::Build, tool: &str) -> String {
    let cc = b.get_compiler();
    let name = cc.path().file_name().unwrap().to_string_lossy();
    match name.rfind('-') {
        Some(i) => format!("{}{tool}", &name[..=i]),
        None => tool.to_string(),
    }
}

fn run(cmd: &mut Command) -> Vec<u8> {
    let out = cmd.output().unwrap_or_else(|e| panic!("{cmd:?}: {e}"));
    assert!(
        out.status.success(),
        "{cmd:?}: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    out.stdout
}

fn main() {
    let src = std::env::var_os("REF26173_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| cache_dir().join("3gpp/ref26173/c-code"));
    let makefile = std::fs::read_to_string(src.join("makefile.gcc")).unwrap_or_else(|e| {
        panic!(
            "3GPP TS 26.173 code not found at {}: run `just vectors` ({e})",
            src.display()
        )
    });
    // `DOBJ =  decoder.o agc2.o ... \` continued lines
    let start = makefile.find("DOBJ =").expect("DOBJ") + "DOBJ =".len();
    let mut objs = Vec::new();
    for line in makefile[start..].lines() {
        objs.extend(
            line.split_whitespace()
                .filter_map(|w| w.strip_suffix(".o"))
                .filter(|&o| o != "decoder")
                .map(String::from),
        );
        if !line.trim_end().ends_with('\\') {
            break;
        }
    }
    let mut b = cc::Build::new();
    // typedefs.h only knows 32-bit `linux && i386` as little-endian PC: the same types here.
    // -fwrapv: basicop2.c detects overflow after a signed add (UB), as 2001 compilers wrapped.
    b // Legacy C: quiet, but what GCC >= 14 and Clang >= 16 reject fails here too (not `-w`,
        // which also hides -Werror; cc adds -Wall -Wextra before these flags).
        .flag("-Wno-all")
        .flag("-Wno-extra")
        .flag_if_supported("-Wno-builtin-declaration-mismatch")
        .flag("-Werror=implicit-function-declaration")
        .flag("-Werror=implicit-int")
        .flag("-Werror=int-conversion")
        .flag("-Werror=incompatible-pointer-types")
        .flag("-fwrapv")
        .opt_level(2)
        .include(&src)
        .define("linux", None)
        .define("i386", None)
        .define("WMOPS", "0");
    for f in &objs {
        b.file(src.join(format!("{f}.c")));
    }
    b.compile("ref26173");

    let out = PathBuf::from(std::env::var_os("OUT_DIR").unwrap());
    let lib = out.join("libref26173.a");
    let nm = run(Command::new(binutil(&b, "nm"))
        .args(["--defined-only", "-g", "-P"])
        .arg(&lib));
    let mut map = String::new();
    let mut seen = std::collections::BTreeSet::new();
    for line in String::from_utf8_lossy(&nm).lines() {
        // archive member headers end in ':'; symbol lines are "name type value size"
        if let Some(name) = line.split_whitespace().next().filter(|n| !n.ends_with(':'))
            && seen.insert(name.to_string())
        {
            map += &format!("{name} refwb_{name}\n");
        }
    }
    let syms = out.join("ref26173.syms");
    std::fs::write(&syms, map).unwrap();
    run(Command::new(binutil(&b, "objcopy"))
        .arg(format!("--redefine-syms={}", syms.display()))
        .arg(&lib));

    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-env-changed=REF26173_DIR");
    // a scratch copy with debug prints (REF26173_DIR) rebuilds when edited
    println!("cargo:rerun-if-changed={}", src.display());
}
