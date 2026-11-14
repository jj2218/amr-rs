//! Compiles the 3GPP TS 26.073 reference C code (AMR-NB fixed point, VAD option 1) from
//! `.cache/3gpp/ref26073/c-code` (run `just vectors`), the source list being its makefile's
//! (`SPC_OBJS` + `FIPOP_SRCS`, no `coder.c`/`decoder.c` mains). Then every global symbol the
//! library defines is renamed `ref_<name>` (`nm` + `objcopy --redefine-syms`), so it links next
//! to opencore, which defines the same names. 3GPP code is "all rights reserved": an oracle only.
use std::path::PathBuf;
use std::process::Command;

/// The gitignored `.cache/` of this checkout, or of the main checkout when building in a
/// `git worktree` (as `opencore-sys`).
fn cache_dir() -> PathBuf {
    let ws = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
    if ws.join(".cache/3gpp/ref26073").exists() {
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
    let src = std::env::var_os("REF26073_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| cache_dir().join("3gpp/ref26073/c-code"));
    let makefile = std::fs::read_to_string(src.join("makefile")).unwrap_or_else(|e| {
        panic!(
            "3GPP TS 26.073 code not found at {}: run `just vectors` ({e})",
            src.display()
        )
    });
    // `SPC_OBJS= a.o b.o \` continued lines, then `FIPOP_SRCS=x.c y.c \` ...
    let list = |var: &str, ext: &str| -> Vec<String> {
        let start = makefile.find(&format!("{var}=")).expect(var) + var.len() + 1;
        let mut names = Vec::new();
        for line in makefile[start..].lines() {
            names.extend(
                line.split_whitespace()
                    .filter_map(|w| w.strip_suffix(ext))
                    .map(String::from),
            );
            if !line.trim_end().ends_with('\\') {
                break;
            }
        }
        names
    };
    let mut b = cc::Build::new();
    // typedefs.h only knows 32-bit `linux && i386` as little-endian PC: the same types here.
    // -fwrapv: basicop2.c detects overflow after a signed add (UB), as 1999 compilers wrapped.
    b.flag("-w")
        .flag("-fwrapv")
        .opt_level(2)
        .include(&src)
        .define("linux", None)
        .define("i386", None)
        .define("WMOPS", "0")
        .define("VAD1", None);
    for f in list("SPC_OBJS", ".o")
        .iter()
        .chain(&list("FIPOP_SRCS", ".c"))
    {
        b.file(src.join(format!("{f}.c")));
    }
    b.compile("ref26073");

    let out = PathBuf::from(std::env::var_os("OUT_DIR").unwrap());
    let lib = out.join("libref26073.a");
    let nm = run(Command::new("nm")
        .args(["--defined-only", "-g", "-P"])
        .arg(&lib));
    let mut map = String::new();
    let mut seen = std::collections::BTreeSet::new();
    for line in String::from_utf8_lossy(&nm).lines() {
        // archive member headers end in ':'; symbol lines are "name type value size"
        if let Some(name) = line.split_whitespace().next().filter(|n| !n.ends_with(':'))
            && seen.insert(name.to_string())
        {
            map += &format!("{name} ref_{name}\n");
        }
    }
    let syms = out.join("ref26073.syms");
    std::fs::write(&syms, map).unwrap();
    run(Command::new("objcopy")
        .arg(format!("--redefine-syms={}", syms.display()))
        .arg(&lib));

    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-env-changed=REF26073_DIR");
}
