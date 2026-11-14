//! Compiles opencore-amr's AMR-NB codec as C99 (as its own `--enable-compile-c` build does).
//! Sources come from `.cache/opencore/` (run `just opencore`); the list mirrors amrnb/Makefile.am,
//! plus opencore's VAD option 2 (`vad2`, `r_fft`), which it ships but doesn't build (oracle for
//! the Rust port's VAD2; the C encoder itself stays VAD1).
use std::path::PathBuf;

const DEC: &[&str] = &[
    "agc",
    "amrdecode",
    "a_refl",
    "b_cn_cod",
    "bgnscd",
    "c_g_aver",
    "d1035pf",
    "d2_11pf",
    "d2_9pf",
    "d3_14pf",
    "d4_17pf",
    "d8_31pf",
    "dec_amr",
    "dec_gain",
    "dec_input_format_tab",
    "dec_lag3",
    "dec_lag6",
    "d_gain_c",
    "d_gain_p",
    "d_plsf_3",
    "d_plsf_5",
    "d_plsf",
    "dtx_dec",
    "ec_gains",
    "ex_ctrl",
    "if2_to_ets",
    "int_lsf",
    "lsp_avg",
    "ph_disp",
    "post_pro",
    "preemph",
    "pstfilt",
    "qgain475_tab",
    "sp_dec",
    "wmf_to_ets",
];
const ENC: &[&str] = &[
    "amrencode",
    "div_32",
    "l_abs",
    "vad1",
    "autocorr",
    "c1035pf",
    "c2_11pf",
    "c2_9pf",
    "c3_14pf",
    "c4_17pf",
    "c8_31pf",
    "calc_cor",
    "calc_en",
    "cbsearch",
    "cl_ltp",
    "cod_amr",
    "convolve",
    "cor_h",
    "cor_h_x2",
    "cor_h_x",
    "corrwght_tab",
    "dtx_enc",
    "enc_lag3",
    "enc_lag6",
    "enc_output_format_tab",
    "ets_to_if2",
    "ets_to_wmf",
    "g_adapt",
    "gain_q",
    "g_code",
    "g_pitch",
    "hp_max",
    "inter_36",
    "inter_36_tab",
    "lag_wind",
    "lag_wind_tab",
    "l_comp",
    "levinson",
    "l_extract",
    "lflg_upd",
    "l_negate",
    "lpc",
    "ol_ltp",
    "pitch_fr",
    "pitch_ol",
    "p_ol_wgh",
    "pre_big",
    "pre_proc",
    "prm2bits",
    "qgain475",
    "qgain795",
    "q_gain_c",
    "q_gain_p",
    "qua_gain",
    "s10_8pf",
    "set_sign",
    "sid_sync",
    "sp_enc",
    "spreproc",
    "spstproc",
    "ton_stab",
];
const COMMON: &[&str] = &[
    "add",
    "az_lsp",
    "bitno_tab",
    "bitreorder_tab",
    "c2_9pf_tab",
    "div_s",
    "extract_h",
    "extract_l",
    "gains_tbl",
    "gc_pred",
    "get_const_tbls",
    "gmed_n",
    "gray_tbl",
    "grid_tbl",
    "int_lpc",
    "inv_sqrt",
    "inv_sqrt_tbl",
    "l_deposit_h",
    "l_deposit_l",
    "log2",
    "log2_norm",
    "log2_tbl",
    "lsfwt",
    "l_shr_r",
    "lsp_az",
    "lsp",
    "lsp_lsf",
    "lsp_lsf_tbl",
    "lsp_tab",
    "mult_r",
    "negate",
    "norm_l",
    "norm_s",
    "overflow_tbl",
    "ph_disp_tab",
    "pow2",
    "pow2_tbl",
    "pred_lt",
    "r_fft",
    "q_plsf_3",
    "q_plsf_3_tbl",
    "q_plsf_5",
    "q_plsf_5_tbl",
    "q_plsf",
    "qua_gain_tbl",
    "reorder",
    "residu",
    "round",
    "set_zero",
    "shr",
    "shr_r",
    "sqrt_l",
    "sqrt_l_tbl",
    "sub",
    "syn_filt",
    "vad2",
    "weight_a",
    "window_tab",
];

/// The gitignored `.cache/` of this checkout, or of the main checkout when building in a
/// `git worktree` (which has no `.cache/` of its own).
fn cache_dir() -> PathBuf {
    let ws = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
    if ws.join(".cache/opencore").exists() {
        return ws.join(".cache");
    }
    let out = std::process::Command::new("git")
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

fn main() {
    let root = std::env::var_os("OPENCORE_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| cache_dir().join("opencore/opencore-amr-0.1.6"));
    assert!(
        root.join("amrnb/wrapper.cpp").exists(),
        "opencore sources not found at {}: run `just opencore`",
        root.display()
    );
    let amr = root.join("opencore/codecs_v2/audio/gsm_amr");
    let nb = amr.join("amr_nb");

    let mut b = cc::Build::new();
    b.flag("-x")
        .flag("c")
        .flag("-std=c99")
        .flag("-w")
        .opt_level(3);
    if std::env::var_os("CARGO_FEATURE_PIN_UB").is_some() {
        b.flag("-ftrivial-auto-var-init=zero").flag("-fwrapv");
    }
    for inc in [
        root.join("oscl"),
        nb.join("dec/src"),
        nb.join("common/include"),
        nb.join("dec/include"),
        amr.join("common/dec/include"),
        nb.join("enc/src"),
        nb.join("enc/include"),
        root.join("amrnb"),
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("compat"),
    ] {
        b.include(inc);
    }
    b.file(root.join("amrnb/wrapper.cpp"));
    b.file(PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("shim.c"));
    b.file(PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("compat/ops.c"));
    for (dir, files) in [("dec/src", DEC), ("enc/src", ENC), ("common/src", COMMON)] {
        for f in files {
            b.file(nb.join(dir).join(format!("{f}.cpp")));
        }
    }
    b.compile("opencore_amrnb");
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed=shim.c");
    println!("cargo:rerun-if-changed=compat");
    println!("cargo:rerun-if-env-changed=OPENCORE_DIR");
}
