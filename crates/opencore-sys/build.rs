//! Compiles opencore-amr's codecs as C99 (as its own `--enable-compile-c` build does), one
//! library per feature. Sources come from `.cache/opencore/` (run `just opencore`).
//! - `nb`: AMR-NB; the list mirrors amrnb/Makefile.am, plus opencore's VAD option 2 (`vad2`,
//!   `r_fft`), which it ships but doesn't build (oracle for the Rust port's VAD2; the C encoder
//!   itself stays VAD1).
//! - `wb`: the AMR-WB decoder; the list mirrors amrwb/Makefile.am, plus `shim_wb.c`.
//! - `wbenc`: the AMR-WB encoder, vo-amrwbenc 0.1.3 (`.cache/opencore/vo-amrwbenc-0.1.3`), plus
//!   `shim_wbenc.c`.
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
const WB_DEC: &[&str] = &[
    "agc2_amr_wb",
    "band_pass_6k_7k",
    "dec_acelp_2p_in_64",
    "dec_acelp_4p_in_64",
    "dec_alg_codebook",
    "dec_gain2_amr_wb",
    "deemphasis_32",
    "dtx_decoder_amr_wb",
    "get_amr_wb_bits",
    "highpass_400hz_at_12k8",
    "highpass_50hz_at_12k8",
    "homing_amr_wb_dec",
    "interpolate_isp",
    "isf_extrapolation",
    "isp_az",
    "isp_isf",
    "lagconceal",
    "low_pass_filt_7k",
    "median5",
    "mime_io",
    "noise_gen_amrwb",
    "normalize_amr_wb",
    "oversamp_12k8_to_16k",
    "phase_dispersion",
    "pit_shrp",
    "pred_lt4",
    "preemph_amrwb_dec",
    "pvamrwbdecoder",
    "pvamrwb_math_op",
    "q_gain2_tab",
    "qisf_ns",
    "qisf_ns_tab",
    "qpisf_2s",
    "qpisf_2s_tab",
    "scale_signal",
    "synthesis_amr_wb",
    "voice_factor",
    "wb_syn_filt",
    "weight_amrwb_lpc",
];
/// vo-amrwbenc 0.1.3's `amrwbenc/src` list (Makefile.am, without the ARM assembly).
const WB_ENC: &[&str] = &[
    "autocorr",
    "az_isp",
    "bits",
    "c2t64fx",
    "c4t64fx",
    "convolve",
    "cor_h_x",
    "decim54",
    "deemph",
    "dtx",
    "g_pitch",
    "gpclip",
    "homing",
    "hp400",
    "hp50",
    "hp6k",
    "hp_wsp",
    "int_lpc",
    "isp_az",
    "isp_isf",
    "lag_wind",
    "levinson",
    "log2",
    "lp_dec2",
    "math_op",
    "mem_align",
    "oper_32b",
    "p_med_ol",
    "pit_shrp",
    "pitch_f4",
    "pred_lt4",
    "preemph",
    "q_gain2",
    "q_pulse",
    "qisf_ns",
    "qpisf_2s",
    "random",
    "residu",
    "scale",
    "stream",
    "syn_filt",
    "updt_tar",
    "util",
    "voAMRWBEnc",
    "voicefac",
    "wb_vad",
    "weight_a",
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
    let here = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let build = || {
        let mut b = cc::Build::new();
        b.flag("-x")
            .flag("c")
            .flag("-std=c99")
            // Legacy C: quiet, but what GCC >= 14 and Clang >= 16 reject fails here too (not `-w`,
            // which also hides -Werror; cc adds -Wall -Wextra before these flags).
            .flag("-Wno-all")
            .flag("-Wno-extra")
            .flag_if_supported("-Wno-builtin-declaration-mismatch")
            .flag("-Werror=implicit-function-declaration")
            .flag("-Werror=implicit-int")
            .flag("-Werror=int-conversion")
            .flag("-Werror=incompatible-pointer-types")
            .opt_level(3);
        if std::env::var_os("CARGO_FEATURE_PIN_UB").is_some() {
            b.flag("-ftrivial-auto-var-init=zero").flag("-fwrapv");
        }
        b
    };

    if std::env::var_os("CARGO_FEATURE_NB").is_some() {
        let nb = amr.join("amr_nb");
        let mut b = build();
        for inc in [
            root.join("oscl"),
            nb.join("dec/src"),
            nb.join("common/include"),
            nb.join("dec/include"),
            amr.join("common/dec/include"),
            nb.join("enc/src"),
            nb.join("enc/include"),
            root.join("amrnb"),
            here.join("compat"),
        ] {
            b.include(inc);
        }
        b.flag("-include")
            .flag(here.join("compat/ops.h").to_str().unwrap());
        b.file(root.join("amrnb/wrapper.cpp"));
        b.file(here.join("shim.c"));
        b.file(here.join("compat/ops.c"));
        for (dir, files) in [("dec/src", DEC), ("enc/src", ENC), ("common/src", COMMON)] {
            for f in files {
                b.file(nb.join(dir).join(format!("{f}.cpp")));
            }
        }
        b.compile("opencore_amrnb");
    }

    if std::env::var_os("CARGO_FEATURE_WB").is_some() {
        let wb = amr.join("amr_wb/dec");
        let mut b = build();
        for inc in [
            root.join("oscl"),
            wb.join("src"),
            wb.join("include"),
            amr.join("common/dec/include"),
            root.join("amrwb"),
        ] {
            b.include(inc);
        }
        b.file(root.join("amrwb/wrapper.cpp"));
        b.file(here.join("shim_wb.c"));
        for f in WB_DEC {
            b.file(wb.join("src").join(format!("{f}.cpp")));
        }
        b.compile("opencore_amrwb");
    }
    if std::env::var_os("CARGO_FEATURE_WBENC").is_some() {
        let vo = std::env::var_os("VO_AMRWBENC_DIR")
            .map(PathBuf::from)
            .unwrap_or_else(|| cache_dir().join("opencore/vo-amrwbenc-0.1.3"));
        assert!(
            vo.join("enc_if.h").exists(),
            "vo-amrwbenc sources not found at {}: run `just opencore`",
            vo.display()
        );
        let mut b = build();
        b.include(vo.join("amrwbenc/inc"))
            .include(vo.join("common/include"));
        b.file(vo.join("common/cmnMemory.c"));
        b.file(here.join("shim_wbenc.c"));
        for f in WB_ENC {
            b.file(vo.join("amrwbenc/src").join(format!("{f}.c")));
        }
        b.compile("vo_amrwbenc");
    }
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed=shim.c");
    println!("cargo:rerun-if-changed=shim_wb.c");
    println!("cargo:rerun-if-changed=shim_wbenc.c");
    println!("cargo:rerun-if-env-changed=VO_AMRWBENC_DIR");
    println!("cargo:rerun-if-changed=compat");
    println!("cargo:rerun-if-env-changed=OPENCORE_DIR");
}
