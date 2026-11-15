//! AMR-WB RFC 4867 readers on arbitrary bytes: no panics, and every frame they accept
//! survives a write → read round trip in the canonical (minimal) encoding.
//!
//! Input: one selector byte, then the body. Selector % 3: 0 = `.awb` file (the magic is
//! prepended unless bit 2 is set, to reach `BadMagic`), 1 = RTP octet-aligned,
//! 2 = RTP bandwidth-efficient.
#![no_main]

use amr_wb::rfc4867::{
    Error, Frame, Packing, STORAGE_MAGIC, read_rtp, read_storage, write_rtp, write_storage,
};
use libfuzzer_sys::fuzz_target;

const PACKINGS: [Packing; 2] = [Packing::OctetAligned, Packing::BandwidthEfficient];

fn check_frame(f: &Frame) {
    assert!(f.bits.iter().all(|&b| b == 0 || b == 1), "non-binary bit");
}

fn storage(file: &[u8]) {
    let Ok(frames) = read_storage(file) else {
        assert!(!file.starts_with(STORAGE_MAGIC));
        return;
    };
    let mut ok = Vec::new();
    let mut ended = false;
    for r in frames {
        assert!(!ended, "item after a terminal error");
        match r {
            Ok(f) => {
                check_frame(&f);
                ok.push(f);
            }
            Err(Error::UnsupportedFrameType(10..=13) | Error::Truncated) => ended = true,
            Err(e) => panic!("unexpected {e:?}"),
        }
    }
    let out = write_storage(&ok).expect("read frames are writable");
    // Each accepted frame took the same octets in the input as in the canonical output.
    assert!(out.len() <= file.len());
    let back: Result<Vec<_>, _> = read_storage(&out).unwrap().collect();
    assert_eq!(back.unwrap(), ok, "storage round trip");
}

fn rtp(payload: &[u8], packing: Packing) {
    let Ok((cmr, frames)) = read_rtp(payload, packing) else {
        return;
    };
    let frames: Vec<Frame> = frames.collect();
    assert!(!frames.is_empty());
    frames.iter().for_each(check_frame);
    for p in PACKINGS {
        let out = write_rtp(cmr, &frames, p).expect("read frames are writable");
        if p == packing {
            assert!(out.len() <= payload.len(), "canonical form is not longer");
        }
        let (c, f) = read_rtp(&out, p).expect("written payload reads back");
        assert_eq!(
            (c, f.collect::<Vec<_>>()),
            (cmr, frames.clone()),
            "{p:?} round trip"
        );
    }
}

fuzz_target!(|data: &[u8]| {
    let Some((&sel, body)) = data.split_first() else {
        return;
    };
    match sel % 3 {
        0 if sel & 4 != 0 => storage(body),
        0 => storage(&[STORAGE_MAGIC.as_slice(), body].concat()),
        1 => rtp(body, Packing::OctetAligned),
        _ => rtp(body, Packing::BandwidthEfficient),
    }
});
