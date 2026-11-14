//! G.711 A-law/µ-law for the TS 26.074 companded vectors (clause 9): a test-side helper, not
//! a port. Linear samples are 16-bit left-justified, as in the codec's `.INP`/`.OUT` files
//! (the encoder homing value 0x0008 maps to A-law 0xD5 / µ-law 0xFE and back). The 8-bit codes
//! are stored one per 16-bit word.

/// Sample encoding of a test vector.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Law {
    Linear,
    A,
    Mu,
}

impl Law {
    /// The 8-bit code in each word of a companded file → 16-bit linear; linear passes through.
    pub fn expand(self, s: i16) -> i16 {
        match self {
            Law::Linear => s,
            Law::A => alaw_expand(s as u8),
            Law::Mu => ulaw_expand(s as u8),
        }
    }

    /// 16-bit linear → the 8-bit code as the files store it; linear passes through.
    pub fn compress(self, s: i16) -> i16 {
        match self {
            Law::Linear => s,
            Law::A => alaw_compress(s) as i16,
            Law::Mu => ulaw_compress(s) as i16,
        }
    }
}

/// A-law: 13-bit magnitude (input >> 3) in 8 segments, even bits inverted (XOR 0x55).
fn alaw_compress(s: i16) -> u8 {
    let mut ix = if s < 0 { !s >> 4 } else { s >> 4 }; // 12-bit magnitude, halved step
    if ix > 15 {
        let mut seg = 1;
        while ix > 16 + 15 {
            ix >>= 1;
            seg += 1;
        }
        ix = ix - 16 + (seg << 4);
    }
    let sign = if s >= 0 { 0x80 } else { 0 };
    (ix as u8 | sign) ^ 0x55
}

fn alaw_expand(c: u8) -> i16 {
    let ix = (c ^ 0x55) & 0x7f;
    let seg = ix >> 4;
    let mut m = (ix & 0x0f) as i16;
    if seg > 0 {
        m += 16;
    }
    m = (m << 4) + 8; // segment midpoint
    if seg > 1 {
        m <<= seg - 1;
    }
    if c & 0x80 != 0 { m } else { -m }
}

/// µ-law: biased (+33) 14-bit magnitude in 8 segments, all bits inverted.
fn ulaw_compress(s: i16) -> u8 {
    let mag = if s < 0 {
        (!s >> 2) as i32
    } else {
        (s >> 2) as i32
    } + 33;
    let mag = mag.min(0x1fff);
    let seg = 32 - (mag >> 6).leading_zeros() as i32 + 1; // 1..=8
    let low = 0x0f - ((mag >> seg) & 0x0f);
    let code = ((0x08 - seg) << 4 | low) as u8;
    if s >= 0 { code | 0x80 } else { code }
}

fn ulaw_expand(c: u8) -> i16 {
    let m = !c;
    let exp = ((m >> 4) & 0x07) as i32;
    let step = 4 << (exp + 1);
    let mag = (0x80 << exp) + step * (m & 0x0f) as i32 + step / 2 - 4 * 33;
    (if c < 0x80 { -mag } else { mag }) as i16
}

#[test]
fn homing_value_and_round_trip() {
    for (law, code) in [(Law::A, 0xd5), (Law::Mu, 0xfe)] {
        assert_eq!(law.expand(code), 8);
        assert_eq!(law.compress(8), code);
        // Every code survives expand → compress (µ-law's two zeros, 0x7F and 0xFF, map to +0).
        for c in 0..=255i16 {
            let back = law.compress(law.expand(c));
            assert!(
                back == c || (law == Law::Mu && c == 0x7f && back == 0xff),
                "{law:?} {c}"
            );
        }
    }
    // Saturation and sign at the extremes.
    assert_eq!(Law::A.compress(i16::MAX), 0xaa);
    assert_eq!(Law::A.compress(i16::MIN), 0x2a);
    assert_eq!(Law::Mu.compress(i16::MAX), 0x80);
    assert_eq!(Law::Mu.compress(i16::MIN), 0x00);
}
