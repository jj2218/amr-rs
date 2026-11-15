//! Stress streams for differential tests: real vector frames with flipped bits, random and
//! all-ones bits, every frame type, mode switching. NO_DATA frames are generated only while
//! the decoder is in DTX (after a SID frame), keeping out CR 26.173-0027..29's rule, the one
//! deliberate difference from the reference (it only applies to NO_DATA during speech).

use super::{Frame, cases};
use amr_wb::{MAX_SERIAL_SIZE, Mode, RxFrameType, SID_BITS};

/// xorshift64*, deterministic per seed.
pub struct Rng(pub u64);

impl Rng {
    pub fn seeded(seed: u64) -> Self {
        Self(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1)
    }
    pub fn next(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }
    pub fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
    pub fn chance(&mut self, percent: u64) -> bool {
        self.below(100) < percent
    }
}

/// The C decoders' form of the bits (BIT_0 = -127, BIT_1 = 127).
pub fn to_prms(bits: &[i16; MAX_SERIAL_SIZE]) -> [i16; MAX_SERIAL_SIZE] {
    bits.map(|b| if b != 0 { 127 } else { -127 })
}

/// Which frame types and bits a stream draws from.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct Mix {
    /// Bad, lost and degraded speech frames.
    pub bad_frames: bool,
    /// SID and NO_DATA frames.
    pub dtx: bool,
    /// Random and all-ones bits (else real frames with a few flipped bits).
    pub random_bits: bool,
}

impl Mix {
    pub const ALL: Mix = Mix {
        bad_frames: true,
        dtx: true,
        random_bits: true,
    };
}

/// Real good frames per mode (the single-mode vectors): the source of realistic bits.
pub fn pool() -> [Vec<Frame>; 9] {
    let cases = cases();
    std::array::from_fn(|m| {
        let dir = format!(
            "T_{}/",
            [
                "660", "885", "1265", "1425", "1585", "1825", "1985", "2305", "2385"
            ][m]
        );
        cases
            .iter()
            .filter(|c| c.name.starts_with(&dir) && !c.name.contains("Dtx"))
            .flat_map(|c| {
                c.frames
                    .iter()
                    .filter(|f| f.frame_type == RxFrameType::SpeechGood)
                    .take(200)
                    .cloned()
            })
            .collect()
    })
}

/// One stream of `n` frames.
pub fn stream(rng: &mut Rng, pool: &[Vec<Frame>; 9], n: usize, mix: Mix) -> Vec<Frame> {
    use RxFrameType::*;
    let switching = rng.chance(30);
    let mut mode = rng.below(9) as usize;
    let mut in_dtx = false;
    let mut dtx_left = 0;
    let mut out = Vec::with_capacity(n);
    while out.len() < n {
        if switching {
            mode = rng.below(9) as usize;
        }
        let mut frame_type = if dtx_left > 0 {
            dtx_left -= 1;
            match rng.below(100) {
                0..=79 => NoData,
                80..=91 => SidUpdate,
                92..=95 => SidBad,
                _ => SpeechLost,
            }
        } else {
            match rng.below(100) {
                0..=69 => SpeechGood,
                70..=74 => SpeechProbablyDegraded,
                75..=79 => SpeechBad,
                80..=84 => SpeechLost,
                85..=89 => SidBad,
                _ => {
                    dtx_left = rng.below(80);
                    SidFirst
                }
            }
        };
        if !mix.dtx && matches!(frame_type, SidFirst | SidUpdate | SidBad | NoData) {
            dtx_left = 0;
            frame_type = SpeechGood;
        }
        if !mix.bad_frames && matches!(frame_type, SpeechBad | SpeechLost | SpeechProbablyDegraded)
        {
            frame_type = SpeechGood;
        }
        // NO_DATA only in DTX (tracked as the decoder does)
        if frame_type == NoData && !in_dtx {
            frame_type = if mix.bad_frames {
                SpeechLost
            } else {
                SpeechGood
            };
        }
        in_dtx = match frame_type {
            SidFirst | SidUpdate | SidBad => true,
            NoData | SpeechBad | SpeechLost => in_dtx,
            SpeechGood | SpeechProbablyDegraded => false,
        };
        let nbits = match frame_type {
            SidFirst | SidUpdate | SidBad => SID_BITS,
            NoData => 0,
            _ => Mode::from_index(mode).unwrap().bits(),
        };
        let mut bits = [0i16; MAX_SERIAL_SIZE];
        let src = &pool[mode][rng.below(pool[mode].len() as u64) as usize];
        match if mix.random_bits { rng.below(10) } else { 0 } {
            0..=5 => {
                bits[..nbits].copy_from_slice(&src.bits[..nbits]);
                for _ in 0..rng.below(4) {
                    if nbits > 0 {
                        let i = rng.below(nbits as u64) as usize;
                        bits[i] ^= 1;
                    }
                }
            }
            6..=8 => bits[..nbits]
                .iter_mut()
                .for_each(|b| *b = (rng.next() >> 40) as i16 & 1),
            _ => bits[..nbits].fill(1),
        }
        out.push(Frame {
            frame_type,
            mode: Mode::from_index(mode).unwrap(),
            bits,
            prms: to_prms(&bits),
        });
    }
    out
}
