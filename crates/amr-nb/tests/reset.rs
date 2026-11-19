//! `reset()` against a fresh instance (the oracle): stream X, `reset`, stream Y must equal a
//! new encoder or decoder on Y. The resets assign state field by field; a field left out
//! would leak X into Y. X drives the stateful paths (loud noise, mode switches, DTX hangover
//! and comfort noise, bad and lost frames); Y has speech-like input and silence. Modes switch
//! every frame, so mode-specific state (6.60/8.85 pitch-gain clipping) is read on speech.
//! Mutants (2026-11-18, one field dropped from a reset each): pre-processing and post-filter
//! states fail. Equivalent after a reset: `sid.reset()` (the VAD's reset hangover sends speech
//! frames first, which re-initialise the SID state) and the decoder's `prev_mode` (NO_DATA
//! concealment from an all-zero state is the same in every mode).
use amr_nb::{Decoder, Dtx, Encoder, L_FRAME, MAX_SERIAL_SIZE, Mode, RxFrameType};

/// xorshift64*: deterministic, no dependencies.
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_f491_4f6c_dd1d)
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next() >> 33) as usize % n
    }
}

/// Frames of a tone with noise at `level` (0: digital silence); `n` frames.
fn signal(rng: &mut Rng, n: usize, level: i32) -> Vec<[i16; L_FRAME]> {
    (0..n)
        .map(|f| {
            std::array::from_fn(|i| {
                let t = (f * L_FRAME + i) as f64;
                let tone = (t * 0.07).sin() * 0.6 + (t * 0.31).sin() * 0.3;
                let noise = (rng.below(2001) as f64 - 1000.0) / 1000.0 * 0.1;
                ((tone + noise) * level as f64) as i16
            })
        })
        .collect()
}

/// Stream X: loud and quiet noise (DTX); it ends quiet or loud.
fn disturb(rng: &mut Rng, end_loud: bool) -> Vec<[i16; L_FRAME]> {
    let (loud, quiet) = (signal(rng, 40, 32_000), signal(rng, 40, 30));
    if end_loud {
        [quiet, loud]
    } else {
        [loud, quiet]
    }
    .concat()
}

/// Stream Y: speech-like at a moderate level and silence (DTX), in either order.
fn clean(rng: &mut Rng, silence_first: bool) -> Vec<[i16; L_FRAME]> {
    let (speech, silence) = (signal(rng, 30, 8_000), signal(rng, 30, 0));
    if silence_first {
        [silence, speech]
    } else {
        [speech, silence]
    }
    .concat()
}

/// Every DTX setting; the last one feeds the decoder test.
const DTX: &[Dtx] = &[Dtx::Off, Dtx::Vad1, Dtx::Vad2];

fn mode(i: usize) -> Mode {
    Mode::ALL[i % Mode::ALL.len()]
}

#[test]
fn encoder_reset_equals_a_fresh_encoder() {
    let cases = DTX.iter().flat_map(|&d| [(d, false), (d, true)]);
    for (dtx, (end_loud, silence_first)) in
        cases.flat_map(|(d, e)| [(d, (e, false)), (d, (e, true))])
    {
        let mut rng = Rng(0x5eed_0001);
        let (x, y) = (disturb(&mut rng, end_loud), clean(&mut rng, silence_first));
        let mut enc = Encoder::new(dtx);
        for (i, f) in x.iter().enumerate() {
            enc.encode(mode(i), f);
        }
        enc.reset();
        let mut fresh = Encoder::new(dtx);
        for (i, f) in y.iter().enumerate() {
            let what = format!(
                "{dtx:?}, X ends loud: {end_loud}, Y silence first: {silence_first}, frame {i}"
            );
            assert_eq!(enc.encode(mode(i), f), fresh.encode(mode(i), f), "{what}");
        }
    }
}

#[test]
fn decoder_reset_equals_a_fresh_decoder() {
    use RxFrameType::*;
    let types = [
        SpeechGood,
        SpeechDegraded,
        Onset,
        SpeechBad,
        SidFirst,
        SidUpdate,
        SidBad,
        NoData,
    ];
    let mut rng = Rng(0x5eed_0002);
    let mut dec = Decoder::new();
    let mut pcm = [0i16; L_FRAME];
    // X: random bits and frame types (bad frames, comfort noise, concealment state).
    for i in 0..200 {
        let bits: [i16; MAX_SERIAL_SIZE] = std::array::from_fn(|_| (rng.next() & 1) as i16);
        let t = if i % 3 == 0 {
            types[rng.below(types.len())]
        } else {
            SpeechGood
        };
        dec.decode(t, mode(rng.below(16)), &bits, &mut pcm);
    }
    dec.reset();
    // Y: NO_DATA first (the previous mode is used), then an encoder's frames as a stream
    // really arrives.
    let mut enc = Encoder::new(DTX[DTX.len() - 1]);
    let mut fresh = Decoder::new();
    let (mut a, mut b) = ([0i16; L_FRAME], [0i16; L_FRAME]);
    let no_data = [0; MAX_SERIAL_SIZE];
    for _ in 0..3 {
        dec.decode(NoData, mode(1), &no_data, &mut a);
        fresh.decode(NoData, mode(1), &no_data, &mut b);
        assert_eq!(a, b, "NO_DATA before Y");
    }
    for (i, f) in clean(&mut rng, false).iter().enumerate() {
        let fr = enc.encode(Mode::MR122, f);
        dec.decode(fr.frame_type, fr.mode, &fr.bits, &mut a);
        fresh.decode(fr.frame_type, fr.mode, &fr.bits, &mut b);
        assert_eq!(a, b, "frame {i} of Y");
    }
}
