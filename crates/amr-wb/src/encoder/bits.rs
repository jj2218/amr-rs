//! The serial bit array the encoder writes (`Parm_serial`, vo-amrwbenc `bits.c`): one word
//! per bit, ±127 (`BIT_1`, `BIT_0`), MSB first.

/// `BIT_0`, `BIT_1`.
pub(crate) const BIT_0: i16 = -127;
pub(crate) const BIT_1: i16 = 127;

/// A cursor over the output bits (`Word16 **prms`).
pub(crate) struct Serial<'a> {
    out: &'a mut [i16],
    pos: usize,
}

impl<'a> Serial<'a> {
    pub(crate) fn new(out: &'a mut [i16]) -> Self {
        Self { out, pos: 0 }
    }

    /// `Parm_serial`: the `n` low bits of `value`, MSB first.
    pub(crate) fn parm(&mut self, value: i16, n: usize) {
        for k in 0..n {
            let bit = (value as i32 >> (n - 1 - k).min(15)) & 1;
            self.out[self.pos + k] = if bit == 0 { BIT_0 } else { BIT_1 };
        }
        self.pos += n;
    }

    /// Bits written so far.
    pub(crate) fn len(&self) -> usize {
        self.pos
    }
}
