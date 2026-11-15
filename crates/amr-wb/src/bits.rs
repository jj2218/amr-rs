//! Reading parameters from the serial bits (`get_amr_wb_bits.cpp`): one `i16` per bit, MSB
//! first, nonzero = 1.

pub(crate) struct Bits<'a> {
    bits: &'a [i16],
    pos: usize,
}

impl<'a> Bits<'a> {
    pub(crate) fn new(bits: &'a [i16]) -> Self {
        Self { bits, pos: 0 }
    }

    /// `Serial_parm`: the next `n` bits (n <= 15) as an unsigned value.
    #[inline]
    pub(crate) fn parm(&mut self, n: usize) -> i16 {
        let mut v: i16 = 0;
        for &b in &self.bits[self.pos..self.pos + n] {
            v = (v << 1) | (b != 0) as i16;
        }
        self.pos += n;
        v
    }

    /// `Serial_parm_1bit`.
    #[inline]
    pub(crate) fn bit(&mut self) -> i16 {
        self.parm(1)
    }
}
