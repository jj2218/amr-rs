//! Shared helpers for differential tests against the C oracle.

/// Deterministic xorshift64* generator: reproducible inputs without extra dependencies.
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Self {
        Self(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1)
    }
    pub fn u64(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }
    /// Mostly uniform i16, with a bias towards edge values (0, ±1, extremes, powers of two).
    pub fn i16(&mut self) -> i16 {
        match self.u64() % 8 {
            0 => [0, 1, -1, i16::MAX, i16::MIN, 0x4000, -0x4000, 0x7ffe][(self.u64() % 8) as usize],
            1 => (1i16 << (self.u64() % 15)).wrapping_mul(if self.u64() & 1 == 0 { 1 } else { -1 }),
            _ => self.u64() as i16,
        }
    }
    pub fn i32(&mut self) -> i32 {
        match self.u64() % 8 {
            0 => [
                0,
                1,
                -1,
                i32::MAX,
                i32::MIN,
                0x4000_0000,
                -0x4000_0000,
                0x7fff_ffff - 1,
            ][(self.u64() % 8) as usize],
            1 => (1i32 << (self.u64() % 31)).wrapping_mul(if self.u64() & 1 == 0 { 1 } else { -1 }),
            _ => self.u64() as i32,
        }
    }
    /// Small shift amount in -40..=40 (covers saturation and out-of-range branches).
    pub fn shift(&mut self) -> i16 {
        (self.u64() % 81) as i16 - 40
    }
}

/// Number of random cases per differential test.
pub const N: usize = 200_000;

/// Runs `f` on the TS 26.073 reference (`ref26073-sys`) holding its lock, with its global
/// overflow flag cleared first and ORed into `o` after (unless null): the sticky opencore-style
/// `Flag *pOverflow` the differential tests were written against.
#[cfg(test)]
#[allow(unsafe_code)]
pub fn with_ref<T>(o: *mut i32, f: impl FnOnce() -> T) -> T {
    let _l = ref26073_sys::lock();
    ref26073_sys::clear_overflow();
    let r = f();
    if !o.is_null() {
        // SAFETY: callers pass a valid pointer or null.
        unsafe { *o |= ref26073_sys::overflow() as i32 };
    }
    r
}
