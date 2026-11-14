//! Table layout for the lane searches: codebooks and gain tables stored component-major, and
//! the "first minimum over the lanes" reduction they share.

/// A table stored component-major: `cols[k][i]` is component `k` of row `i`, so a block of
/// consecutive rows is one contiguous load per component.
pub(crate) type Soa<const D: usize, const N: usize> = [[i16; N]; D];

/// The first `D` components of `N` rows of `table` (a row every `stride` values) as a [`Soa`].
pub(crate) const fn soa<const D: usize, const N: usize>(table: &[i16], stride: usize) -> Soa<D, N> {
    let mut t = [[0i16; N]; D];
    let mut i = 0;
    while i < N {
        let mut k = 0;
        while k < D {
            t[k][i] = table[stride * i + k];
            k += 1;
        }
        i += 1;
    }
    t
}

/// A [`Soa`]'s columns as slices (the dispatched kernels take any table size).
pub(crate) fn cols<const D: usize, const N: usize>(t: &'static Soa<D, N>) -> [&'static [i16]; D] {
    t.each_ref().map(|c| &c[..])
}

/// The C's "first strict minimum" from per-lane results: lane `l` kept the first strict
/// minimum `best[l]` (at `idx[l]`) of the candidates it saw in order, starting from `MAX_32` and
/// index 0. The overall first minimum is the smallest distance, ties to the lowest index; 0 if
/// no candidate was below `MAX_32` (lanes never updated keep index 0). Lane-width independent.
#[inline(always)]
pub(crate) fn first_min<const W: usize>(best: &[i32; W], idx: &[i32; W]) -> i32 {
    let (mut dist_min, mut index) = (i32::MAX, 0i32);
    for l in 0..W {
        if best[l] < dist_min || (best[l] == dist_min && idx[l] < index) {
            (dist_min, index) = (best[l], idx[l]);
        }
    }
    index
}
