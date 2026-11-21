//! RFC 4867 packed AMR-WB: the `.awb` storage format (§5) and the RTP payload format (§4),
//! octet-aligned and bandwidth-efficient, converted to and from the decoder's serial bits.
//!
//! Supported: single channel; RTP without interleaving, CRCs or robust sorting (the SDP
//! defaults). Speech bits are carried in the class-sorted order of TS 26.201 (opencore's
//! `mime_io.cpp` sorting tables); a SID frame (frame type 9) carries 35 comfort-noise bits,
//! the SID type indicator (STI) and a 4-bit mode indication (most significant bit first);
//! a SID_FIRST frame's comfort-noise bits are written as zeros, as the TS 26.173 reference's
//! `.awb` writer does (decoders don't use them).
//! Frame type 14 is SPEECH_LOST, 15 NO_DATA; 10..=13 are reserved.
//!
//! ```
//! use amr_wb::{Decoder, rfc4867};
//! # let file = rfc4867::write_storage(&[rfc4867::Frame::no_data()]).unwrap();
//! let mut dec = Decoder::new();
//! let mut pcm = [0i16; 320];
//! for f in rfc4867::read_storage(&file)? {
//!     let f = f?; // frames before a bad one are still delivered
//!     dec.decode(f.frame_type, f.mode, &f.bits, &mut pcm);
//! }
//! # Ok::<(), rfc4867::Error>(())
//! ```

use crate::{MAX_SERIAL_SIZE, Mode, RxFrameType, SID_BITS, tables};

/// Magic number that starts a single-channel `.awb` file.
pub const STORAGE_MAGIC: &[u8; 9] = b"#!AMR-WB\n";

const FT_SID: u8 = 9;
const FT_SPEECH_LOST: u8 = 14;
const FT_NO_DATA: u8 = 15;
/// SID frame: comfort-noise bits, STI, 4-bit mode indication.
const SID_FRAME_BITS: usize = SID_BITS + 1 + 4;
/// CMR value for "no mode request".
const CMR_NONE: u32 = 15;

/// One frame in the decoder's terms: the arguments of [`crate::Decoder::decode`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Frame {
    /// Receive type. Reading yields `SpeechGood`/`SpeechBad`, `SidFirst`/`SidUpdate`/`SidBad`
    /// (`*Bad` when the frame's quality bit Q is 0), `SpeechLost` or `NoData`.
    pub frame_type: RxFrameType,
    /// Codec mode. For SID frames the mode indication (the encoder's mode; `MD66` when it is
    /// not a mode, 9..=15), which the decoder ignores, as the TS 26.173 reference does: it
    /// uses the previous mode. Ignored (`MD66`) for `SpeechLost` and `NoData`.
    pub mode: Mode,
    /// Serial bits, one per element (0 or 1), in the decoder's (TS 26.173) order. Bits past the
    /// frame's size are zero when read and ignored when written.
    pub bits: [i16; MAX_SERIAL_SIZE],
}

impl Frame {
    /// A `NoData` frame (nothing received).
    pub fn no_data() -> Self {
        Self {
            frame_type: RxFrameType::NoData,
            mode: Mode::MD66,
            bits: [0; MAX_SERIAL_SIZE],
        }
    }
}

/// RTP payload layout, negotiated out of band (SDP `octet-align=1` or not).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Packing {
    /// Every header field and frame padded to whole octets (`octet-align=1`).
    OctetAligned,
    /// Fields and frames packed back to back (the RFC 4867 default).
    BandwidthEfficient,
}

/// Why packed data could not be read or written.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Error {
    /// The file doesn't start with `#!AMR-WB\n` (multi-channel files are not supported).
    BadMagic,
    /// The data ends inside a header or frame.
    Truncated,
    /// Frame type 10..=13 (reserved; size unknown). RFC 4867 §4.3.2: discard the packet. In a
    /// file, reading ends there.
    UnsupportedFrameType(u8),
    /// Writing: a `SpeechProbablyDegraded` frame, or an RTP packet with no frames; neither
    /// has a packed form.
    NotPackable,
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::BadMagic => f.write_str("not a single-channel AMR-WB file (no #!AMR-WB header)"),
            Error::Truncated => f.write_str("data ends inside a header or frame"),
            Error::UnsupportedFrameType(ft) => write!(f, "unsupported AMR-WB frame type {ft}"),
            Error::NotPackable => f.write_str("frame has no RFC 4867 representation"),
        }
    }
}

impl std::error::Error for Error {}

/// Reads a single-channel `.awb` file (magic, then one header octet + padded frame per frame)
/// as a stream of frames. A truncated frame or a frame of unknown size (types 10..=13) ends
/// the stream with its error.
///
/// ```no_run
/// use amr_wb::{Decoder, rfc4867};
///
/// let file = std::fs::read("speech.awb")?;
/// let mut dec = Decoder::new();
/// let mut pcm = [0i16; 320];
/// for frame in rfc4867::read_storage(&file)? {
///     let f = frame?; // truncated or reserved: the stream ends
///     dec.decode(f.frame_type, f.mode, &f.bits, &mut pcm);
/// }
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn read_storage(data: &[u8]) -> Result<impl Iterator<Item = Result<Frame, Error>> + '_, Error> {
    let body = data.strip_prefix(STORAGE_MAGIC).ok_or(Error::BadMagic)?;
    let mut r = BitReader { data: body, pos: 0 };
    let mut done = false;
    Ok(std::iter::from_fn(move || {
        if done || r.remaining() == 0 {
            return None;
        }
        let h = r.get(8);
        let ft = (h >> 3) as u8 & 15;
        let res = match frame_bits(ft) {
            Some(n) if r.remaining() < n => Err(Error::Truncated),
            Some(_) => Ok(read_frame(&mut r, ft, h & 4 != 0)),
            None => Err(Error::UnsupportedFrameType(ft)),
        };
        r.align();
        done = matches!(res, Err(Error::Truncated)) || frame_bits(ft).is_none();
        Some(res)
    }))
}

/// Writes a single-channel `.awb` file.
pub fn write_storage(frames: &[Frame]) -> Result<Vec<u8>, Error> {
    let mut w = BitWriter::default();
    w.out.extend_from_slice(STORAGE_MAGIC);
    for f in frames {
        let (ft, q) = frame_header(f)?;
        w.put(u32::from(ft) << 3 | u32::from(q) << 2, 8);
        write_frame(&mut w, f, ft);
        w.align();
    }
    Ok(w.out)
}

/// Reads one RTP payload: the codec mode request (`None` = no request; invalid values are
/// ignored, RFC 4867 §4.3.1) and the frames, in order. The whole payload is validated first;
/// octets after the last frame are ignored.
pub fn read_rtp(
    payload: &[u8],
    packing: Packing,
) -> Result<(Option<Mode>, impl Iterator<Item = Frame> + '_), Error> {
    let oa = packing == Packing::OctetAligned;
    let mut r = BitReader {
        data: payload,
        pos: 0,
    };
    let need = |r: &BitReader, n| {
        if r.remaining() < n {
            Err(Error::Truncated)
        } else {
            Ok(())
        }
    };
    need(&r, 4)?;
    let cmr = Mode::from_index(r.get(4) as usize);
    let mut toc = Vec::new();
    loop {
        if oa {
            r.align();
        }
        need(&r, 6)?;
        let e = r.get(6);
        let ft = (e >> 1) as u8 & 15;
        if frame_bits(ft).is_none() {
            return Err(Error::UnsupportedFrameType(ft));
        }
        toc.push((ft, e & 1 != 0));
        if e & 32 == 0 {
            break;
        }
    }
    let mut end = r.clone();
    for &(ft, _) in &toc {
        if oa {
            end.align();
        }
        let n = frame_bits(ft).unwrap_or(0);
        need(&end, n)?;
        end.pos += n;
    }
    Ok((
        cmr,
        toc.into_iter().map(move |(ft, q)| {
            if oa {
                r.align();
            }
            read_frame(&mut r, ft, q)
        }),
    ))
}

/// Writes one RTP payload carrying `frames` (at least one) and the codec mode request `cmr`
/// (`None`: no request).
///
/// ```
/// use amr_wb::{Mode, RxFrameType, rfc4867::{self, Frame, Packing}};
///
/// // Two 20 ms frames per packet, asking the far end to send 12.65 kbit/s.
/// let mut speech = Frame::no_data();
/// speech.frame_type = RxFrameType::SpeechGood;
/// speech.mode = Mode::MD2385;
/// let frames = [speech, Frame::no_data()];
/// let payload = rfc4867::write_rtp(Some(Mode::MD1265), &frames, Packing::BandwidthEfficient)?;
///
/// let (cmr, received) = rfc4867::read_rtp(&payload, Packing::BandwidthEfficient)?;
/// assert_eq!(cmr, Some(Mode::MD1265));
/// assert!(received.eq(frames));
/// # Ok::<(), rfc4867::Error>(())
/// ```
pub fn write_rtp(cmr: Option<Mode>, frames: &[Frame], packing: Packing) -> Result<Vec<u8>, Error> {
    if frames.is_empty() {
        return Err(Error::NotPackable);
    }
    let oa = packing == Packing::OctetAligned;
    let headers = frames
        .iter()
        .map(frame_header)
        .collect::<Result<Vec<_>, _>>()?;
    let mut w = BitWriter::default();
    w.put(cmr.map_or(CMR_NONE, |m| m as u32), 4);
    for (i, &(ft, q)) in headers.iter().enumerate() {
        if oa {
            w.align();
        }
        let more = i + 1 < headers.len();
        w.put(u32::from(more) << 5 | u32::from(ft) << 1 | u32::from(q), 6);
    }
    for (f, &(ft, _)) in frames.iter().zip(&headers) {
        if oa {
            w.align();
        }
        write_frame(&mut w, f, ft);
    }
    w.align();
    Ok(w.out)
}

/// Frame type (FT) and quality bit (Q) for `f`.
fn frame_header(f: &Frame) -> Result<(u8, bool), Error> {
    use RxFrameType::*;
    match f.frame_type {
        SpeechGood => Ok((f.mode as u8, true)),
        SpeechBad => Ok((f.mode as u8, false)),
        SidFirst | SidUpdate => Ok((FT_SID, true)),
        // Written with STI = 1 (a bad SID's type is unknown).
        SidBad => Ok((FT_SID, false)),
        SpeechLost => Ok((FT_SPEECH_LOST, true)),
        NoData => Ok((FT_NO_DATA, true)),
        SpeechProbablyDegraded => Err(Error::NotPackable),
    }
}

/// Speech bit `k` of the sorted payload is serial bit `SORT[ft][k]` (`mime_io.cpp`).
static SORT: [&[i16]; 9] = [
    &tables::SORT_660,
    &tables::SORT_885,
    &tables::SORT_1265,
    &tables::SORT_1425,
    &tables::SORT_1585,
    &tables::SORT_1825,
    &tables::SORT_1985,
    &tables::SORT_2305,
    &tables::SORT_2385,
];

/// Payload bits of frame type `ft` (`None`: reserved types 10..=13, size unknown).
fn frame_bits(ft: u8) -> Option<usize> {
    match ft {
        0..=8 => Some(SORT[ft as usize].len()),
        FT_SID => Some(SID_FRAME_BITS),
        FT_SPEECH_LOST | FT_NO_DATA => Some(0),
        _ => None,
    }
}

/// Reads a frame of type 0..=9, 14 or 15 whose `frame_bits` the caller has checked are present.
fn read_frame(r: &mut BitReader, ft: u8, q: bool) -> Frame {
    let mut f = Frame::no_data();
    match ft {
        0..=8 => {
            for &i in SORT[ft as usize] {
                f.bits[i as usize] = r.get(1) as i16;
            }
            f.mode = Mode::ALL[ft as usize];
            f.frame_type = if q {
                RxFrameType::SpeechGood
            } else {
                RxFrameType::SpeechBad
            };
        }
        FT_SID => {
            // `sort_SID` is the identity
            for b in &mut f.bits[..SID_BITS] {
                *b = r.get(1) as i16;
            }
            let sti = r.get(1);
            let mi = r.get(4);
            f.mode = Mode::from_index(mi as usize).unwrap_or(Mode::MD66);
            f.frame_type = match (q, sti) {
                (false, _) => RxFrameType::SidBad,
                (true, 0) => RxFrameType::SidFirst,
                (true, _) => RxFrameType::SidUpdate,
            };
        }
        FT_SPEECH_LOST => f.frame_type = RxFrameType::SpeechLost,
        FT_NO_DATA => {}
        _ => unreachable!("frame type {ft} is rejected before reading"),
    }
    f
}

fn write_frame(w: &mut BitWriter, f: &Frame, ft: u8) {
    match ft {
        0..=8 => {
            for &i in SORT[ft as usize] {
                w.put(f.bits[i as usize] as u32 & 1, 1);
            }
        }
        FT_SID => {
            // SID_FIRST: zeros, as the reference's MIME writer (`Write_serial`) sends.
            let first = f.frame_type == RxFrameType::SidFirst;
            for &b in &f.bits[..SID_BITS] {
                w.put(if first { 0 } else { b as u32 & 1 }, 1);
            }
            w.put(u32::from(!first), 1);
            w.put(f.mode as u32, 4);
        }
        _ => {}
    }
}

/// MSB-first bit reader; callers check `remaining()` first.
#[derive(Clone)]
struct BitReader<'a> {
    data: &'a [u8],
    pos: usize,
}

impl BitReader<'_> {
    fn remaining(&self) -> usize {
        self.data.len() * 8 - self.pos
    }

    /// The next `n` (≤ 8) bits as an unsigned number (panics past the end).
    fn get(&mut self, n: usize) -> u32 {
        let mut v = 0;
        for _ in 0..n {
            let bit = self.data[self.pos / 8] >> (7 - self.pos % 8) & 1;
            v = v << 1 | u32::from(bit);
            self.pos += 1;
        }
        v
    }

    /// Skips to the next octet boundary (padding bits are ignored).
    fn align(&mut self) {
        self.pos = self.pos.next_multiple_of(8);
    }
}

/// MSB-first bit writer; padding bits are zero.
#[derive(Default)]
struct BitWriter {
    out: Vec<u8>,
    bits: usize,
}

impl BitWriter {
    fn put(&mut self, v: u32, n: usize) {
        for k in (0..n).rev() {
            if self.bits.is_multiple_of(8) {
                self.out.push(0);
            }
            let last = self.out.len() - 1;
            self.out[last] |= ((v >> k & 1) as u8) << (7 - self.bits % 8);
            self.bits += 1;
        }
    }

    fn align(&mut self) {
        self.bits = self.bits.next_multiple_of(8);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testutil::Rng;
    use RxFrameType::*;

    pub(super) fn storage(file: &[u8]) -> Result<Vec<Frame>, Error> {
        read_storage(file)?.collect()
    }

    fn rtp(payload: &[u8], p: Packing) -> Result<(Option<Mode>, Vec<Frame>), Error> {
        let (cmr, frames) = read_rtp(payload, p)?;
        Ok((cmr, frames.collect()))
    }

    /// A frame of `frame_type`/`mode` with random bits in exactly the bits its format carries.
    pub(super) fn random_frame(g: &mut Rng, frame_type: RxFrameType, mode: Mode) -> Frame {
        let n = match frame_type {
            NoData | SpeechLost => 0,
            SidFirst | SidUpdate | SidBad => SID_BITS,
            _ => mode.bits(),
        };
        let mut f = Frame::no_data();
        f.frame_type = frame_type;
        if !matches!(frame_type, NoData | SpeechLost) {
            f.mode = mode;
        }
        for b in &mut f.bits[..n] {
            *b = (g.u64() & 1) as i16;
        }
        f
    }

    /// `f` as the writers send it: a SID_FIRST frame's comfort-noise bits zeroed.
    pub(super) fn as_written(f: &Frame) -> Frame {
        let mut f = f.clone();
        if f.frame_type == SidFirst {
            f.bits[..SID_BITS].fill(0);
        }
        f
    }

    /// Every packable frame: 9 speech modes × good/bad, SID first/update/bad × 9 mode
    /// indications, SPEECH_LOST and NO_DATA.
    pub(super) fn all_kinds(g: &mut Rng) -> Vec<Frame> {
        let mut v = vec![Frame::no_data(), random_frame(g, SpeechLost, Mode::MD66)];
        for m in Mode::ALL {
            for t in [SpeechGood, SpeechBad, SidFirst, SidUpdate, SidBad] {
                v.push(random_frame(g, t, m));
            }
        }
        v
    }

    #[test]
    fn sort_tables_are_permutations() {
        for (m, s) in SORT.iter().enumerate() {
            let mut v: Vec<i16> = s.to_vec();
            v.sort();
            assert!(
                v.iter().copied().eq(0..Mode::ALL[m].bits() as i16),
                "mode {m}"
            );
        }
        assert!(tables::SORT_SID.iter().copied().eq(0..SID_BITS as i16));
    }

    #[test]
    fn round_trip_every_mode_sid_lost_and_no_data() {
        let mut g = Rng::new(4867);
        for _ in 0..200 {
            let frames = all_kinds(&mut g);
            let sent: Vec<Frame> = frames.iter().map(as_written).collect();
            let file = write_storage(&frames).unwrap();
            assert_eq!(storage(&file).unwrap(), sent, "storage");
            for packing in [Packing::OctetAligned, Packing::BandwidthEfficient] {
                // Single frames, then packets of 1..=4 mixed frames.
                for f in &frames {
                    let p = write_rtp(Some(f.mode), std::slice::from_ref(f), packing).unwrap();
                    assert_eq!(
                        rtp(&p, packing).unwrap(),
                        (Some(f.mode), vec![as_written(f)])
                    );
                }
                for chunk in frames.chunks(1 + (g.u64() % 4) as usize) {
                    let p = write_rtp(None, chunk, packing).unwrap();
                    assert_eq!(
                        rtp(&p, packing).unwrap(),
                        (None, chunk.iter().map(as_written).collect()),
                        "{packing:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn sizes_match_opencore_and_rfc() {
        let mut g = Rng::new(1);
        for f in all_kinds(&mut g) {
            let (ft, _) = frame_header(&f).unwrap();
            let bits = frame_bits(ft).unwrap();
            let storage = write_storage(std::slice::from_ref(&f)).unwrap();
            let want = 1 + tables::PACKED_SIZE[ft as usize] as usize;
            assert_eq!(storage.len() - STORAGE_MAGIC.len(), want, "storage FT {ft}");
            assert_eq!(
                storage[STORAGE_MAGIC.len()] & !4,
                tables::TOC_BYTE[ft as usize] & !4,
                "ToC FT {ft}"
            );
            let oa = write_rtp(None, std::slice::from_ref(&f), Packing::OctetAligned).unwrap();
            assert_eq!(oa.len(), 2 + bits.div_ceil(8), "octet-aligned FT {ft}");
            let be =
                write_rtp(None, std::slice::from_ref(&f), Packing::BandwidthEfficient).unwrap();
            assert_eq!(
                be.len(),
                (4 + 6 + bits).div_ceil(8),
                "bandwidth-efficient FT {ft}"
            );
        }
    }

    /// Header layouts worked out by hand from RFC 4867 §4.3–4.4 and §5.3.
    #[test]
    fn header_layout() {
        let mut f = Frame::no_data();
        f.frame_type = SpeechGood;
        f.mode = Mode::MD2385;
        // Serial bit 0 is sorted payload bit 0 for 23.85 kbit/s (sort_2385[0] = 0).
        assert_eq!(tables::SORT_2385[0], 0);
        f.bits[0] = 1;
        let file = write_storage(std::slice::from_ref(&f)).unwrap();
        assert_eq!(file[..12], *b"#!AMR-WB\n\x44\x80\x00");
        // CMR 15, ToC F=0 FT=8 Q=1, then d(0)=1: 1111 0|100 0|1 + zeros.
        let be = write_rtp(None, std::slice::from_ref(&f), Packing::BandwidthEfficient).unwrap();
        assert_eq!(be[..3], [0xf4, 0x60, 0x00]);
        assert_eq!(be.len(), (4 + 6 + 477usize).div_ceil(8));
        // CMR 2 (12.65) | ToC F=0 FT=8 Q=1 | d(0)=1.
        let oa = write_rtp(
            Some(Mode::MD1265),
            std::slice::from_ref(&f),
            Packing::OctetAligned,
        )
        .unwrap();
        assert_eq!(oa[..3], [0x20, 0x44, 0x80]);
        assert_eq!(oa.len(), 2 + 60);
        // Two frames: the first ToC entry has F=1.
        let two = [f.clone(), Frame::no_data()];
        let oa = write_rtp(None, &two, Packing::OctetAligned).unwrap();
        assert_eq!(oa[..4], [0xf0, 0xc4, 0x7c, 0x80]);
        // SID update, MI = 18.25 (5 = 0101, MSB first), all CN bits 0.
        let mut sid = Frame::no_data();
        sid.frame_type = SidUpdate;
        sid.mode = Mode::MD1825;
        let file = write_storage(&[sid]).unwrap();
        // ToC FT 9 Q 1, 35 zero bits, STI=1, MI 0101: byte 4 = 000 1 0101.
        assert_eq!(file[9..], [0x4c, 0, 0, 0, 0, 0x15]);
        // SID_FIRST: the same with STI=0 and its comfort-noise bits written as zeros.
        let mut first = Frame::no_data();
        first.frame_type = SidFirst;
        first.mode = Mode::MD1825;
        first.bits[..SID_BITS].fill(1);
        assert_eq!(
            write_storage(&[first]).unwrap()[9..],
            [0x4c, 0, 0, 0, 0, 0x05]
        );
        // SPEECH_LOST: a ToC octet alone.
        let mut lost = Frame::no_data();
        lost.frame_type = SpeechLost;
        assert_eq!(write_storage(&[lost]).unwrap()[9..], [0x74]);
    }

    #[test]
    fn quality_bit_marks_bad_frames() {
        // Storage header for 6.60 with Q=0, then SID with Q=0 (MI 0101 in the last octet).
        let mut file = STORAGE_MAGIC.to_vec();
        file.push(0x00);
        file.extend([0; 17]);
        file.push(FT_SID << 3);
        file.extend([0xff, 0xff, 0xff, 0xff, 0xf5]);
        let frames = storage(&file).unwrap();
        assert_eq!(frames[0].frame_type, SpeechBad);
        assert_eq!(frames[1].frame_type, SidBad);
        assert_eq!(frames[1].mode, Mode::MD1825);
        assert_eq!(frames[1].bits[..SID_BITS], [1; SID_BITS]);
        assert_eq!(frames[1].bits[SID_BITS..], [0; MAX_SERIAL_SIZE - SID_BITS]);
    }

    /// A file reader keeps the frames before an error; a SID's mode indication that is not a
    /// mode reads as `MD66` (the decoder ignores it).
    #[test]
    fn storage_streams_past_bad_frames() {
        let mut g = Rng::new(3);
        let a = random_frame(&mut g, SpeechGood, Mode::MD66);
        let b = random_frame(&mut g, SidUpdate, Mode::MD1425);
        let body = |f: &Frame| write_storage(std::slice::from_ref(f)).unwrap()[9..].to_vec();
        for mi in 9..=15u8 {
            let mut bad = body(&b);
            bad[5] = bad[5] & 0xf0 | mi;
            let file = [STORAGE_MAGIC.as_slice(), &body(&a), &bad, &body(&b)].concat();
            let mut b0 = b.clone();
            b0.mode = Mode::MD66;
            assert_eq!(storage(&file), Ok(vec![a.clone(), b0, b.clone()]));
        }
        for ft in 10..=13u8 {
            let file = [
                STORAGE_MAGIC.as_slice(),
                &body(&a),
                &[ft << 3 | 4],
                &body(&b),
            ]
            .concat();
            let got: Vec<_> = read_storage(&file).unwrap().collect();
            assert_eq!(got, [Ok(a.clone()), Err(Error::UnsupportedFrameType(ft))]);
        }
        let file = [STORAGE_MAGIC.as_slice(), &body(&a), &body(&b)].concat();
        let got: Vec<_> = read_storage(&file[..file.len() - 1]).unwrap().collect();
        assert_eq!(got, [Ok(a), Err(Error::Truncated)]);
    }

    #[test]
    fn errors() {
        let mut g = Rng::new(2);
        let f = random_frame(&mut g, SpeechGood, Mode::MD1985);
        assert_eq!(storage(b"#!AMR-WB_MC1.0\n"), Err(Error::BadMagic));
        assert_eq!(storage(b"#!AMR\n"), Err(Error::BadMagic));
        assert_eq!(storage(b""), Err(Error::BadMagic));
        assert_eq!(storage(STORAGE_MAGIC), Ok(vec![]));
        let file = write_storage(std::slice::from_ref(&f)).unwrap();
        for cut in 10..file.len() {
            assert_eq!(storage(&file[..cut]), Err(Error::Truncated), "cut {cut}");
        }
        for ft in 10..=13u8 {
            for p in [Packing::OctetAligned, Packing::BandwidthEfficient] {
                // Unsupported FT in the second ToC entry discards the whole packet.
                let oa = [0xf0, 0x80 | 8 << 3 | 4, ft << 3 | 4, 0];
                let be = [0xf0 | 1 << 3 | 4, 0x40 | ft << 1 | 1, 0, 0];
                let pkt = if p == Packing::OctetAligned { oa } else { be };
                assert_eq!(rtp(&pkt, p), Err(Error::UnsupportedFrameType(ft)), "{p:?}");
            }
        }
        for p in [Packing::OctetAligned, Packing::BandwidthEfficient] {
            let pkt = write_rtp(None, std::slice::from_ref(&f), p).unwrap();
            for cut in 0..pkt.len() {
                assert_eq!(
                    rtp(&pkt[..cut], p),
                    Err(Error::Truncated),
                    "{p:?} cut {cut}"
                );
            }
            let long = [pkt.as_slice(), &[0]].concat();
            // Octets after the last frame are ignored.
            assert_eq!(rtp(&long, p), rtp(&pkt, p), "{p:?}");
            assert_eq!(write_rtp(None, &[], p), Err(Error::NotPackable));
            // A SID's mode indication 9..=15 reads as MD66; the packet is kept.
            let sid = random_frame(&mut g, SidUpdate, Mode::MD2385);
            let pkt = write_rtp(None, &[f.clone(), sid.clone()], p).unwrap();
            // MI 1000 → 1111. Octet-aligned: 8 + 16 header bits, 397 + 3 padded, SID's MI at
            // bits 460..464 (the last octet's low nibble). Bandwidth-efficient: 4 + 12 + 397,
            // MI at bits 449..453 (the last octet's bits 6..3).
            let mut bad = pkt.clone();
            *bad.last_mut().unwrap() |= if p == Packing::OctetAligned {
                0x0f
            } else {
                0x78
            };
            let mut sid0 = sid;
            sid0.mode = Mode::MD66;
            assert_eq!(rtp(&bad, p), Ok((None, vec![f.clone(), sid0])), "{p:?}");
        }
        // Invalid CMR values are ignored.
        assert_eq!(rtp(&[0x90, 0x7c], Packing::OctetAligned).unwrap().0, None);
        let mut f = Frame::no_data();
        f.frame_type = SpeechProbablyDegraded;
        assert_eq!(write_storage(&[f]), Err(Error::NotPackable));
    }
}

/// Against the TS 26.173 reference's storage reader (`Read_serial`, format 2 = MIME, as its
/// `decoder -mime` reads `.awb` files), run on the same bytes through a C `FILE`.
#[cfg(all(test, oracle))]
#[allow(unsafe_code)]
mod vs_ref {
    use super::tests::{all_kinds, as_written, storage};
    use super::*;
    use crate::testutil::{Rng, with_ref};
    use std::ffi::c_void;

    unsafe extern "C" {
        fn tmpfile() -> *mut c_void;
        fn fwrite(p: *const c_void, size: usize, n: usize, f: *mut c_void) -> usize;
        fn rewind(f: *mut c_void);
        fn fclose(f: *mut c_void) -> i32;
        fn refwb_Init_read_serial(st: *mut *mut c_void) -> i16;
        fn refwb_Close_read_serial(st: *mut c_void) -> i16;
        fn refwb_Read_serial(
            fp: *mut c_void,
            prms: *mut i16,
            frame_type: *mut i16,
            mode: *mut i16,
            st: *mut c_void,
            format: i16,
        ) -> i16;
    }

    /// One frame as the reference reads it: frame type, mode, bits (0/1, the frame's size).
    type RefFrame = (i16, i16, Vec<i16>);

    /// The frames of a storage body (no magic: `decoder.c` reads that itself).
    fn ref_read(body: &[u8]) -> Vec<RefFrame> {
        let mut out = Vec::new();
        with_ref(|| unsafe {
            let fp = tmpfile();
            assert!(!fp.is_null());
            assert_eq!(fwrite(body.as_ptr().cast(), 1, body.len(), fp), body.len());
            rewind(fp);
            let mut st = std::ptr::null_mut();
            refwb_Init_read_serial(&mut st);
            loop {
                let mut prms = [0i16; MAX_SERIAL_SIZE];
                let (mut ft, mut mode) = (0, 0);
                // SAFETY: prms holds the largest frame; the state and file are live.
                if refwb_Read_serial(fp, prms.as_mut_ptr(), &mut ft, &mut mode, st, 2) == 0 {
                    break;
                }
                let n = match ft {
                    4..=6 => SID_BITS,
                    2 | 7 => 0,
                    _ => Mode::ALL[mode as usize].bits(),
                };
                out.push((
                    ft,
                    mode,
                    prms[..n].iter().map(|&b| (b > 0) as i16).collect(),
                ));
            }
            refwb_Close_read_serial(st);
            fclose(fp);
        });
        out
    }

    /// Our frame in the reference's terms. The reference gives a SID frame the previous
    /// frame's mode (`prev_mode`, as opencore), ignoring its mode indication, as our decoder
    /// does (we report the indication in `Frame::mode`): `sid_mode`.
    fn ours(f: &Frame, sid_mode: i16) -> RefFrame {
        let (n, mode) = match f.frame_type {
            RxFrameType::SidFirst | RxFrameType::SidUpdate | RxFrameType::SidBad => {
                (SID_BITS, sid_mode)
            }
            RxFrameType::NoData | RxFrameType::SpeechLost => (0, sid_mode),
            _ => (f.mode.bits(), f.mode as i16),
        };
        (f.frame_type as i16, mode, f.bits[..n].to_vec())
    }

    /// Our frames vs the reference's; the mode the reference carries over follows its rule.
    fn compare(file: &[u8], frames: &[Frame]) {
        let theirs = ref_read(&file[STORAGE_MAGIC.len()..]);
        assert_eq!(theirs.len(), frames.len());
        let mut prev_mode = 0;
        for (i, (f, t)) in frames.iter().zip(&theirs).enumerate() {
            assert_eq!(ours(f, prev_mode), *t, "frame {i}");
            prev_mode = t.1;
        }
    }

    /// Our writer, read back by the reference: every packable frame kind.
    #[test]
    fn reference_reads_what_we_write() {
        let mut g = Rng::new(26_201);
        for _ in 0..100 {
            let frames = all_kinds(&mut g);
            let sent: Vec<Frame> = frames.iter().map(as_written).collect();
            compare(&write_storage(&frames).unwrap(), &sent);
        }
    }

    /// Random frame bodies (every frame type of known size, random Q, padding and mode
    /// indications, 9..=15 included) read by both readers.
    #[test]
    fn we_read_random_frames_like_the_reference() {
        let mut g = Rng::new(7);
        let mut body = Vec::new();
        let mut kinds = [0usize; 16];
        for _ in 0..20_000 {
            let ft = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9, 14, 15][(g.u64() % 13) as usize];
            let q = (g.u64() & 1) as u8;
            body.push(ft << 3 | q << 2 | (g.u64() & 0x83) as u8);
            let n = tables::PACKED_SIZE[ft as usize] as usize;
            body.extend((0..n).map(|_| g.u64() as u8));
            kinds[ft as usize] += 1;
        }
        assert!(kinds.iter().filter(|&&k| k > 0).count() == 12);
        let file = [STORAGE_MAGIC.as_slice(), &body].concat();
        compare(&file, &storage(&file).unwrap());
    }
}

#[cfg(all(test, oracle))]
#[allow(unsafe_code)]
mod vs_c {
    use super::tests::storage;
    use super::*;
    use crate::Decoder;
    use crate::testutil::Rng;
    use opencore_sys::{D_IF_decode, D_IF_exit, D_IF_init};

    /// RESULTS.md defect: opencore's `D_IF_decode` ignores the quality bit (`quality = 1`),
    /// so frames marked damaged (Q = 0) decode exactly like good ones instead of being
    /// concealed (as its AMR-NB interface).
    #[test]
    fn opencore_ignores_the_quality_bit() {
        let mut g = Rng::new(11);
        // SAFETY: D_IF_init allocates its own state.
        let (good, bad) = unsafe { (D_IF_init(), D_IF_init()) };
        let (mut dg, mut db) = (Decoder::new(), Decoder::new());
        let mut differs_in_rust = 0;
        for _ in 0..500 {
            let ft = (g.u64() % 9) as u8;
            let mut buf = [0u8; 64];
            buf.iter_mut().for_each(|b| *b = g.u64() as u8);
            let n = 1 + tables::PACKED_SIZE[ft as usize] as usize;
            let (mut a, mut b) = ([0i16; 320], [0i16; 320]);
            buf[0] = ft << 3 | 4;
            let fg = storage(&[STORAGE_MAGIC.as_slice(), &buf[..n]].concat()).unwrap();
            // SAFETY: `buf` (64 bytes) holds the largest frame (1 + 60), `a` 320 samples.
            unsafe { D_IF_decode(good, buf.as_ptr(), a.as_mut_ptr(), 0) };
            buf[0] = ft << 3;
            let fb = storage(&[STORAGE_MAGIC.as_slice(), &buf[..n]].concat()).unwrap();
            // SAFETY: as above.
            unsafe { D_IF_decode(bad, buf.as_ptr(), b.as_mut_ptr(), 0) };
            assert_eq!(a, b, "FT {ft}");
            let (mut ra, mut rb) = ([0i16; 320], [0i16; 320]);
            dg.decode(fg[0].frame_type, fg[0].mode, &fg[0].bits, &mut ra);
            db.decode(fb[0].frame_type, fb[0].mode, &fb[0].bits, &mut rb);
            differs_in_rust += usize::from(ra != rb);
        }
        // SAFETY: each state freed once.
        unsafe { (D_IF_exit(good), D_IF_exit(bad)) };
        assert!(
            differs_in_rust > 400,
            "Rust conceals Q = 0 frames ({differs_in_rust}/500 differ)"
        );
    }
}
