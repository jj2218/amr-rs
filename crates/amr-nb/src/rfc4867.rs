//! RFC 4867 packed AMR-NB: the `.amr` storage format (§5) and the RTP payload format (§4),
//! octet-aligned and bandwidth-efficient, converted to and from the decoder's serial bits.
//!
//! Supported: single channel; RTP without interleaving, CRCs or robust sorting (the SDP
//! defaults). Speech bits are carried in the class-sorted order of TS 26.101 (opencore's
//! `wmf_to_ets`/`ets_to_wmf` tables); a SID frame (frame type 8) carries 35 comfort-noise
//! bits, the SID type indicator (STI) and a 3-bit mode indication (least significant bit
//! first); a SID_FIRST frame's comfort-noise bits are written as zeros (TS 26.101 4.2.3;
//! decoders don't use them).
//!
//! ```
//! use amr_nb::{Decoder, rfc4867};
//! # let file = rfc4867::write_storage(&[rfc4867::Frame::no_data()]).unwrap();
//! let mut dec = Decoder::new();
//! let mut pcm = [0i16; 160];
//! for f in rfc4867::read_storage(&file)? {
//!     let f = f?; // frames before a bad one are still delivered
//!     dec.decode(f.frame_type, f.mode, &f.bits, &mut pcm);
//! }
//! # Ok::<(), rfc4867::Error>(())
//! ```

use crate::{MAX_SERIAL_SIZE, Mode, RxFrameType, SID_BITS, tables};

/// Magic number that starts a single-channel `.amr` file.
pub const STORAGE_MAGIC: &[u8; 6] = b"#!AMR\n";

const FT_SID: u8 = 8;
const FT_NO_DATA: u8 = 15;
/// CMR value for "no mode request".
const CMR_NONE: u32 = 15;

/// One frame in the decoder's terms: the arguments of [`crate::Decoder::decode`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Frame {
    /// Receive type. Reading yields `SpeechGood`/`SpeechBad`, `SidFirst`/`SidUpdate`/`SidBad`
    /// (`*Bad` when the frame's quality bit Q is 0) or `NoData`.
    pub frame_type: RxFrameType,
    /// Codec mode; for SID frames the mode indication; ignored (`MR475`) for `NoData`.
    pub mode: Mode,
    /// Serial bits, one per element (0 or 1), in the decoder's (TS 26.073) order. Bits past the
    /// frame's size are zero when read and ignored when written.
    pub bits: [i16; MAX_SERIAL_SIZE],
}

impl Frame {
    /// A `NoData` frame (nothing received).
    pub fn no_data() -> Self {
        Self {
            frame_type: RxFrameType::NoData,
            mode: Mode::MR475,
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
    /// The file doesn't start with `#!AMR\n` (multi-channel files are not supported).
    BadMagic,
    /// The data ends inside a header or frame.
    Truncated,
    /// Frame type 9..=14 (other codecs' SID, or reserved). RFC 4867 §4.3.2: discard the packet.
    /// In a file, frames of types 9..=11 (known size) are skipped and reading continues.
    UnsupportedFrameType(u8),
    /// Writing: an `Onset` or `SpeechDegraded` frame, or an RTP packet with no frames; none
    /// of these has a packed form.
    NotPackable,
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::BadMagic => f.write_str("not a single-channel AMR file (no #!AMR header)"),
            Error::Truncated => f.write_str("data ends inside a header or frame"),
            Error::UnsupportedFrameType(ft) => write!(f, "unsupported AMR frame type {ft}"),
            Error::NotPackable => f.write_str("frame has no RFC 4867 representation"),
        }
    }
}

impl std::error::Error for Error {}

/// Reads a single-channel `.amr` file (magic, then one header octet + padded frame per frame)
/// as a stream of frames. A frame of type 9..=11 yields `Err(UnsupportedFrameType)` and reading
/// continues; a truncated frame or a frame of unknown size (types 12..=14) ends the stream
/// with its error.
///
/// ```no_run
/// use amr_nb::{Decoder, MAX_SERIAL_SIZE, Mode, RxFrameType, rfc4867};
///
/// let file = std::fs::read("speech.amr")?;
/// let mut dec = Decoder::new();
/// let mut pcm = [0i16; 160];
/// for frame in rfc4867::read_storage(&file)? {
///     match frame {
///         Ok(f) => dec.decode(f.frame_type, f.mode, &f.bits, &mut pcm),
///         // Another codec's frame: conceal it like a lost one, keeping the timing.
///         Err(rfc4867::Error::UnsupportedFrameType(_)) => {
///             dec.decode(RxFrameType::NoData, Mode::MR475, &[0; MAX_SERIAL_SIZE], &mut pcm)
///         }
///         Err(e) => return Err(e.into()), // truncated: the stream ends
///     }
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
            Some(n) if (9..=11).contains(&ft) => {
                r.pos += n;
                Err(Error::UnsupportedFrameType(ft))
            }
            Some(_) => Ok(read_frame(&mut r, ft, h & 4 != 0)),
            None => Err(Error::UnsupportedFrameType(ft)),
        };
        r.align();
        done = matches!(res, Err(Error::Truncated)) || frame_bits(ft).is_none();
        Some(res)
    }))
}

/// Writes a single-channel `.amr` file.
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
        if (9..=14).contains(&ft) {
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
        need(&end, frame_bits(ft).unwrap_or(0))?;
        end.pos += frame_bits(ft).unwrap_or(0);
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
/// use amr_nb::{Encoder, Dtx, Mode, rfc4867::{self, Packing}};
///
/// // Three 20 ms frames per packet, asking the far end to send 7.4 kbit/s.
/// let mut enc = Encoder::new(Dtx::Off);
/// let frames: Vec<_> = (0..3).map(|_| enc.encode(Mode::MR122, &[0; 160])).collect();
/// let payload = rfc4867::write_rtp(Some(Mode::MR74), &frames, Packing::BandwidthEfficient)?;
///
/// let (cmr, received) = rfc4867::read_rtp(&payload, Packing::BandwidthEfficient)?;
/// assert_eq!(cmr, Some(Mode::MR74));
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
        NoData => Ok((FT_NO_DATA, true)),
        Onset | SpeechDegraded => Err(Error::NotPackable),
    }
}

/// Speech bit `k` of the sorted payload is serial bit `REORDER[ft][k]` (`wmf_to_ets`).
static REORDER: [&[i16]; 8] = [
    &tables::REORDERBITS_MR475,
    &tables::REORDERBITS_MR515,
    &tables::REORDERBITS_MR59,
    &tables::REORDERBITS_MR67,
    &tables::REORDERBITS_MR74,
    &tables::REORDERBITS_MR795,
    &tables::REORDERBITS_MR102,
    &tables::REORDERBITS_MR122,
];

/// Payload bits of frame type `ft` (`None`: reserved types 12..=14, size unknown).
fn frame_bits(ft: u8) -> Option<usize> {
    match ft {
        0..=11 => Some(tables::NUMOFBITS[ft as usize] as usize),
        FT_NO_DATA => Some(0),
        _ => None,
    }
}

/// Reads a frame of type 0..=8 or 15 whose `frame_bits` the caller has checked are present.
fn read_frame(r: &mut BitReader, ft: u8, q: bool) -> Frame {
    let mut f = Frame::no_data();
    match ft {
        0..=7 => {
            for &i in REORDER[ft as usize] {
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
            for b in &mut f.bits[..SID_BITS] {
                *b = r.get(1) as i16;
            }
            let sti = r.get(1);
            let mi = r.get(3);
            // Mode indication is sent least significant bit first (opencore `AMRDecode`).
            f.mode = Mode::ALL[(mi.reverse_bits() >> 29) as usize];
            f.frame_type = match (q, sti) {
                (false, _) => RxFrameType::SidBad,
                (true, 0) => RxFrameType::SidFirst,
                (true, _) => RxFrameType::SidUpdate,
            };
        }
        FT_NO_DATA => {}
        _ => unreachable!("frame type {ft} is rejected before reading"),
    }
    f
}

fn write_frame(w: &mut BitWriter, f: &Frame, ft: u8) {
    match ft {
        0..=7 => {
            for &i in REORDER[ft as usize] {
                w.put(f.bits[i as usize] as u32 & 1, 1);
            }
        }
        FT_SID => {
            // TS 26.101 4.2.3: a SID_FIRST frame's comfort-noise bits are set to 0.
            let first = f.frame_type == RxFrameType::SidFirst;
            for &b in &f.bits[..SID_BITS] {
                w.put(if first { 0 } else { b as u32 & 1 }, 1);
            }
            w.put(u32::from(f.frame_type != RxFrameType::SidFirst), 1);
            w.put((f.mode as u32).reverse_bits() >> 29, 3);
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
            NoData => 0,
            SidFirst | SidUpdate | SidBad => SID_BITS,
            _ => mode.bits(),
        };
        let mut f = Frame::no_data();
        if frame_type != NoData {
            f.frame_type = frame_type;
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

    /// Every packable frame: 8 speech modes × good/bad, SID first/update/bad × 8 mode
    /// indications, and NO_DATA.
    fn all_kinds(g: &mut Rng) -> Vec<Frame> {
        let mut v = vec![Frame::no_data()];
        for m in Mode::ALL {
            for t in [SpeechGood, SpeechBad, SidFirst, SidUpdate, SidBad] {
                v.push(random_frame(g, t, m));
            }
        }
        v
    }

    #[test]
    fn round_trip_every_mode_sid_and_no_data() {
        let mut g = Rng::new(4867);
        for _ in 0..200 {
            let frames = all_kinds(&mut g);
            let file = write_storage(&frames).unwrap();
            let sent: Vec<Frame> = frames.iter().map(as_written).collect();
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
        // TS 26.101 table 3: a SID frame is the comfort-noise bits, STI and mode indication.
        assert_eq!(
            tables::NUMOFBITS[FT_SID as usize] as usize,
            SID_BITS + 1 + 3
        );
        let speech = [95, 103, 118, 134, 148, 159, 204, 244]; // TS 26.101 table 2
        assert!(Mode::ALL.iter().map(|m| m.bits()).eq(speech));
        assert!(
            (0..9)
                .map(Mode::from_index)
                .eq(Mode::ALL.map(Some).into_iter().chain([None]))
        );
        let mut g = Rng::new(1);
        for f in all_kinds(&mut g) {
            let (ft, _) = frame_header(&f).unwrap();
            let bits = if ft == FT_NO_DATA {
                0
            } else {
                tables::NUMOFBITS[ft as usize] as usize
            };
            let storage = write_storage(std::slice::from_ref(&f)).unwrap();
            let want = 1 + tables::WMFDECBYTESPERFRAME[ft as usize] as usize;
            assert_eq!(storage.len() - STORAGE_MAGIC.len(), want, "storage FT {ft}");
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
        f.mode = Mode::MR122;
        // Serial bit 0 is sorted payload bit 0 for MR122 (reorderBits_MR122[0] = 0).
        f.bits[0] = 1;
        let file = write_storage(std::slice::from_ref(&f)).unwrap();
        assert_eq!(file[..9], *b"#!AMR\n\x3c\x80\x00");
        // CMR 15, ToC F=0 FT=7 Q=1, then d(0)=1: 1111 0|011 1|1 + zeros.
        let be = write_rtp(None, std::slice::from_ref(&f), Packing::BandwidthEfficient).unwrap();
        assert_eq!(be[..3], [0xf3, 0xe0, 0x00]);
        assert_eq!(be.len(), 32);
        // CMR 4 (MR74) | ToC F=0 FT=7 Q=1 | d(0)=1.
        let oa = write_rtp(
            Some(Mode::MR74),
            std::slice::from_ref(&f),
            Packing::OctetAligned,
        )
        .unwrap();
        assert_eq!(oa[..3], [0x40, 0x3c, 0x80]);
        assert_eq!(oa.len(), 33);
        // Two frames: the first ToC entry has F=1.
        let two = [f.clone(), Frame::no_data()];
        let oa = write_rtp(None, &two, Packing::OctetAligned).unwrap();
        assert_eq!(oa[..4], [0xf0, 0xbc, 0x7c, 0x80]);
        // SID update, MI = MR59 (2 = 010, sent LSB first as 0,1,0), all CN bits 0.
        let mut sid = Frame::no_data();
        sid.frame_type = SidUpdate;
        sid.mode = Mode::MR59;
        let file = write_storage(&[sid]).unwrap();
        // 35 zero bits, STI=1, MI 0 1 0, pad: byte 4 = 000 1 010 0.
        assert_eq!(file[6..], [0x44, 0, 0, 0, 0, 0x14]);
    }

    #[test]
    fn quality_bit_marks_bad_frames() {
        // Storage header for MR475 with Q=0, then SID with Q=0.
        let mut file = STORAGE_MAGIC.to_vec();
        file.push(0x00);
        file.extend([0; 12]);
        file.push(FT_SID << 3);
        file.extend([0xff; 5]);
        let frames = storage(&file).unwrap();
        assert_eq!(frames[0].frame_type, SpeechBad);
        assert_eq!(frames[1].frame_type, SidBad);
        assert_eq!(frames[1].mode, Mode::MR122);
        // Padding bits are ignored.
        assert_eq!(frames[1].bits[35..], [0; MAX_SERIAL_SIZE - 35]);
    }

    /// A file reader keeps the frames before an error, and skips frame types of known size.
    #[test]
    fn storage_streams_past_bad_frames() {
        let mut g = Rng::new(3);
        let a = random_frame(&mut g, SpeechGood, Mode::MR475);
        let b = random_frame(&mut g, SidUpdate, Mode::MR67);
        let body = |f: &Frame| write_storage(std::slice::from_ref(f)).unwrap()[6..].to_vec();
        for ft in 9..=11u8 {
            let n = 1 + tables::WMFDECBYTESPERFRAME[ft as usize] as usize;
            let mut other = vec![0xff; n];
            other[0] = ft << 3 | 4;
            let file = [STORAGE_MAGIC.as_slice(), &body(&a), &other, &body(&b)].concat();
            let got: Vec<_> = read_storage(&file).unwrap().collect();
            assert_eq!(
                got,
                [
                    Ok(a.clone()),
                    Err(Error::UnsupportedFrameType(ft)),
                    Ok(b.clone())
                ]
            );
        }
        let file = [STORAGE_MAGIC.as_slice(), &body(&a), &[0x6c], &body(&b)].concat();
        let got: Vec<_> = read_storage(&file).unwrap().collect();
        assert_eq!(got, [Ok(a.clone()), Err(Error::UnsupportedFrameType(13))]);
        let file = [STORAGE_MAGIC.as_slice(), &body(&a), &body(&b)].concat();
        let got: Vec<_> = read_storage(&file[..file.len() - 1]).unwrap().collect();
        assert_eq!(got, [Ok(a), Err(Error::Truncated)]);
    }

    #[test]
    fn errors() {
        let mut g = Rng::new(2);
        let f = random_frame(&mut g, SpeechGood, Mode::MR102);
        assert_eq!(storage(b"#!AMR_MC1.0\n"), Err(Error::BadMagic));
        assert_eq!(storage(b""), Err(Error::BadMagic));
        assert_eq!(storage(STORAGE_MAGIC), Ok(vec![]));
        let file = write_storage(std::slice::from_ref(&f)).unwrap();
        for cut in 7..file.len() {
            assert_eq!(storage(&file[..cut]), Err(Error::Truncated), "cut {cut}");
        }
        for ft in 9..=14u8 {
            let file = [STORAGE_MAGIC.as_slice(), &[ft << 3 | 4, 0, 0, 0, 0, 0, 0]].concat();
            assert_eq!(storage(&file), Err(Error::UnsupportedFrameType(ft)));
            for p in [Packing::OctetAligned, Packing::BandwidthEfficient] {
                // Unsupported FT in the second ToC entry discards the whole packet.
                let oa = [0xf0, 0x80 | 7 << 3 | 4, ft << 3 | 4, 0];
                let be = [0xf0 | 1 << 3 | 3, 0xc0 | ft << 1 | 1, 0, 0];
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
        }
        // Invalid CMR values are ignored.
        assert_eq!(rtp(&[0x80, 0x7c], Packing::OctetAligned).unwrap().0, None);
        for (t, m) in [(Onset, Mode::MR122), (SpeechDegraded, Mode::MR122)] {
            let mut f = Frame::no_data();
            f.frame_type = t;
            f.mode = m;
            assert_eq!(write_storage(&[f]), Err(Error::NotPackable), "{t:?} {m:?}");
        }
    }
}

#[cfg(all(test, oracle))]
#[allow(unsafe_code)]
mod vs_c {
    use super::*;
    use crate::Decoder;
    use crate::testutil::Rng;
    use opencore_sys::{
        Decoder_Interface_Decode, Decoder_Interface_exit, Decoder_Interface_init,
        Encoder_Interface_Encode, Encoder_Interface_exit, Encoder_Interface_init,
    };
    use std::ffi::c_void;
    use std::ptr::null_mut;

    unsafe extern "C" {
        fn AMREncodeInit(enc: *mut *mut c_void, sid: *mut *mut c_void, dtx: i32) -> i16;
        fn GSMInitDecode(state: *mut *mut c_void, id: *mut i8) -> i16;
        fn GSMDecodeFrameExit(state: *mut *mut c_void);
        fn AMRDecode(st: *mut c_void, ft: i32, bits: *mut u8, pcm: *mut i16, fmt: i32) -> i16;
        fn AMREncodeExit(enc: *mut *mut c_void, sid: *mut *mut c_void);
        fn AMREncode(
            enc: *mut c_void,
            sid: *mut c_void,
            mode: i32,
            input: *mut i16,
            out: *mut u8,
            frame_type: *mut i32,
            format: i16,
        ) -> i16;
    }
    const AMR_TX_ETS: i16 = 2;

    /// Speech-like bursts alternating with near-silence, so DTX produces SID and NO_DATA frames.
    fn signal(g: &mut Rng, frame: usize) -> [i16; 160] {
        let loud = (frame / 40).is_multiple_of(2);
        let f0 = 0.02 + 0.01 * ((frame / 7) % 5) as f64;
        std::array::from_fn(|i| {
            let t = (frame * 160 + i) as f64;
            let noise = (g.u64() % 64) as f64 - 32.0;
            if loud {
                (4000.0 * (f0 * t).sin() + 1500.0 * (3.1 * f0 * t).sin() + 10.0 * noise) as i16
            } else {
                (noise / 4.0) as i16
            }
        })
    }

    /// opencore's encoder (DTX on, modes cycling) run twice in lockstep: ETS serial output
    /// (as our `Frame`) and its packed `Encoder_Interface_Encode` output.
    fn opencore_stream(frames: usize) -> Vec<(Frame, Vec<u8>)> {
        let mut g = Rng::new(26_101);
        let (mut enc, mut sid) = (null_mut(), null_mut());
        let mut v = Vec::with_capacity(frames);
        unsafe {
            AMREncodeInit(&mut enc, &mut sid, 1);
            let packed_enc = Encoder_Interface_init(1);
            for n in 0..frames {
                let mode = ((n / 3) % 8) as i32;
                let pcm = signal(&mut g, n);
                let mut ets = [0i16; MAX_SERIAL_SIZE + 2];
                let (mut a, mut ft) = (pcm, 0);
                AMREncode(
                    enc,
                    sid,
                    mode,
                    a.as_mut_ptr(),
                    ets.as_mut_ptr().cast(),
                    &mut ft,
                    AMR_TX_ETS,
                );
                let mut packed = [0u8; 64];
                let b = pcm; // the C overwrites its "const" input
                let len =
                    Encoder_Interface_Encode(packed_enc, mode, b.as_ptr(), packed.as_mut_ptr(), 0);
                let mut f = Frame::no_data();
                f.frame_type = match ets[0] {
                    0 => RxFrameType::SpeechGood,
                    1 => RxFrameType::SidFirst,
                    2 => RxFrameType::SidUpdate,
                    3 => RxFrameType::NoData,
                    t => panic!("TX type {t}"),
                };
                if f.frame_type != RxFrameType::NoData {
                    f.mode = Mode::ALL[ets[MAX_SERIAL_SIZE + 1] as usize];
                    f.bits.copy_from_slice(&ets[1..=MAX_SERIAL_SIZE]);
                }
                v.push((f, packed[..len as usize].to_vec()));
            }
            Encoder_Interface_exit(packed_enc);
            AMREncodeExit(&mut enc, &mut sid);
        }
        let count = |t| v.iter().filter(|(f, _)| f.frame_type == t).count();
        for t in [
            RxFrameType::SpeechGood,
            RxFrameType::SidFirst,
            RxFrameType::SidUpdate,
            RxFrameType::NoData,
        ] {
            assert!(count(t) > 0, "stream has no {t:?} frames");
        }
        v
    }

    /// Our writer reproduces opencore's packed encoder output byte for byte, except that a
    /// SID_FIRST frame's comfort-noise bits are zeros (TS 26.101 4.2.3; opencore sends them).
    #[test]
    fn write_storage_matches_opencore_encoder() {
        for (i, (f, packed)) in opencore_stream(4000).iter().enumerate() {
            let ours = write_storage(std::slice::from_ref(f)).unwrap();
            let mut expected = packed.clone();
            if f.frame_type == RxFrameType::SidFirst {
                // After the header octet: 35 comfort-noise bits, most significant first.
                for k in 0..SID_BITS {
                    expected[1 + k / 8] &= !(0x80 >> (k % 8));
                }
            }
            assert_eq!(
                &ours[STORAGE_MAGIC.len()..],
                expected,
                "frame {i}: {:?} {:?}",
                f.frame_type,
                f.mode
            );
        }
    }

    /// Decodes a storage body with opencore's `Decoder_Interface_Decode` and with
    /// `read_storage` + the Rust decoder; returns the number of frames compared.
    fn compare_decoders(frames: &[Vec<u8>]) -> usize {
        let file = [STORAGE_MAGIC.as_slice(), &frames.concat()].concat();
        let ours = super::tests::storage(&file).unwrap();
        assert_eq!(ours.len(), frames.len());
        let mut dec = Decoder::new();
        let st = unsafe { Decoder_Interface_init() };
        for (i, (packed, f)) in frames.iter().zip(&ours).enumerate() {
            let mut buf = [0u8; 64]; // opencore reads a fixed size per frame type
            buf[..packed.len()].copy_from_slice(packed);
            let (mut c, mut r) = ([0i16; 160], [0i16; 160]);
            unsafe { Decoder_Interface_Decode(st, buf.as_ptr(), c.as_mut_ptr(), 0) };
            dec.decode(f.frame_type, f.mode, &f.bits, &mut r);
            assert_eq!(r, c, "frame {i}: {:?} {:?}", f.frame_type, f.mode);
        }
        unsafe { Decoder_Interface_exit(st) };
        frames.len()
    }

    #[test]
    fn read_storage_decodes_like_opencore_on_encoder_output() {
        let frames: Vec<Vec<u8>> = opencore_stream(4000).into_iter().map(|(_, p)| p).collect();
        compare_decoders(&frames);
    }

    /// The reader on random payloads, against opencore's own unpacking: opencore decodes each
    /// stream twice, from the packed bytes (`Decoder_Interface_Decode`) and from our parsed
    /// frames (`AMRDecode`, ETS format); equal PCM throughout means equal frames. (Our decoder is
    /// checked against the TS 26.073 reference separately, `decoder::vs_ref`: on random frames
    /// opencore's arithmetic deviates from it.)
    fn compare_readers(frames: &[Vec<u8>]) -> usize {
        let file = [STORAGE_MAGIC.as_slice(), &frames.concat()].concat();
        let ours = super::tests::storage(&file).unwrap();
        assert_eq!(ours.len(), frames.len());
        let st = unsafe { Decoder_Interface_init() };
        let mut ets_st = std::ptr::null_mut();
        unsafe { GSMInitDecode(&mut ets_st, std::ptr::null_mut()) };
        for (i, (packed, f)) in frames.iter().zip(&ours).enumerate() {
            let mut buf = [0u8; 64]; // opencore reads a fixed size per frame type
            buf[..packed.len()].copy_from_slice(packed);
            let (mut c, mut r) = ([0i16; 160], [0i16; 160]);
            unsafe { Decoder_Interface_Decode(st, buf.as_ptr(), c.as_mut_ptr(), 0) };
            let mut ets = [0i16; MAX_SERIAL_SIZE + 2];
            ets[0] = f.frame_type as i16;
            ets[1..=MAX_SERIAL_SIZE].copy_from_slice(&f.bits);
            ets[MAX_SERIAL_SIZE + 1] = f.mode as i16;
            unsafe { AMRDecode(ets_st, 0, ets.as_mut_ptr().cast(), r.as_mut_ptr(), 0) };
            assert_eq!(r, c, "frame {i}: {:?} {:?}", f.frame_type, f.mode);
        }
        unsafe {
            Decoder_Interface_exit(st);
            GSMDecodeFrameExit(&mut ets_st);
        }
        frames.len()
    }

    /// Random payloads: every speech mode, SID (random STI and mode indication) and NO_DATA,
    /// with random padding bits.
    #[test]
    fn read_storage_reads_like_opencore_on_random_frames() {
        let mut g = Rng::new(7);
        let frames: Vec<Vec<u8>> = (0..20_000)
            .map(|_| {
                let ft = [0, 1, 2, 3, 4, 5, 6, 7, 7, 8, 15][(g.u64() % 11) as usize];
                let n = tables::WMFDECBYTESPERFRAME[ft as usize] as usize;
                let mut p = vec![ft << 3 | 4];
                p.extend((0..n).map(|_| g.u64() as u8));
                p
            })
            .collect();
        compare_readers(&frames);
    }

    /// RESULTS.md defect: opencore ignores the quality bit, so frames marked damaged (Q = 0)
    /// decode exactly like good ones instead of being concealed.
    #[test]
    fn opencore_ignores_the_quality_bit() {
        let mut g = Rng::new(11);
        let (good, bad) = unsafe { (Decoder_Interface_init(), Decoder_Interface_init()) };
        let mut differs_in_rust = 0;
        let (mut dg, mut db) = (Decoder::new(), Decoder::new());
        for _ in 0..500 {
            let ft = (g.u64() % 9) as u8;
            let mut buf = [0u8; 64];
            buf.iter_mut().for_each(|b| *b = g.u64() as u8);
            let (mut a, mut b) = ([0i16; 160], [0i16; 160]);
            buf[0] = ft << 3 | 4;
            let fg = super::tests::storage(
                &[
                    STORAGE_MAGIC.as_slice(),
                    &buf[..1 + tables::WMFDECBYTESPERFRAME[ft as usize] as usize],
                ]
                .concat(),
            )
            .unwrap();
            unsafe { Decoder_Interface_Decode(good, buf.as_ptr(), a.as_mut_ptr(), 0) };
            buf[0] = ft << 3;
            let fb = super::tests::storage(
                &[
                    STORAGE_MAGIC.as_slice(),
                    &buf[..1 + tables::WMFDECBYTESPERFRAME[ft as usize] as usize],
                ]
                .concat(),
            )
            .unwrap();
            unsafe { Decoder_Interface_Decode(bad, buf.as_ptr(), b.as_mut_ptr(), 0) };
            assert_eq!(a, b, "FT {ft}");
            let (mut ra, mut rb) = ([0i16; 160], [0i16; 160]);
            dg.decode(fg[0].frame_type, fg[0].mode, &fg[0].bits, &mut ra);
            db.decode(fb[0].frame_type, fb[0].mode, &fb[0].bits, &mut rb);
            differs_in_rust += usize::from(ra != rb);
        }
        unsafe { (Decoder_Interface_exit(good), Decoder_Interface_exit(bad)) };
        assert!(
            differs_in_rust > 400,
            "Rust conceals Q = 0 frames ({differs_in_rust}/500 differ)"
        );
    }

    /// RESULTS.md defect: for frame types 9..=14 opencore decodes nothing and leaves the
    /// caller's PCM buffer untouched, with no way to tell (the function returns `void`).
    #[test]
    fn opencore_leaves_pcm_unwritten_for_unsupported_frame_types() {
        let st = unsafe { Decoder_Interface_init() };
        for ft in 9..=14u8 {
            let buf = [ft << 3 | 4; 64];
            let mut pcm = [0x5555i16; 160];
            unsafe { Decoder_Interface_Decode(st, buf.as_ptr(), pcm.as_mut_ptr(), 0) };
            assert_eq!(pcm, [0x5555; 160], "FT {ft}");
        }
        unsafe { Decoder_Interface_exit(st) };
    }
}
