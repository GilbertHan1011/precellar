use std::ops::Range;

use noodles::fastq::{self, record::Definition};

use crate::{utils::hamming_distance, Region, RegionType, SequenceType};

#[derive(Debug, Clone)]
pub enum SegmentType<'a> {
    Single {
        region_id: &'a str,
        region_type: RegionType,
    },
    Joined(Vec<RegionType>),
}

impl core::fmt::Display for SegmentType<'_> {
    fn fmt(&self, f: &mut core::fmt::Formatter) -> core::fmt::Result {
        match self {
            SegmentType::Single { region_type, .. } => write!(f, "{}", region_type),
            SegmentType::Joined(region_types) => {
                let joined = region_types
                    .iter()
                    .map(|r| r.to_string())
                    .collect::<Vec<String>>()
                    .join("");
                write!(f, "{}", joined)
            }
        }
    }
}

/// A segment of a read, consisting of a region type, sequence, and quality scores.
#[derive(Debug, Clone)]
pub struct Segment<'a> {
    pub region_type: SegmentType<'a>,
    pub seq: &'a [u8],
    pub qual: &'a [u8],
}

impl core::fmt::Display for Segment<'_> {
    fn fmt(&self, f: &mut core::fmt::Formatter) -> core::fmt::Result {
        write!(
            f,
            "[{}]{}",
            self.region_type,
            std::str::from_utf8(self.seq).unwrap()
        )
    }
}

impl<'a> Segment<'a> {
    pub fn new(region_type: SegmentType<'a>, seq: &'a [u8], qual: &'a [u8]) -> Self {
        Self {
            region_type,
            seq,
            qual,
        }
    }

    pub fn region_id(&self) -> &str {
        match &self.region_type {
            SegmentType::Single { region_id, .. } => region_id,
            SegmentType::Joined(_) => panic!("Joined segments have no region id"),
        }
    }

    pub fn into_fq(&self, defi: &Definition) -> fastq::Record {
        fastq::Record::new(defi.clone(), self.seq, self.qual)
    }

    pub fn is_barcode(&self) -> bool {
        match &self.region_type {
            SegmentType::Single { region_type, .. } => region_type.is_barcode(),
            _ => false,
        }
    }

    pub fn is_umi(&self) -> bool {
        match &self.region_type {
            SegmentType::Single { region_type, .. } => region_type.is_umi(),
            _ => false,
        }
    }

    pub fn contains_target(&self) -> bool {
        match &self.region_type {
            SegmentType::Single { region_type, .. } => region_type.is_target(),
            SegmentType::Joined(region_types) => region_types.iter().any(|r| r.is_target()),
        }
    }
}

#[derive(Debug, Clone)]
pub struct SegmentInfo {
    elems: Vec<SegmentInfoElem>,
    is_reverse: bool,
}

#[derive(Debug, Clone)]
pub struct SplitResult<'a> {
    pub segments: Vec<Segment<'a>>,
    pub mismatches: usize,
    pub is_reverse: bool,
}

/// Result of trying both orientations when parsing an unstranded read.
/// Used for two-phase strand selection (e.g. with barcode validation).
#[derive(Debug)]
pub enum SplitAutoRcBothResult<'a> {
    /// Only one orientation parsed successfully.
    One(SplitResult<'a>),
    /// Both forward and reverse parsed; caller can choose by mismatch count or barcode.
    Both {
        forward: SplitResult<'a>,
        reverse: SplitResult<'a>,
    },
}

#[derive(Debug, Copy, Clone)]
pub enum SplitError {
    AnchorNotFound,
    PatternMismatch(usize),
}

impl std::fmt::Display for SplitError {
    fn fmt(&self, f: &mut core::fmt::Formatter) -> core::fmt::Result {
        match self {
            SplitError::AnchorNotFound => write!(f, "Anchor not found"),
            SplitError::PatternMismatch(d) => write!(f, "Pattern mismatch ({} mismatches)", d),
        }
    }
}

impl std::error::Error for SplitError {}

fn make_rc_record(read: &fastq::Record) -> fastq::Record {
    let rc_seq = crate::utils::rev_compl(read.sequence());
    let mut rc_qual = read.quality_scores().to_vec();
    rc_qual.reverse();
    fastq::Record::new(
        fastq::record::Definition::new("", ""),
        rc_seq,
        rc_qual,
    )
}

/// Internal: mapping from RC segment coordinates to original read (for building SplitResult).
struct SegmentMapping<'a> {
    region_id: &'a str,
    region_type_variant: RegionType,
    is_joined: bool,
    joined_types: Vec<RegionType>,
    orig_start: usize,
    orig_end: usize,
}

impl SegmentInfo {
    pub fn new<T, S>(iter: T, is_reverse: bool) -> Self
    where
        T: IntoIterator<Item = S>,
        S: Into<SegmentInfoElem>,
    {
        let elems = iter.into_iter().map(|x| x.into()).collect::<Vec<_>>();
        Self { elems, is_reverse }
    }

    pub fn len(&self) -> usize {
        self.elems.len()
    }

    pub fn is_reverse(&self) -> bool {
        self.is_reverse
    }

    pub fn iter<'a>(&'a self) -> impl Iterator<Item = &'a SegmentInfoElem> {
        self.elems.iter()
    }

    pub fn into_iter(self) -> impl Iterator<Item = SegmentInfoElem> {
        self.elems.into_iter()
    }

    /// Truncate the segment information to the given length L, such that L is larger
    /// than the sum of the maximum lengths of the n-1 segments + minimum length of the last segment.
    pub fn truncate_max(self, len: usize) -> Self {
        let is_reverse = self.is_reverse;
        let mut sum = 0;
        let mut result = Vec::new();
        let mut segments = self.into_iter();
        while let Some(elem) = segments.next() {
            if sum + elem.max_len() > len {
                if sum + elem.min_len() <= len {
                    result.push(elem);
                }
                break;
            } else {
                sum += elem.max_len();
                result.push(elem);
            }
        }
        Self {
            elems: result,
            is_reverse: is_reverse,
        }
    }

    /// Split a read into segments according to the segment information.
    /// This function uses a default tolerance of 1.0 for the linker sequence and 0.2 for the anchor sequence,
    /// which means that we do not check the sequence of the linker and allow up to 20% mismatches for the anchor sequence.
    pub fn split<'a>(&'a self, read: &'a fastq::Record) -> Result<Vec<Segment<'a>>, SplitError> {
        self.split_with_tolerance(read, 1.0, 0.2)
    }

    /// Split a read into segments according to the segment information.
    ///
    /// Backward-compatible wrapper over split_with_score.
    pub fn split_with_tolerance<'a>(
        &'a self,
        read: &'a fastq::Record,
        linker_tolerance: f64,
        anchor_tolerance: f64,
    ) -> Result<Vec<Segment<'a>>, SplitError> {
        self.split_with_score(read, linker_tolerance, anchor_tolerance)
            .map(|res| res.segments)
    }

    /// Split a read into segments according to the segment information.
    ///
    /// # Arguments
    ///
    /// * `read` - The read to split
    /// * `linker_tolerance` - The fraction of mismatches allowed for the linker sequence
    /// * `anchor_tolerance` - The fraction of mismatches allowed for the anchor sequence
    pub fn split_with_score<'a>(
        &'a self,
        read: &'a fastq::Record,
        linker_tolerance: f64,
        anchor_tolerance: f64,
    ) -> Result<SplitResult<'a>, SplitError> {
        fn consume_buf<'a>(
            buf: &mut Vec<&'a SegmentInfoElem>,
            mut seq: &'a [u8],
            mut qual: &'a [u8],
            result: &mut Vec<Segment<'a>>,
        ) {
            let n1 = result.len();
            while let Some(elem) = buf.pop() {
                if !elem.is_fixed_len() {
                    buf.push(elem);
                    break;
                }

                let l = seq.len();
                let i = l - elem.max_len();
                let s = &seq[i..l];
                let q = &qual[i..l];
                seq = &seq[..i];
                qual = &qual[..i];
                result.push(Segment::new(
                    SegmentType::Single {
                        region_id: elem.region_id.as_str(),
                        region_type: elem.region_type,
                    },
                    s,
                    q,
                ));
            }

            if !buf.is_empty() {
                let segment_type = if buf.len() == 1 {
                    SegmentType::Single {
                        region_id: buf[0].region_id.as_str(),
                        region_type: buf[0].region_type,
                    }
                } else {
                    SegmentType::Joined(buf.iter().map(|s| s.region_type).collect())
                };
                buf.clear();
                result.push(Segment::new(segment_type, seq, qual));
            }

            let n2 = result.len();
            if n2 > n1 {
                result[n1..n2].reverse();
            }
        }

        let mut result = Vec::new();
        let mut min_offset = 0;
        let mut max_offset = 0;
        let mut total_mismatches = 0usize;
        let mut buffer: Vec<&SegmentInfoElem> = Vec::new();
        let len = read.sequence().len();

        for segment in self.iter() {
            let min_len = segment.min_len();
            let max_len = segment.max_len();

            // Check if the segment is within the read bounds
            if max_offset >= len || min_offset + min_len > len {
                //if max_offset + min_len > len {
                break;
            }

            if min_len == max_len {
                if !buffer.is_empty() {
                    // if the sequcence contains fixed patterns, we can try to use them as anchors
                    // to find the location of the next segment
                    if segment.sequence_type.is_fixed() {
                        // by definition, the pattern can only be found within the window defined by [min_offset, max_offset+max_len]
                        let search_start = min_offset;
                        let search_end = len.min(max_offset + max_len);
                        let seq = read
                            .sequence()
                            .get(search_start..search_end)
                            .unwrap();

                        let match_result = if self.is_reverse {
                            find_best_pattern_match(seq, segment.rc_sequence.as_slice())
                        } else {
                            find_best_pattern_match(seq, segment.sequence.as_bytes())
                        };

                        if let (Some(pos), mis) = match_result
                        {
                            if mis as f64 <= anchor_tolerance * segment.sequence.len() as f64 {
                                total_mismatches += mis;
                                // [offset_left, offset_right] is the region of the read that
                                // contains segments prior to the matched pattern
                                let offset_left =
                                    min_offset - buffer.iter().map(|s| s.min_len()).sum::<usize>();
                                let offset_right = min_offset + pos;

                                consume_buf(
                                    &mut buffer,
                                    read.sequence().get(offset_left..offset_right).unwrap(),
                                    read.quality_scores()
                                        .get(offset_left..offset_right)
                                        .unwrap(),
                                    &mut result,
                                );

                                // as the lengths of the previous segments are identified, we can now
                                // update the min_offset and max_offset
                                min_offset = offset_right;
                                max_offset = offset_right;

                                let seq = read
                                    .sequence()
                                    .get(offset_right..offset_right + max_len)
                                    .unwrap();
                                let qual = read
                                    .quality_scores()
                                    .get(offset_right..offset_right + max_len)
                                    .unwrap();

                                result.push(Segment::new(
                                    SegmentType::Single {
                                        region_id: segment.region_id.as_str(),
                                        region_type: segment.region_type,
                                    },
                                    seq,
                                    qual,
                                ))
                            } else if max_offset + max_len > len {
                                // if the pattern is not found near read end, treat as structural failure
                                // instead of silently returning a partial split
                                return Err(SplitError::AnchorNotFound);
                            } else {
                                return Err(SplitError::PatternMismatch(mis));
                            }
                        } else {
                            return Err(SplitError::AnchorNotFound);
                        }
                    } else {
                        buffer.push(segment);
                    }
                } else {
                    let seq = read
                        .sequence()
                        .get(max_offset..max_offset + max_len)
                        .unwrap();
                    let qual = read
                        .quality_scores()
                        .get(max_offset..max_offset + max_len)
                        .unwrap();

                    if linker_tolerance < 1.0 && segment.sequence_type.is_fixed() {
                        let d = if self.is_reverse {
                            hamming_distance(seq, segment.rc_sequence.as_slice()).unwrap()
                        } else {
                            hamming_distance(seq, segment.sequence.as_bytes()).unwrap()
                        };

                        if d as f64 > linker_tolerance * seq.len() as f64 {
                            return Err(SplitError::PatternMismatch(d as usize));
                        }
                        total_mismatches += d as usize;
                    }

                    result.push(Segment::new(
                        SegmentType::Single {
                            region_id: segment.region_id.as_str(),
                            region_type: segment.region_type,
                        },
                        seq,
                        qual,
                    ))
                }
            } else {
                buffer.push(segment);
            }

            min_offset += min_len;
            max_offset += max_len;
        }

        if !buffer.is_empty() {
            let offset_left = min_offset - buffer.iter().map(|s| s.min_len()).sum::<usize>();
            consume_buf(
                &mut buffer,
                read.sequence()
                    .get(offset_left..max_offset.min(len))
                    .unwrap(),
                read.quality_scores()
                    .get(offset_left..max_offset.min(len))
                    .unwrap(),
                &mut result,
            );
        }

        Ok(SplitResult {
            segments: result,
            mismatches: total_mismatches,
            is_reverse: self.is_reverse,
        })
    }

    /// Automatically try both orientations using sequence RC and coordinate mapping.
    ///
    /// Delegates to `split_auto_rc_both` and selects one result: fewer mismatches wins; on tie, forward.
    pub fn split_auto_rc<'a>(
        fwd_info: &'a SegmentInfo,
        read: &'a fastq::Record,
        linker_tolerance: f64,
        anchor_tolerance: f64,
    ) -> Result<SplitResult<'a>, SplitError> {
        let both_res = Self::split_auto_rc_both(fwd_info, read, linker_tolerance, anchor_tolerance)?;
        match both_res {
            SplitAutoRcBothResult::One(res) => Ok(res),
            SplitAutoRcBothResult::Both { forward, reverse } => {
                if reverse.mismatches < forward.mismatches {
                    Ok(reverse)
                } else {
                    Ok(forward)
                }
            }
        }
    }

    /// Try both forward and reverse orientations; return both results when both parse.
    /// Used for two-phase strand selection (e.g. with barcode validation).
    pub fn split_auto_rc_both<'a>(
        fwd_info: &'a SegmentInfo,
        read: &'a fastq::Record,
        linker_tolerance: f64,
        anchor_tolerance: f64,
    ) -> Result<SplitAutoRcBothResult<'a>, SplitError> {
        let res_fwd = fwd_info.split_with_score(read, linker_tolerance, anchor_tolerance);
        let rc_read = make_rc_record(read);
        let res_rev = fwd_info.split_with_score(&rc_read, linker_tolerance, anchor_tolerance);
        let read_len = read.sequence().len();

        match (res_fwd, res_rev) {
            (Ok(fwd), Ok(rev)) => {
                let (mismatches, mappings) = Self::rev_to_mappings(&rev, &rc_read, read_len, fwd_info);
                let reverse = Self::build_reverse_mapped_result(read, fwd_info, mappings, mismatches);
                Ok(SplitAutoRcBothResult::Both { forward: fwd, reverse })
            }
            (Ok(fwd), Err(_)) => Ok(SplitAutoRcBothResult::One(fwd)),
            (Err(_), Ok(rev)) => {
                let (mismatches, mappings) = Self::rev_to_mappings(&rev, &rc_read, read_len, fwd_info);
                let reverse = Self::build_reverse_mapped_result(read, fwd_info, mappings, mismatches);
                Ok(SplitAutoRcBothResult::One(reverse))
            }
            (Err(e), Err(_)) => Err(e.clone()),
        }
    }

    /// Build (mismatches, Vec<SegmentMapping>) from reverse split result and RC read.
    /// Uses fwd_info to get region_id slices so mappings borrow from fwd_info, not rc_read.
    fn rev_to_mappings<'a>(
        rev: &SplitResult,
        rc_read: &fastq::Record,
        read_len: usize,
        fwd_info: &'a SegmentInfo,
    ) -> (usize, Vec<SegmentMapping<'a>>) {
        let rc_seq_ptr = rc_read.sequence().as_ptr() as usize;
        let mut mappings = Vec::with_capacity(rev.segments.len());
        for seg in &rev.segments {
            let seg_ptr = seg.seq.as_ptr() as usize;
            let start = seg_ptr - rc_seq_ptr;
            let end = start + seg.seq.len();
            let orig_start = read_len - end;
            let orig_end = read_len - start;
            match &seg.region_type {
                SegmentType::Single { region_id, region_type } => {
                    let region_id: &'a str = fwd_info
                        .iter()
                        .find(|e| e.region_id == *region_id)
                        .map(|e| e.region_id.as_str())
                        .expect("segment region_id should exist in fwd_info");
                    mappings.push(SegmentMapping {
                        region_id,
                        region_type_variant: *region_type,
                        is_joined: false,
                        joined_types: Vec::new(),
                        orig_start,
                        orig_end,
                    });
                }
                SegmentType::Joined(types) => {
                    mappings.push(SegmentMapping {
                        region_id: "",
                        region_type_variant: RegionType::Linker,
                        is_joined: true,
                        joined_types: types.clone(),
                        orig_start,
                        orig_end,
                    });
                }
            }
        }
        (rev.mismatches, mappings)
    }

    /// Build SplitResult (reverse orientation) from original read and mappings.
    fn build_reverse_mapped_result<'a>(
        read: &'a fastq::Record,
        fwd_info: &'a SegmentInfo,
        mappings: Vec<SegmentMapping<'a>>,
        mismatches: usize,
    ) -> SplitResult<'a> {
        let mut mapping_map: std::collections::HashMap<&'a str, Vec<SegmentMapping<'a>>> =
            std::collections::HashMap::new();
        let mut joined_mappings = Vec::new();
        for mapping in mappings {
            if mapping.is_joined {
                joined_mappings.push(mapping);
            } else if !mapping.region_id.is_empty() {
                mapping_map
                    .entry(mapping.region_id)
                    .or_insert_with(Vec::new)
                    .push(mapping);
            }
        }
        let mut mapped_segments = Vec::with_capacity(fwd_info.len() + joined_mappings.len());
        for elem in fwd_info.iter() {
            if let Some(mapping_vec) = mapping_map.get_mut(elem.region_id.as_str()) {
                if let Some(mapping) = mapping_vec.pop() {
                    let region_type = SegmentType::Single {
                        region_id: elem.region_id.as_str(),
                        region_type: mapping.region_type_variant,
                    };
                    let seq_slice = &read.sequence()[mapping.orig_start..mapping.orig_end];
                    let qual_slice = &read.quality_scores()[mapping.orig_start..mapping.orig_end];
                    mapped_segments.push(Segment::new(region_type, seq_slice, qual_slice));
                }
            }
        }
        for mapping in joined_mappings {
            let region_type = SegmentType::Joined(mapping.joined_types);
            let seq_slice = &read.sequence()[mapping.orig_start..mapping.orig_end];
            let qual_slice = &read.quality_scores()[mapping.orig_start..mapping.orig_end];
            mapped_segments.push(Segment::new(region_type, seq_slice, qual_slice));
        }
        SplitResult {
            segments: mapped_segments,
            mismatches,
            is_reverse: true,
        }
    }
}

#[derive(Debug, Clone)]
pub struct SegmentInfoElem {
    pub region_id: String,
    pub region_type: RegionType,
    pub sequence_type: SequenceType,
    pub sequence: String,
    /// Precomputed reverse complement of sequence; avoids allocs in split_with_score when is_reverse.
    pub rc_sequence: Vec<u8>,
    pub len: Range<usize>, // length of the segment
}

impl From<Region> for SegmentInfoElem {
    fn from(region: Region) -> Self {
        let rc_sequence = crate::utils::rev_compl(region.sequence.as_bytes());
        Self {
            region_id: region.region_id,
            region_type: region.region_type,
            sequence_type: region.sequence_type,
            sequence: region.sequence,
            rc_sequence,
            len: region.min_len as usize..region.max_len as usize,
        }
    }
}

impl From<&Region> for SegmentInfoElem {
    fn from(region: &Region) -> Self {
        region.clone().into()
    }
}

impl SegmentInfoElem {
    pub fn new(
        region_id: &str,
        region_type: RegionType,
        sequence_type: SequenceType,
        sequence: &str,
        min_len: usize,
        max_len: usize,
    ) -> Self {
        let sequence = sequence.to_string();
        let rc_sequence = crate::utils::rev_compl(sequence.as_bytes());
        Self {
            region_id: region_id.to_string(),
            region_type,
            sequence_type,
            sequence,
            rc_sequence,
            len: min_len..max_len,
        }
    }

    pub fn min_len(&self) -> usize {
        self.len.start
    }

    pub fn max_len(&self) -> usize {
        self.len.end
    }

    pub fn is_fixed_len(&self) -> bool {
        self.len.start == self.len.end
    }

    pub fn is_barcode(&self) -> bool {
        self.region_type.is_barcode()
    }

    pub fn is_umi(&self) -> bool {
        self.region_type.is_umi()
    }
}

/// Zero-allocation pattern search. For short DNA sequences, slice windows are faster than KMP.
/// I have benchmarked these two methods: this implementation is 2x faster than KMP.
fn find_pattern(haystack: &[u8], needle: &[u8]) -> Option<usize> {
    if needle.is_empty() || haystack.is_empty() || needle.len() > haystack.len() {
        return None;
    }
    haystack.windows(needle.len()).position(|w| w == needle)
}

/// Find the best match for a pattern within a haystack, always returning the position
/// with the minimal number of mismatches. Highly optimized for short patterns (< 10 characters).
///
/// # Arguments
///
/// * `haystack` - The sequence to search in
/// * `needle` - The pattern to search for
///
/// # Returns
///
/// A tuple of (Option<position>, mismatch_count) with the best match position and its mismatch count
fn find_best_pattern_match(haystack: &[u8], needle: &[u8]) -> (Option<usize>, usize) {
    if needle.is_empty() || haystack.is_empty() || needle.len() > haystack.len() {
        return (None, usize::MAX);
    }
    // Try exact matching first - it's very efficient
    if let Some(pos) = find_pattern(haystack, needle) {
        return (Some(pos), 0);
    }
    let n = haystack.len();
    let m = needle.len();

    // Ultra-optimized special case for 1-character patterns
    if m == 1 {
        // For single character patterns, we either have an exact match or nothing
        // We already checked for exact matches above
        return (None, 1);
    }
    let mut min_mismatches = usize::MAX;
    let mut best_pos = None;
    // Simple and efficient approach for very short patterns
    for i in 0..=n.saturating_sub(m) {
        let mut mismatches = 0;
        // Count mismatches at this position
        for j in 0..m {
            if haystack[i + j] != needle[j] {
                mismatches += 1;
                // Early termination if we exceed current best
                if mismatches >= min_mismatches && min_mismatches != usize::MAX {
                    break;
                }
            }
        }
        // Update if this is a better match
        if mismatches < min_mismatches {
            min_mismatches = mismatches;
            best_pos = Some(i);
        }
    }
    (best_pos, min_mismatches)
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;

    use super::*;

    fn new_fq(seq: &str) -> fastq::Record {
        fastq::Record::new(
            fastq::record::Definition::default(),
            seq.chars()
                .filter(|c| !c.is_whitespace())
                .collect::<String>()
                .as_bytes(),
            '@'.to_string().repeat(seq.len()).as_bytes(),
        )
    }

    fn split_seq_by<I: IntoIterator<Item = SegmentInfoElem>>(seq: &str, elems: I) -> String {
        let segment_info = SegmentInfo::new(elems, false);
        let fq = new_fq(seq);
        segment_info
            .split(&fq)
            .unwrap()
            .into_iter()
            .map(|s| s.to_string())
            .join(",")
    }

    fn bc(n: usize) -> SegmentInfoElem {
        SegmentInfoElem::new("bc", RegionType::Barcode, SequenceType::Random, "N", n, n)
    }

    fn umi(n: usize) -> SegmentInfoElem {
        SegmentInfoElem::new("umi", RegionType::Umi, SequenceType::Random, "N", n, n)
    }

    fn cdna(min: usize, max: usize) -> SegmentInfoElem {
        SegmentInfoElem::new(
            "cdna",
            RegionType::Cdna,
            SequenceType::Random,
            "N",
            min,
            max,
        )
    }

    fn poly(min: usize, max: usize) -> SegmentInfoElem {
        SegmentInfoElem::new(
            "polya",
            RegionType::PolyA,
            SequenceType::Random,
            "N",
            min,
            max,
        )
    }

    fn linker(seq: &str) -> SegmentInfoElem {
        SegmentInfoElem::new(
            "linker",
            RegionType::Linker,
            SequenceType::Fixed,
            seq,
            seq.len(),
            seq.len(),
        )
    }

    fn var(min: usize, max: usize) -> SegmentInfoElem {
        SegmentInfoElem::new(
            "var",
            RegionType::Linker,
            SequenceType::Random,
            "N",
            min,
            max,
        )
    }

    #[test]
    fn test_fixed() {
        assert_eq!(
            split_seq_by("AAAAAA", [cdna(1, 5), linker("ATGT")]),
            "[T]AAAAA",
        );

        assert_eq!(
            split_seq_by("AAAA TT CCC", [bc(4), umi(2), cdna(2, 4), poly(2, 4)]),
            "[B]AAAA,[U]TT,[T]CCC",
        );

        assert_eq!(
            split_seq_by("AATT", [bc(4), umi(2), cdna(2, 4), poly(2, 4)]),
            "[B]AATT",
        );

        assert_eq!(
            split_seq_by("AATT GG T", [bc(4), umi(2), cdna(2, 4), poly(2, 4)]),
            "[B]AATT,[U]GG",
        );

        assert_eq!(
            split_seq_by("AAAA TT CCG GG", [bc(4), umi(2), cdna(2, 4), poly(2, 4)]),
            "[B]AAAA,[U]TT,[TO]CCGGG",
        );

        assert_eq!(
            split_seq_by(
                "AAAA TT CCCCCCC GGGGGGG AAAA TT",
                [bc(4), umi(2), cdna(2, 10), poly(2, 10), bc(4), umi(2)]
            ),
            "[B]AAAA,[U]TT,[TO]CCCCCCCGGGGGGGAAAATT",
        );
    }

    #[test]
    fn test_variable() {
        assert_eq!(
            split_seq_by(
                "AAAA A TT ACTG CCC",
                [bc(4), var(1, 3), umi(2), linker("ACTG")]
            ),
            "[B]AAAA,[O]A,[U]TT,[O]ACTG",
        );

        assert_eq!(
            split_seq_by(
                "AAAA A TT ACCGG CCCACTGG",
                [bc(4), var(1, 3), umi(2), linker("ACTGG")]
            ),
            "[B]AAAA,[O]A,[U]TT,[O]ACCGG",
        );

        assert_eq!(
            split_seq_by(
                "AAAA A TT ACCGG CCCACTG",
                [bc(4), var(1, 3), umi(2), linker("ACTGG"), cdna(10, 100)]
            ),
            "[B]AAAA,[O]A,[U]TT,[O]ACCGG",
        );

        assert_eq!(
            split_seq_by(
                "AAAA A TT ACCGG CCCACTG",
                [bc(4), var(1, 3), umi(2), linker("ACTGG"), cdna(1, 100)]
            ),
            "[B]AAAA,[O]A,[U]TT,[O]ACCGG,[T]CCCACTG",
        );

        assert_eq!(
            split_seq_by(
                "AAAA A TT ACCGG CC CA GGGG TTTT AT GGGG ACTGGGGGACTG",
                [
                    bc(4),
                    var(1, 3),
                    umi(2),
                    linker("ACTGG"),
                    var(2, 4),
                    bc(2),
                    linker("GGGG"),
                    var(2, 4),
                    bc(2),
                    linker("GGGG"),
                    cdna(1, 100),
                ]
            ),
            "[B]AAAA,[O]A,[U]TT,[O]ACCGG,[O]CC,[B]CA,[O]GGGG,[O]TTTT,[B]AT,[O]GGGG,[T]ACTGGGGGACTG",
        );
    }

    #[test]
    fn test_truncate() {
        let info = [
            bc(4),
            var(2, 4),
            bc(2),
            linker("GGGG"),
            cdna(1, 1000),
        ];
        let info = SegmentInfo::new(info, false);
        println!("{:?}", info.truncate_max(50));
    }
}
