use super::aligners::{Aligner, MultiMap, MultiMapR};

use crate::adapter::trim_poly_nucleotide;
use crate::barcode::{BarcodeCorrector, OligoFrequncy, Whitelist, filter_cellular_barcodes_ordmag_advanced};
use crate::qc::{AlignQC, Metrics};
use crate::utils::rev_compl_fastq_record;
use anyhow::{bail, Result};
use bstr::BString;
use indexmap::IndexMap;
use itertools::Itertools;
use kdam::{tqdm, BarExt};
use log::{debug, info};
use noodles::{bam, fastq};
use rayon::iter::ParallelIterator;
use rayon::slice::ParallelSlice;
use seqspec::{Assay, FastqReader, Modality, Read, RegionId, SegmentInfo, SegmentInfoElem};
use smallvec::SmallVec;
use std::collections::{HashMap, HashSet};

pub struct FastqProcessor {
    assay: Assay,
    current_modality: Option<Modality>,
    mito_dna: HashSet<String>,
    metrics: HashMap<Modality, Metrics>,
    align_qc: HashMap<Modality, AlignQC>,
    barcode_correct_prob: f64, // if the posterior probability of a correction
    // exceeds this threshold, the barcode will be corrected.
    // cellrange uses 0.975 for ATAC and 0.9 for multiome.
    mismatch_in_barcode: usize, // The number of mismatches allowed in barcode
    expected_cells: Option<usize>, // Expected number of cells
    use_advanced_barcode_filtering: bool, // Whether to use advanced barcode filtering
    barcode_filtering_quantile: f64, // Quantile for barcode filtering
    barcode_bootstrap_samples: usize, // Number of bootstrap samples for filtering
}

impl FastqProcessor {
    pub fn new(assay: Assay) -> Self {
        debug!("Creating new FastqProcessor");
        Self {
            assay,
            current_modality: None,
            metrics: HashMap::new(),
            align_qc: HashMap::new(),
            mito_dna: HashSet::new(),
            barcode_correct_prob: 0.975,
            mismatch_in_barcode: 1,
            expected_cells: None,
            use_advanced_barcode_filtering: true, // Default to true for backward compatibility
            barcode_filtering_quantile: 0.99,
            barcode_bootstrap_samples: 100,
        }
    }

    pub fn with_barcode_correct_prob(mut self, prob: f64) -> Self {
        debug!("Setting barcode correction probability to {}", prob);
        self.barcode_correct_prob = prob;
        self
    }

    pub fn with_expected_cells(mut self, cells: usize) -> Self {
        debug!("Setting expected number of cells to {}", cells);
        self.expected_cells = Some(cells);
        self
    }

    pub fn with_advanced_barcode_filtering(mut self, enable: bool) -> Self {
        debug!("Setting advanced barcode filtering to {}", enable);
        self.use_advanced_barcode_filtering = enable;
        self
    }
    
    pub fn with_barcode_filtering_params(mut self, quantile: f64, bootstrap_samples: usize) -> Self {
        debug!("Setting barcode filtering parameters: quantile={}, bootstrap_samples={}", 
               quantile, bootstrap_samples);
        self.barcode_filtering_quantile = quantile;
        self.barcode_bootstrap_samples = bootstrap_samples;
        self
    }

    pub fn modality(&self) -> Modality {
        self.current_modality
            .expect("modality not set, please call set_modality first")
    }

    pub fn add_mito_dna(&mut self, mito_dna: impl Into<String>) {
        self.mito_dna.insert(mito_dna.into());
    }

    pub fn with_modality(mut self, modality: Modality) -> Self {
        debug!("Setting modality to {:?}", modality);
        self.current_modality = Some(modality);
        self
    }

    pub fn get_report(&self) -> Metrics {
        debug!("Generating metrics report for modality {:?}", self.modality());
        let mut metrics = self
            .metrics
            .get(&self.modality())
            .map_or(Metrics::default(), |x| x.clone());
        if let Some(align_qc) = self.align_qc.get(&self.modality()) {
            debug!("Adding alignment QC metrics to report");
            align_qc.report(&mut metrics);
        }
        metrics
    }

    /// Align reads and return the alignments.
    /// If the fastq file is paired-end, the alignments will be returned as a tuple.
    /// Otherwise, the alignments will be returned as a single vector.
    ///
    /// # Arguments
    ///
    /// * `num_threads` - The number of threads to use for alignment.
    /// * `chunk_size` - The maximum number of bases in a chunk.
    ///
    /// # Returns
    ///
    /// An iterator of alignments. If the fastq file is paired-end, the alignments will be returned as a tuple.
    /// Otherwise, the alignments will be returned as a single vector.
    pub fn gen_barcoded_alignments<'a, A: Aligner>(
        &'a mut self,
        aligner: &'a mut A,
        num_threads: u16,
        chunk_size: usize,
    ) -> impl Iterator<Item = Vec<(Option<MultiMapR>, Option<MultiMapR>)>> + 'a {
        info!("Starting gen_barcoded_alignments with chunk_size={}", chunk_size);
        let fq_reader = self.gen_barcoded_fastq(true).with_chunk_size(chunk_size);
        
        // Log reader state
        info!("FastqReader created. Is paired-end: {}", fq_reader.is_paired_end());
        info!("Number of annotators: {}", fq_reader.annotators.len());
        info!("Number of readers: {}", fq_reader.readers.len());
        
        // Log barcode and UMI information
        let barcodes = fq_reader.get_all_barcodes();
        info!("Barcodes found: {:?}", barcodes);
        let umis = fq_reader.get_all_umi();
        info!("UMIs found: {:?}", umis);

        let total_reads = fq_reader.total_reads.unwrap_or(0);
        info!("Total reads reported: {}", total_reads);

        let modality = self.modality();
        info!("Aligning {} reads...", total_reads);
        let header = aligner.header();
        let mut qc = AlignQC::default();
        self.mito_dna.iter().for_each(|mito| {
            header
                .reference_sequences()
                .get_index_of(&BString::from(mito.as_str()))
                .map(|x| qc.mito_dna.insert(x));
        });
        self.align_qc.insert(modality, qc);

        let mut progress_bar = tqdm!(total = total_reads);
        fq_reader.map(move |data| {
            debug!("Processing chunk with {} records", data.len());
            // Log details of first few records in chunk
            if !data.is_empty() {
                let sample = &data[0];
                debug!("Sample record - Barcode present: {}, UMI present: {}, Read1 present: {}, Read2 present: {}", 
                    sample.barcode.is_some(),
                    sample.umi.is_some(),
                    sample.read1.is_some(),
                    sample.read2.is_some()
                );
                
                // Check if at least one read is present
                if sample.read1.is_none() && sample.read2.is_none() {
                    panic!("Neither Read1 nor Read2 is present. Please provide at least one read.");
                }
            }
            
            let align_qc = self.align_qc.get_mut(&modality).unwrap();
            let results: Vec<_> = aligner.align_reads(num_threads, data);
            results.iter().for_each(|ali| match ali {
                (Some(ali1), Some(ali2)) => {
                    align_qc.add_pair(&header, ali1, ali2).unwrap();
                }
                (Some(ali1), None) => {
                    align_qc.add_read1(&header, ali1).unwrap();
                }
                (None, Some(ali2)) => {
                    align_qc.add_read2(&header, ali2).unwrap();
                }
                _ => {
                    debug!("No alignment found for read");
                }
            });
            progress_bar.update(results.len()).unwrap();
            results
        })
    }

    pub fn gen_barcoded_fastq(&mut self, correct_barcode: bool) -> AnnotatedFastqReader {
        let modality = self.modality();
        info!("Starting gen_barcoded_fastq for modality: {:?}", modality);

        let whitelists = if correct_barcode {
            info!("Counting barcodes...");
            debug!("Calling count_barcodes()");
            match self.count_barcodes() {
                Ok(wl) => {
                    debug!("Successfully counted barcodes. Found {} whitelists", wl.len());
                    // Log details for each whitelist
                    for (id, whitelist) in wl.iter() {
                        debug!(
                            "Whitelist '{}': exact match rate={:.2}%, unique barcodes={}",
                            id,
                            whitelist.frac_exact_match() * 100.0,
                            whitelist.num_seen_barcodes()
                        );
                    }
                    wl
                }
                Err(e) => {
                    debug!("Error counting barcodes: {}", e);
                    IndexMap::new()
                }
            }
        } else {
            debug!("Skipping barcode correction");
            IndexMap::new()
        };

        debug!("Creating BarcodeCorrector with threshold={}", self.barcode_correct_prob);
        let corrector = BarcodeCorrector::default()
            .with_max_missmatch(self.mismatch_in_barcode)
            .with_bc_confidence_threshold(self.barcode_correct_prob);

        debug!("Creating FastqAnnotators from segments");
        let mut fq_reader: AnnotatedFastqReader = self
            .assay
            .get_segments_by_modality(modality)
            .filter(|(read, _)| read.open().is_some())
            .filter_map(|(read, index)| {
                info!(
                    "Processing read {} with {} segments, is_reverse: {}", 
                    read.read_id,
                    index.segments.len(),
                    read.is_reverse()
                );
                
                // Log each segment's type
                for seg in &index.segments {
                    info!(
                        "Segment in read {}: type={:?}, range={:?}, id={}", 
                        read.read_id, 
                        seg.region_type, 
                        seg.range,
                        seg.region_id
                    );
                }

                let annotator = FastqAnnotator::new(read, index, &whitelists, corrector.clone())?;
                info!(
                    "Created annotator for read {} with {} subregions",
                    annotator.id,
                    annotator.subregions.len()
                );
                Some((annotator, read.open().unwrap()))
            })
            .collect();

        if !whitelists.is_empty() {
            fq_reader.total_reads = Some(whitelists[0].total_count);
            debug!("Set total_reads to {}", whitelists[0].total_count);
        }

        debug!("Completed gen_barcoded_fastq setup");
        fq_reader
    }

    fn count_barcodes(&mut self) -> Result<IndexMap<RegionId, Whitelist>> {
        let modality = self.modality();
        let mut whitelists = self.get_whitelists();
        let mut total_filtered_bcs = 0;

        fn count(
            read: &Read,
            barcode_region_index: SegmentInfo,
            whitelist: &mut Whitelist,
        ) -> Result<()> {
            let range = &barcode_region_index
                .segments
                .iter()
                .find(|x| x.region_type.is_barcode())
                .unwrap()
                .range;
            read.open().unwrap().records().for_each(|record| {
                let mut record = record.unwrap();
                record = slice_fastq_record(&record, range.start as usize, range.end as usize);
                if read.is_reverse() {
                    record = rev_compl_fastq_record(record);
                }
                whitelist.count_barcode(record.sequence(), record.quality_scores());
            });
            Ok(())
        }

        for (i, (read, barcode_region_index)) in self
            .assay
            .get_segments_by_modality(modality)
            .filter(|(_, region_index)| {
                region_index
                    .segments
                    .iter()
                    .any(|x| x.region_type.is_barcode())
            })
            .enumerate()
        {
            count(read, barcode_region_index, &mut whitelists[i])?;
        }

        // Apply advanced filtering to each whitelist if enabled and no predefined whitelist exists
        if self.use_advanced_barcode_filtering {
            for (id, whitelist) in whitelists.iter_mut() {
                // Use get_barcode_counts() to check if we have a predefined whitelist
                // We consider a whitelist "predefined" if it has entries but no counts
                let has_predefined_whitelist = !whitelist.get_barcode_counts().is_empty() && 
                                              whitelist.num_seen_barcodes() == 0;
                
                if !has_predefined_whitelist {
                    info!("No predefined whitelist for '{}', applying order-of-magnitude filtering", id);
                    
                    let results = filter_cellular_barcodes_ordmag_advanced(
                        whitelist,
                        self.expected_cells,
                        None,
                        None,
                        Some(self.barcode_filtering_quantile),
                        Some(self.barcode_bootstrap_samples),
                    );
                    
                    info!(
                        "Found {} cellular barcodes (95% CI: {}-{}) for '{}' with counts >= {}",
                        results.filtered_bcs, 
                        results.filtered_bcs_lb,
                        results.filtered_bcs_ub,
                        id,
                        results.filtered_bcs_cutoff
                    );
                    
                    total_filtered_bcs += results.filtered_bcs;
                    
                    // Store filtering metrics
                    self.metrics.entry(modality).or_default().insert(
                        format!("filtered_bcs_{}", id),
                        results.filtered_bcs as f64,
                    );
                    self.metrics.entry(modality).or_default().insert(
                        format!("filtered_bcs_cutoff_{}", id),
                        results.filtered_bcs_cutoff as f64,
                    );
                }
            }
            
            if total_filtered_bcs > 0 {
                self.metrics.entry(modality).or_default().insert(
                    "total_filtered_bcs".to_string(),
                    total_filtered_bcs as f64,
                );
            }
        } else {
            debug!("Advanced barcode filtering is disabled");
        }

        self.metrics.entry(modality).or_default().insert(
            "frac_q30_bases_barcode".to_string(),
            whitelists.values().map(|x| x.frac_q30_bases()).sum::<f64>() / whitelists.len() as f64,
        );
        
        Ok(whitelists)
    }

    fn get_whitelists(&self) -> IndexMap<RegionId, Whitelist> {
        debug!("Getting whitelists for modality {:?}", self.modality());
        let regions = self
            .assay
            .library_spec
            .get_modality(&self.modality())
            .unwrap()
            .read()
            .unwrap();
        
        let whitelists: IndexMap<RegionId, Whitelist> = regions
            .subregions
            .iter()
            .filter_map(|r| {
                let r = r.read().unwrap();
                if r.region_type.is_barcode() {
                    let id = r.region_id.to_string();
                    let list = if let Some(onlist) = r.onlist.as_ref() {
                        debug!("Found whitelist for barcode region {}", id);
                        Whitelist::new(onlist.read().unwrap())
                    } else {
                        debug!("No whitelist found for barcode region {}, using empty list", id);
                        Whitelist::empty()
                    };
                    Some((id, list))
                } else {
                    None
                }
            })
            .collect();

        debug!("Found {} whitelists", whitelists.len());
        whitelists
    }
}

pub struct AnnotatedFastqReader {
    buffer: fastq::Record,
    total_reads: Option<usize>,
    trim_poly_a: bool,
    annotators: Vec<FastqAnnotator>,
    readers: Vec<FastqReader>,
    chunk_size: usize,
    chunk: Vec<SmallVec<[fastq::Record; 4]>>,
}

impl AnnotatedFastqReader {
    pub fn with_chunk_size(mut self, chunk_size: usize) -> Self {
        self.chunk_size = chunk_size;
        self
    }

    pub fn with_polya_trimmed(mut self) -> Self {
        self.trim_poly_a = true;
        self
    }

    pub fn get_all_barcodes(&self) -> Vec<(&str, usize)> {
        self.annotators
            .iter()
            .flat_map(|annotator| {
                annotator
                    .subregions
                    .iter()
                    .filter(|info| info.region_type.is_barcode())
                    .map(|info| (info.region_id.as_str(), info.range.len()))
            })
            .collect()
    }

    pub fn get_all_umi(&self) -> Vec<(&str, usize)> {
        self.annotators
            .iter()
            .flat_map(|annotator| {
                annotator
                    .subregions
                    .iter()
                    .filter(|info| info.region_type.is_umi())
                    .map(|info| (info.region_id.as_str(), info.range.len()))
            })
            .collect()
    }

    pub fn is_paired_end(&self) -> bool {
        let mut has_read1 = false;
        let mut has_read2 = false;
        self.annotators.iter().for_each(|x| {
            x.subregions.iter().for_each(|info| {
                if info.region_type.is_target() {
                    if x.is_reverse {
                        has_read1 = true;
                    } else {
                        has_read2 = true;
                    }
                }
            });
        });
        has_read1 && has_read2
    }

    /// Read a chunk of records from the fastq files.
    fn read_chunk(&mut self) -> usize {
        self.chunk.clear();

        let mut accumulated_length = 0;
        let mut _total_records = 0;  // Prefix with underscore to show it's intentionally unused

        while accumulated_length < self.chunk_size {
            let mut max_read = 0;
            let mut min_read = usize::MAX;
            let records: SmallVec<[_; 4]> = self
                .readers
                .iter_mut()
                .enumerate()
                .flat_map(|(_i, reader)| {  // Prefix with underscore to show it's intentionally unused
                    let n = reader
                        .read_record(&mut self.buffer)
                        .expect("error reading fastq record");
                    min_read = min_read.min(n);
                    max_read = max_read.max(n);
                    if n > 0 {
                        accumulated_length += self.buffer.sequence().len();
                        Some(self.buffer.clone())
                    } else {
                        None
                    }
                })
                .collect();
            
            if max_read == 0 {
                debug!("No more records to read");
                break;
            } else if min_read == 0 {
                panic!("Unequal number of reads in the chunk");
            } else {
                _total_records += 1;
                assert!(
                    records.iter().map(|r| r.name()).all_equal(),
                    "read names mismatch"
                );
                self.chunk.push(records);
            }
        }
        
        self.chunk.len()
    }
}

impl FromIterator<(FastqAnnotator, FastqReader)> for AnnotatedFastqReader {
    fn from_iter<T: IntoIterator<Item = (FastqAnnotator, FastqReader)>>(iter: T) -> Self {
        let (annotators, readers): (Vec<_>, Vec<_>) = iter.into_iter().unzip();
        let chunk = Vec::new();
        Self {
            buffer: fastq::Record::default(),
            total_reads: None,
            annotators,
            readers,
            trim_poly_a: false,
            chunk_size: 10000000,
            chunk,
        }
    }
}

impl Iterator for AnnotatedFastqReader {
    type Item = Vec<AnnotatedFastq>;

    fn next(&mut self) -> Option<Self::Item> {
        let n = self.read_chunk();
        if n == 0 {
            None
        } else {
            let n = (n / 256).max(256);
            let annotators = &self.annotators;
            let result = self
                .chunk
                .par_chunks(n)
                .flat_map_iter(|chunk| {
                    chunk.into_iter().map(move |records| {
                        records
                            .iter()
                            .enumerate()
                            .map(|(i, record)| annotators[i].annotate(record).unwrap())
                            .reduce(|mut this, other| {
                                this.join(other);
                                this
                            })
                            .unwrap()
                    })
                })
                .collect();
            Some(result)
        }
    }
}

/// A FastqAnnotator that splits the reads into subregions, e.g., barcode, UMI, and
/// return annotated reads.
#[derive(Debug)]
struct FastqAnnotator {
    whitelists: IndexMap<String, OligoFrequncy>,
    corrector: BarcodeCorrector,
    id: String,
    is_reverse: bool,
    subregions: Vec<SegmentInfoElem>,
    min_len: usize,
    max_len: usize,
}

impl FastqAnnotator {
    pub fn new(
        read: &Read,
        index: SegmentInfo,
        whitelists: &IndexMap<String, Whitelist>,
        corrector: BarcodeCorrector,
    ) -> Option<Self> {
        let subregions: Vec<_> = index
            .segments
            .into_iter()
            .filter(|x| {
                x.region_type.is_barcode() || x.region_type.is_umi() || x.region_type.is_target()
            }) // only barcode and target regions
            .collect();
        if subregions.is_empty() {
            None
        } else {
            let whitelists = subregions
                .iter()
                .flat_map(|info| {
                    let v = whitelists.get(&info.region_id)?;
                    Some((info.region_id.clone(), v.get_barcode_counts().clone()))
                })
                .collect();
            let anno = Self {
                whitelists,
                corrector,
                id: read.read_id.clone(),
                is_reverse: read.is_reverse(),
                subregions,
                min_len: read.min_len as usize,
                max_len: read.max_len as usize,
            };
            Some(anno)
        }
    }

    fn annotate(&self, record: &fastq::Record) -> Result<AnnotatedFastq> {
        let n = record.sequence().len();
        if n < self.min_len || n > self.max_len {
            bail!(
                "Read length ({}) out of range: {}-{}",
                n,
                self.min_len,
                self.max_len
            );
        }

        let mut barcode: Option<Barcode> = None;
        let mut umi = None;
        let mut read1 = None;
        let mut read2 = None;

        self.subregions.iter().for_each(|info| {
            let mut fq =
                slice_fastq_record(record, info.range.start as usize, info.range.end as usize);
            if self.is_reverse && (info.region_type.is_barcode() || info.region_type.is_umi()) {
                fq = rev_compl_fastq_record(fq);
            }
            if info.region_type.is_umi() {
                umi = Some(fq);
            } else if info.region_type.is_barcode() {
                let corrected = self.whitelists.get(&info.region_id).map_or(
                    Some(fq.sequence().to_vec()),
                    |counts| {
                        self.corrector
                            .correct(counts, fq.sequence(), fq.quality_scores())
                            .ok()
                            .map(|x| x.to_vec())
                    },
                );
                if let Some(bc) = &mut barcode {
                    bc.extend(&Barcode { raw: fq, corrected });
                } else {
                    barcode = Some(Barcode { raw: fq, corrected });
                }
            } else if info.region_type.is_target() {
                if read1.is_some() || read2.is_some() {
                    panic!("Both Read1 and Read2 are set");
                } else {
                    if let Some(nucl) = info.region_type.poly_nucl() {
                        if let Some(idx) = trim_poly_nucleotide(nucl, fq.sequence().iter().copied())
                        {
                            fq = slice_fastq_record(&fq, idx, fq.sequence().len());
                        }
                    }
                    // Only keep reads with length >= 8
                    if fq.sequence().len() >= 8 {
                        if self.is_reverse {
                            read2 = Some(fq);
                        } else {
                            read1 = Some(fq);
                        }
                    }
                }
            }
        });
        Ok(AnnotatedFastq {
            barcode,
            umi,
            read1,
            read2,
        })
    }
}

#[derive(Debug)]
pub struct Barcode {
    pub raw: fastq::Record,
    pub corrected: Option<Vec<u8>>,
}

impl Barcode {
    pub fn extend(&mut self, other: &Self) {
        extend_fastq_record(&mut self.raw, &other.raw);
        if let Some(c2) = &other.corrected {
            if let Some(c1) = &mut self.corrected {
                c1.extend_from_slice(c2);
            }
        } else {
            self.corrected = None;
        }
    }
}

pub type UMI = fastq::Record;

/// An annotated fastq record with barcode, UMI, and sequence.
#[derive(Debug)]
pub struct AnnotatedFastq {
    pub barcode: Option<Barcode>,
    pub umi: Option<UMI>,
    pub read1: Option<fastq::Record>,
    pub read2: Option<fastq::Record>,
}

impl AnnotatedFastq {
    /// The total number of bases, including read1 and read2, in the record.
    pub fn len(&self) -> usize {
        self.read1.as_ref().map_or(0, |x| x.sequence().len())
            + self.read2.as_ref().map_or(0, |x| x.sequence().len())
    }
    pub fn is_empty(&self) -> bool {
        self.read1.is_none() && self.read2.is_none()
    }
}

impl AnnotatedFastq {
    pub fn join(&mut self, other: Self) {
        if let Some(bc) = &mut self.barcode {
            if let Some(x) = other.barcode.as_ref() {
                bc.extend(x)
            }
        } else {
            self.barcode = other.barcode;
        }

        if let Some(umi) = &mut self.umi {
            if let Some(x) = other.umi.as_ref() {
                extend_fastq_record(umi, x)
            }
        } else {
            self.umi = other.umi;
        }

        if self.read1.is_some() {
            if other.read1.is_some() {
                panic!("Read1 already exists");
            }
        } else {
            self.read1 = other.read1;
        }

        if self.read2.is_some() {
            if other.read2.is_some() {
                panic!("Read2 already exists");
            }
        } else {
            self.read2 = other.read2;
        }
    }
}

fn slice_fastq_record(record: &fastq::Record, start: usize, end: usize) -> fastq::Record {
    let end = end.min(record.sequence().len());
    fastq::Record::new(
        record.definition().clone(),
        &record.sequence()[start..end],
        record.quality_scores().get(start..end).unwrap(),
    )
}

pub fn extend_fastq_record(this: &mut fastq::Record, other: &fastq::Record) {
    this.sequence_mut().extend_from_slice(other.sequence());
    this.quality_scores_mut()
        .extend_from_slice(other.quality_scores());
}

pub struct NameCollatedRecords<'a, R> {
    records: bam::io::reader::Records<'a, R>,
    prev_record: Option<(BString, bam::Record)>,
    checker: HashSet<BString>,
}

impl<'a, R: std::io::Read> NameCollatedRecords<'a, R> {
    pub fn new(records: bam::io::reader::Records<'a, R>) -> Self {
        Self {
            records,
            prev_record: None,
            checker: HashSet::new(),
        }
    }

    fn check(&mut self, name: &BString) {
        assert!(
            !self.checker.contains(name),
            "bam file must be name collated or name sorted"
        );
        self.checker.insert(name.to_owned());
    }
}

impl<'a, R: std::io::Read> Iterator for NameCollatedRecords<'a, R> {
    type Item = (MultiMap<bam::Record>, MultiMap<bam::Record>);

    fn next(&mut self) -> Option<Self::Item> {
        let record = self.records.next()?.unwrap();
        let name = record.name().unwrap().to_owned();
        if let Some((prev_name, prev_record)) = self.prev_record.take() {
            if name == prev_name {
                Some((prev_record.into(), record.into()))
            } else {
                panic!(
                    "Expecting paired end reads with the same name, found {} and {}",
                    prev_name, name
                );
            }
        } else {
            self.check(&name);
            self.prev_record = Some((name, record));
            self.next()
        }
    }
}

#[cfg(test)]
mod tests {
    use bwa_mem2::{AlignerOpts, BurrowsWheelerAligner, FMIndex, PairedEndStats};
    use env_logger;

    use super::*;

    #[test]
    fn test_seqspec_io() {
        // Initialize logging with debug level
        env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("debug"))
            .init();

        let bwa_index = "/data/Public/BWA_MEM2_index/GRCh38";
        let seqspec = "/data/kzhang/dev/PreCellar/test/seqspec.yaml";
        let spec = Assay::from_path(seqspec).unwrap();
        let mut aligner = BurrowsWheelerAligner::new(
            FMIndex::read(bwa_index).unwrap(),
            AlignerOpts::default(),
            PairedEndStats::default(),
        );
        let mut processor = FastqProcessor::new(spec).with_modality(Modality::ATAC);

        processor
            .gen_barcoded_alignments(&mut aligner, 4, 40000)
            .take(6)
            .for_each(|x| {
                println!("{:?}", x);
            });

        println!("{}", processor.get_report());
    }
    
    #[test]
    fn test_barcode_filtering_control() {
        // Create a test whitelist with no predefined barcodes
        let mut whitelist = Whitelist::empty();
        
        // Add a mix of high and low count barcodes
        for i in 0..20 {
            let count = 1000 - i * 50;
            let barcode = format!("CELL_{:03}", i).into_bytes();
            let quality = vec![b'F'; barcode.len()];
            
            for _ in 0..count {
                whitelist.count_barcode(&barcode, &quality);
            }
        }
        
        // Add background noise
        for i in 0..100 {
            let count = 10;
            let barcode = format!("BG_{:03}", i).into_bytes();
            let quality = vec![b'F'; barcode.len()];
            
            for _ in 0..count {
                whitelist.count_barcode(&barcode, &quality);
            }
        }
        
        // Check barcode count before filtering
        let initial_count = whitelist.num_seen_barcodes();
        
        // Create a copy for comparison
        let mut whitelist_no_filtering = whitelist.clone();
        
        // Apply the filtering
        let results = filter_cellular_barcodes_ordmag_advanced(
            &mut whitelist,
            Some(20),
            None,
            None,
            Some(0.99),
            Some(10)
        );
        
        // Verify filtering worked
        assert!(whitelist.num_seen_barcodes() < initial_count, 
            "Filtering should reduce the number of barcodes");
        assert!(whitelist.num_seen_barcodes() <= results.filtered_bcs,
            "Number of barcodes should match the filtering results");
        
        // Verify we can skip filtering by not calling the function
        assert_eq!(whitelist_no_filtering.num_seen_barcodes(), initial_count,
            "Whitelist without filtering should maintain all barcodes");
    }
}