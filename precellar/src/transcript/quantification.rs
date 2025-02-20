use anndata::{
    data::utils::{from_csr_data, to_csr_data},
    AnnData, AnnDataOp,
};
use anndata_hdf5::H5;
use anyhow::Result;
use bed_utils::extsort::ExternalSorterBuilder;
use indexmap::IndexMap;
use itertools::Itertools;
use noodles::sam::Header;
use polars::df;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use serde::{Deserialize, Serialize};
use std::{collections::HashSet, path::PathBuf, sync::atomic::{AtomicUsize, Ordering}};

use crate::{align::MultiMapR, qc::GeneQuantQC, transcript::Gene};

use super::{
    annotate::AnnotationRegion, de_dups::count_unique_umi, AlignmentAnnotator, AnnotatedAlignment,
};

#[derive(Debug, Serialize, Deserialize)]
pub struct GeneAlignment {
    pub idx: usize,
    pub umi: Option<String>,
    pub align_type: AnnotationRegion,
}

#[derive(Debug)]
pub struct Quantifier {
    annotator: AlignmentAnnotator,
    genes: IndexMap<String, Gene>,
    temp_dir: Option<PathBuf>,
    chunk_size: usize,
    mito_genes: HashSet<usize>,
}

impl Quantifier {
    pub fn new(annotator: AlignmentAnnotator) -> Self {
        let genes = annotator
            .transcripts
            .iter()
            .map(|(_, t)| {
                let g = t.gene.clone();
                (g.id.clone(), g)
            })
            .collect();
        Self {
            annotator,
            genes,
            temp_dir: None,
            chunk_size: 50000000,
            mito_genes: HashSet::new(),
        }
    }

    pub fn add_mito_dna(&mut self, mito_chr: &str) {
        let iter = self.annotator.transcripts.iter().flat_map(|(_, t)| {
            if t.chrom == mito_chr {
                Some(self.genes.get_full(&t.gene.id).unwrap().0)
            } else {
                None
            }
        });
        self.mito_genes.extend(iter);
    }

    pub fn quantify<'a, I, P>(
        &'a self,
        header: &'a Header,
        records: I,
        output: P,
    ) -> Result<GeneQuantQC>
    where
        I: Iterator<Item = Vec<(Option<MultiMapR>, Option<MultiMapR>)>> + 'a,
        P: AsRef<std::path::Path>,
    {
        eprintln!("Starting quantification...");
        eprintln!("Number of genes: {}", self.genes.len());
        eprintln!("Number of mitochondrial genes: {}", self.mito_genes.len());

        // Add debug counter for initial records
        let mut debug_record_count = 0;
        let records = records.inspect(|recs| {
            if debug_record_count < 5 {
                eprintln!("\nInput record batch {}:", debug_record_count);
                for (i, (r1, r2)) in recs.iter().take(3).enumerate() {
                    eprintln!("  Record {}: r1: {}, r2: {}", 
                        i,
                        r1.as_ref().map_or("None", |_| "Some"),
                        r2.as_ref().map_or("None", |_| "Some"));
                }
                debug_record_count += 1;
            }
        });

        let mut qc = GeneQuantQC::default();
        let adata: AnnData<H5> = AnnData::new(output)?;
        let num_cols = self.genes.len();

        let mut barcodes = Vec::new();
        let mut exon_count: Vec<u64> = Vec::new();
        let mut intron_count: Vec<u64> = Vec::new();
        let mut mito_count: Vec<u64> = Vec::new();

        {
            // Use atomic counters for thread-safe counting
            let record_count = AtomicUsize::new(0);
            let alignment_count = AtomicUsize::new(0);
            
            eprintln!("Processing alignments...");
            let mut debug_alignment_count = 0;
            let tx_alignments = records.flat_map(|recs| {
                let batch_size = recs.len();
                record_count.fetch_add(batch_size, Ordering::Relaxed);
                qc.total_reads += batch_size as u64;
                
                let alignments: Vec<_> = recs.into_par_iter()
                    .filter_map(|(r1, r2)| {
                        let result = self.make_gene_alignment(header, r1, r2);
                        if result.is_some() {
                            alignment_count.fetch_add(1, Ordering::Relaxed);
                        }
                        result
                    })
                    .collect();
                
                // Debug print first few valid alignments
                if debug_alignment_count < 5 && !alignments.is_empty() {
                    eprintln!("\nValid alignments batch {}:", debug_alignment_count);
                    for (i, (barcode, align)) in alignments.iter().take(3).enumerate() {
                        eprintln!("  Alignment {}: barcode: {}, gene_idx: {}, align_type: {:?}, umi: {:?}", 
                            i, barcode, align.idx, align.align_type, align.umi);
                    }
                    debug_alignment_count += 1;
                }

                alignments
            });

            let final_record_count = record_count.load(Ordering::Relaxed);
            let final_alignment_count = alignment_count.load(Ordering::Relaxed);
            eprintln!("Total records processed: {}", final_record_count);
            eprintln!("Valid alignments found: {}", final_alignment_count);

            if final_alignment_count == 0 {
                return Err(anyhow::anyhow!("No valid alignments found. Check if:
                    1. Input records are not empty
                    2. Reads have valid barcodes
                    3. Reads align to annotated genes
                    4. UMI information is present"));
            }

            eprintln!("Sorting alignments by barcode...");
            let tx_alignments_chunks =
                sort_alignments(tx_alignments, self.temp_dir.as_ref(), self.chunk_size)
                    .chunk_by(|x| x.0.clone());

            eprintln!("Processing barcodes...");
            let barcode_count = AtomicUsize::new(0);
            let mut debug_barcode_count = 0;
            let tx_alignments_chunks = tx_alignments_chunks
                .into_iter()
                .map(|(barcode, group)| {
                    let current_count = barcode_count.fetch_add(1, Ordering::Relaxed);
                    if current_count % 10000 == 0 {
                        eprintln!("Processed {} barcodes", current_count);
                    }
                    barcodes.push(barcode);
                    group.map(|x| x.1).collect::<Vec<_>>()
                })
                .chunks(500);

            let final_barcode_count = barcode_count.load(Ordering::Relaxed);
            eprintln!("Found {} unique barcodes", final_barcode_count);
            if final_barcode_count == 0 {
                return Err(anyhow::anyhow!("No valid barcodes found in the alignments"));
            }

            eprintln!("Building count matrix...");
            let counts = tx_alignments_chunks.into_iter().map(|chunk| {
                let results: Vec<_> = chunk
                    .collect::<Vec<_>>()
                    .into_par_iter()
                    .map(|alignments| count_unique_umi(alignments, &self.mito_genes))
                    .collect();

                results.iter().for_each(|r| {
                    qc.total_umi += r.total_umi;
                    qc.unique_umi += r.unique_umi;
                    exon_count.push(r.uniq_exon);
                    intron_count.push(r.uniq_intron);
                    mito_count.push(r.uniq_mito);
                });

                let (nrows, ncols, indptr, indices, data) = to_csr_data(
                    results
                        .into_iter()
                        .map(|r| r.into_counts().collect()),
                    num_cols,
                );
                from_csr_data(nrows, ncols, indptr, indices, data).unwrap()
            });

            eprintln!("Writing count matrix to file...");
            adata.set_x_from_iter(counts)?;
        }

        eprintln!("Writing metadata...");
        let num_barcodes = barcodes.len();  // Store length before moving
        adata.set_obs_names(barcodes.into())?;
        adata.set_obs(df!(
            "exon_count" => exon_count,
            "intron_count" => intron_count,
            "mitochondrial_count" => mito_count
        )?)?;

        adata.set_var_names(self.genes.values().map(|g| g.id.clone()).collect())?;
        adata.set_var(df!(
            "gene_name" => self.genes.values().map(|g| g.name.clone()).collect::<Vec<_>>()
        )?)?;

        eprintln!("Finalizing output...");
        adata.close()?;

        eprintln!("Quantification complete!");
        eprintln!("Total reads processed: {}", qc.total_reads);
        eprintln!("Total UMIs: {}", qc.total_umi);
        eprintln!("Unique UMIs: {}", qc.unique_umi);
        eprintln!("Number of barcodes: {}", num_barcodes);

        Ok(qc)
    }

    fn make_gene_alignment(
        &self,
        header: &Header,
        rec1: Option<MultiMapR>,
        rec2: Option<MultiMapR>,
    ) -> Option<(String, GeneAlignment)> {
        let barcode;
        let umi;
        let anno = if rec1.is_some() && rec2.is_some() {
            let rec1 = rec1.unwrap();
            let rec2 = rec2.unwrap();
            barcode = rec1.barcode().unwrap()?;
            umi = rec1.umi().unwrap();
            self.annotator.annotate_alignments_pe(header, rec1, rec2)
        } else {
            let rec = rec1.or(rec2).unwrap();
            barcode = rec.barcode().unwrap()?;
            umi = rec.umi().unwrap();
            self.annotator.annotate_alignments_se(header, rec)
        }?;

        let gene_id;
        let align_type;

        match anno {
            AnnotatedAlignment::PeMapped(a1, a2, anno) => {
                let genes = anno.genes;
                if genes.len() != 1 {
                    return None;
                }
                let gene = genes.iter().next().unwrap();
                align_type = match (a1.region, a2.region) {
                    (AnnotationRegion::Intronic, _) => AnnotationRegion::Intronic,
                    (_, AnnotationRegion::Intronic) => AnnotationRegion::Intronic,
                    _ => AnnotationRegion::Exonic,
                };
                gene_id = self.genes.get_full(&gene.id).unwrap().0;
            }
            AnnotatedAlignment::SeMapped(anno) => {
                let genes = anno.genes;
                if genes.len() != 1 {
                    return None;
                }
                let gene = genes.iter().next().unwrap();
                align_type = anno.region;
                gene_id = self.genes.get_full(&gene.id).unwrap().0;
            }
        }

        let alignment = GeneAlignment {
            idx: gene_id,
            umi,
            align_type,
        };
        Some((barcode, alignment))
    }
}

fn sort_alignments<I, P>(
    alignments: I,
    temp_dir: Option<P>,
    chunk_size: usize,
) -> impl ExactSizeIterator<Item = (String, GeneAlignment)>
where
    I: Iterator<Item = (String, GeneAlignment)>,
    P: AsRef<std::path::Path>,
{
    let mut sorter = ExternalSorterBuilder::new()
        .with_chunk_size(chunk_size)
        .with_compression(2);
    if let Some(tmp) = temp_dir {
        sorter = sorter.with_tmp_dir(tmp);
    }
    sorter
        .build()
        .unwrap()
        .sort_by(alignments, |a, b| a.0.cmp(&b.0))
        .unwrap()
        .map(|x| x.unwrap())
}
