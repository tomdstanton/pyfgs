use std::path::PathBuf;
use numpy::PyArray1;

use pyo3::prelude::*;
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pyclass_enum, gen_stub_pymethods};
use pyo3::types::PyBytes;

use frag_gene_scan_rs::hmm::{self, Global, Local};
use frag_gene_scan_rs::viterbi::viterbi;
use frag_gene_scan_rs::dna::{Nuc, count_cg_content, trinucleotide, CODON_CODE, ANTI_CODON_CODE};

use rayon::prelude::*;

use crate::hmm::{HmmGlobal, HmmLocal};

/// The available sequencing error models for FragGeneScanRs.
#[gen_stub_pyclass_enum]
#[pyclass(eq, eq_int, from_py_object)]
#[derive(PartialEq, Clone, Debug)]
pub enum Model {
    Illumina1,
    Illumina5,
    Illumina10,
    Sanger5,
    Sanger10,
    Pyro454_5,
    Pyro454_10,
    Pyro454_30,
    Complete,
}

impl Model {
    fn as_str(&self) -> &'static str {
        match self {
            Model::Illumina1 => "illumina_1",
            Model::Illumina5 => "illumina_5",
            Model::Illumina10 => "illumina_10",
            Model::Sanger5 => "sanger_5",
            Model::Sanger10 => "sanger_10",
            Model::Pyro454_5 => "454_5",
            Model::Pyro454_10 => "454_10",
            Model::Pyro454_30 => "454_30",
            Model::Complete => "complete",
        }
    }
}

// 1. The pure data struct (No methods allowed in here!)
#[gen_stub_pyclass]
#[pyclass]
pub struct Gene {
    #[pyo3(get)]
    pub start: usize,
    #[pyo3(get)]
    pub end: usize,
    #[pyo3(get)]
    pub strand: i8,
    #[pyo3(get)]
    pub frame: usize,
    #[pyo3(get)]
    pub score: f64,

    // Expose the crate's built-in indel vectors directly!
    #[pyo3(get)]
    pub insertions: Vec<usize>,
    #[pyo3(get)]
    pub deletions: Vec<usize>,

    // Internal fields hidden from Python used for lazy evaluation
    pub(crate) dna_nucs: Vec<Nuc>,
    pub(crate) forward_strand: bool,
    pub(crate) whole_genome: bool,
}

// 2. Internal Rust methods (Not exposed directly to Python)
impl Gene {
    #[inline]
    pub fn get_dna_bytes(&self) -> Vec<u8> {
        if self.forward_strand {
            self.dna_nucs.iter().map(|&n| u8::from(n)).collect()
        } else {
            self.dna_nucs.iter().rev().map(|&n| u8::from(n.rc())).collect()
        }
    }

    #[inline]
    pub fn get_aa_bytes(&self) -> Vec<u8> {
        let mut protein: Vec<u8> = if self.forward_strand {
            self.dna_nucs.chunks_exact(3)
                .map(|c| trinucleotide(c).map(|i| CODON_CODE[i]).unwrap_or(b'X'))
                .collect()
        } else {
            self.dna_nucs.rchunks_exact(3)
                .map(|c| trinucleotide(c).map(|i| ANTI_CODON_CODE[i]).unwrap_or(b'X'))
                .collect()
        };

        if protein.last() == Some(&b'*') {
            protein.pop();
        }

        if self.whole_genome {
            if self.forward_strand {
                let s = trinucleotide(&self.dna_nucs);
                if s == trinucleotide(&[Nuc::G, Nuc::T, Nuc::G]) || s == trinucleotide(&[Nuc::T, Nuc::T, Nuc::G]) {
                    protein[0] = b'M';
                }
            } else if self.dna_nucs.len() >= 3 {
                let s = trinucleotide(&self.dna_nucs[self.dna_nucs.len() - 3..]);
                if s == trinucleotide(&[Nuc::C, Nuc::A, Nuc::C]) || s == trinucleotide(&[Nuc::C, Nuc::A, Nuc::A]) {
                    protein[0] = b'M';
                }
            }
        }
        protein
    }

    pub fn extract_mutations(&self, sequence: &[u8]) -> Vec<Mutation> {
        let mut mutations = Vec::new();

        if self.insertions.is_empty() && self.deletions.is_empty() {
            return mutations;
        }

        // Helper to calculate the 1-based codon index dynamically
        let get_codon_idx = |pos: usize| -> usize {
            if self.forward_strand {
                let ins_count = self.insertions.partition_point(|&p| p < pos);
                let del_count = self.deletions.partition_point(|&p| p <= pos);
                let effective_cds_pos = (pos.saturating_sub(self.start)) - ins_count + del_count;
                (effective_cds_pos / 3) + 1
            } else {
                let ins_count = self.insertions.len() - self.insertions.partition_point(|&p| p <= pos);
                let del_count = self.deletions.len() - self.deletions.partition_point(|&p| p < pos);
                let effective_cds_pos = (self.end.saturating_sub(1).saturating_sub(pos)) - ins_count + del_count;
                (effective_cds_pos / 3) + 1
            }
        };

        // 1. Process Insertions
        for &pos in &self.insertions {
            if pos == 0 { continue; }
            let anchor_pos = pos - 1;

            let ref_allele = String::from_utf8_lossy(&sequence[anchor_pos..=pos]).to_ascii_uppercase();
            let alt_allele = (sequence[anchor_pos] as char).to_ascii_uppercase().to_string();
            let codon_idx = get_codon_idx(pos);
            let annotation = format!("AA:X->fs|DNA:{}->{}", ref_allele, alt_allele);

            mutations.push(Mutation {
                pos: anchor_pos + 1,
                mut_type: "ins".to_string(),
                ref_allele,
                alt_allele,
                codon_idx,
                annotation,
            });
        }

        // 2. Process Deletions
        let mut sorted_del = self.deletions.clone();
        sorted_del.sort_unstable();
        for &pos in &sorted_del {
            if pos == 0 { continue; }
            let anchor_pos = pos - 1;

            let ref_allele = (sequence[anchor_pos] as char).to_ascii_uppercase().to_string();
            let alt_allele = format!("{}N", ref_allele);
            let codon_idx = get_codon_idx(pos);
            let annotation = format!("AA:X->fs|DNA:{}->{}", ref_allele, alt_allele);

            mutations.push(Mutation {
                pos: anchor_pos + 1,
                mut_type: "del".to_string(),
                ref_allele,
                alt_allele,
                codon_idx,
                annotation,
            });
        }

        mutations
    }
}

// 3. The Python-facing API (wraps the internal methods)
#[gen_stub_pymethods]
#[pymethods]
impl Gene {
    fn __repr__(&self) -> String {
        format!("<Gene start={} end={} strand={} frame={} score={:.2}>",
                self.start, self.end, self.strand, self.frame, self.score)
    }

    pub fn sequence<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        PyBytes::new(py, &self.get_dna_bytes())
    }

    pub fn translation<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        PyBytes::new(py, &self.get_aa_bytes())
    }

    /// Extracts all frameshift mutations as structured objects.
    /// Requires the original raw sequence to determine VCF anchored alleles.
    #[pyo3(name = "mutations")]
    pub fn py_mutations(&self, sequence: &Bound<'_, PyBytes>) -> PyResult<Vec<Mutation>> {
        let sequence = sequence.as_bytes();
        Ok(self.extract_mutations(sequence))
    }
}

/// A vectorized structure representing a batch of gene predictions.
///
/// This structure uses a Structure-of-Arrays (SoA) layout. Fields that have variable 
/// counts per gene (e.g. insertions and deletions) are exposed as flat arrays alongside 
/// parallel "offsets" arrays defining the slice boundaries (i.e. ragged arrays).
#[gen_stub_pyclass]
#[pyclass]
pub struct GeneBatch {
    /// Array of shape (N,) containing the index of the sequence in the batch for each gene.
    #[pyo3(get)] pub sequence_indices: Py<PyArray1<u64>>,
    
    /// Array of shape (N,) containing the start coordinates of each gene.
    #[pyo3(get)] pub starts: Py<PyArray1<u64>>,
    
    /// Array of shape (N,) containing the end coordinates of each gene.
    #[pyo3(get)] pub ends: Py<PyArray1<u64>>,
    
    /// Array of shape (N,) containing the strands of each gene (1 for forward, -1 for reverse).
    #[pyo3(get)] pub strands: Py<PyArray1<i8>>,
    
    /// Array of shape (N,) containing the frame of each gene.
    #[pyo3(get)] pub frames: Py<PyArray1<u64>>,
    
    /// Array of shape (N,) containing the prediction scores of each gene.
    #[pyo3(get)] pub scores: Py<PyArray1<f64>>,
    
    /// Flat array containing all insertion coordinates across all genes.
    #[pyo3(get)] pub insertions_flat: Py<PyArray1<u64>>,
    
    /// Array of shape (N+1,) containing the slice offsets for `insertions_flat`.
    #[pyo3(get)] pub insertions_offsets: Py<PyArray1<u64>>,
    
    /// Flat array containing all deletion coordinates across all genes.
    #[pyo3(get)] pub deletions_flat: Py<PyArray1<u64>>,
    
    /// Array of shape (N+1,) containing the slice offsets for `deletions_flat`.
    #[pyo3(get)] pub deletions_offsets: Py<PyArray1<u64>>,
}

/// The main engine for finding genes, holding the HMM in memory.
#[gen_stub_pyclass]
#[pyclass]
pub struct GeneFinder {
    pub(crate) global: Box<Global>,
    pub(crate) locals: Vec<Local>,
    pub(crate) whole_genome: bool,
}

unsafe impl Send for GeneFinder {}
unsafe impl Sync for GeneFinder {}

#[gen_stub_pymethods]
#[pymethods]
impl GeneFinder {
    /// Initialize the GeneFinder.
    ///
    /// Args:
    ///     model (Model): The sequencing error model to use (e.g., pyfgs.Model.Illumina5).
    ///     whole_genome (bool, optional): Set to True if analyzing complete genomic sequences
    ///         rather than short reads. Defaults to False.
    #[new]
    #[pyo3(signature = (model, whole_genome=None))]
    fn new(model: Model, whole_genome: Option<bool>) -> PyResult<Self> {
        let model_name = model.as_str();

        // Infer whole_genome if the user didn't explicitly provide it
        let is_whole_genome = whole_genome.unwrap_or(model == Model::Complete);

        let (global, locals) = hmm::get_train_from_file(
            PathBuf::from(""),
            PathBuf::from(model_name)
        ).map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(
            format!("Failed to load internal model '{}'. Error: {}", model_name, e)
        ))?;

        Ok(GeneFinder {
            global,
            locals,
            whole_genome: is_whole_genome
        })
    }

    /// Exposes the underlying global HMM parameters (read-only).
    #[pyo3(name = "global_model")]
    pub fn py_global_model(&self) -> HmmGlobal {
        HmmGlobal::new(&self.global)
    }

    /// Exposes the underlying local (GC-specific) HMM parameters (read-only).
    #[pyo3(name = "local_models")]
    pub fn py_local_models(&self) -> Vec<HmmLocal> {
        self.locals.iter().map(|l| HmmLocal::new(l)).collect()
    }

    /// Predict open reading frames in a given DNA sequence.
    ///
    /// This method releases the GIL, allowing for safe multi-threading
    /// across multiple CPU cores.
    ///
    /// Args:
    ///     sequence (bytes): The raw nucleotide bytes.
    ///
    /// Returns:
    ///     List[Gene]: A list of predicted Gene objects.
    fn find_genes(&self, py: Python, sequence: &Bound<'_, PyBytes>) -> PyResult<Vec<Gene>> {
        let sequence = sequence.as_bytes();
        let nuc_seq: Vec<Nuc> = sequence
            .iter()
            .map(|&b| Nuc::from(b))
            .collect();

        // 2. Run the heavy HMM inference without the GIL
        let prediction = py.detach(|| {
            let cg = count_cg_content(&nuc_seq);
            let local_model = &self.locals[cg];

            // Pass an empty vector to satisfy the crate's signature without allocating memory
            viterbi(
                &self.global,
                local_model,
                Vec::new(), 
                nuc_seq,
                self.whole_genome
            )
        });

        let genes = prediction.genes.into_iter().map(|g| {
            Gene {
                start: g.start.saturating_sub(1),
                end: g.end,
                strand: if g.forward_strand { 1 } else { -1 },
                frame: g.frame,
                score: g.score,
                insertions: g.inserted.into_iter().map(|i| i.saturating_sub(1)).collect(),
                deletions: g.deleted.into_iter().map(|i| i.saturating_sub(1)).collect(),
                dna_nucs: g.dna,
                forward_strand: g.forward_strand,
                whole_genome: self.whole_genome,
            }
        }).collect();

        Ok(genes)
    }

    /// Predict open reading frames for a batch of DNA sequences using Rayon.
    ///
    /// Args:
    ///     sequences (list[bytes]): A list of raw nucleotide bytes.
    ///
    /// Returns:
    ///     List[List[Gene]]: A list of predicted Gene objects for each sequence.
    #[pyo3(name = "find_genes_batch")]
    fn find_genes_batch(
        &self,
        py: Python<'_>,
        sequences: Vec<Bound<'_, PyBytes>>,
    ) -> PyResult<GeneBatch> {
        struct RawQuery {
            seq_ptr: *const u8,
            seq_len: usize,
        }

        unsafe impl Send for RawQuery {}
        unsafe impl Sync for RawQuery {}

        let mut raw_queries = Vec::with_capacity(sequences.len());
        for seq in sequences {
            raw_queries.push(RawQuery {
                seq_ptr: seq.as_bytes().as_ptr(),
                seq_len: seq.as_bytes().len(),
            });
        }

        let predictions: Vec<Vec<_>> = py.detach(|| {
            raw_queries
                .par_iter()
                .map(|raw_q| {
                    let sequence = unsafe { std::slice::from_raw_parts(raw_q.seq_ptr, raw_q.seq_len) };
                    let nuc_seq: Vec<Nuc> = sequence
                        .iter()
                        .map(|&b| Nuc::from(b))
                        .collect();

                    let cg = count_cg_content(&nuc_seq);
                    let local_model = &self.locals[cg];

                    viterbi(
                        &self.global,
                        local_model,
                        Vec::new(), 
                        nuc_seq,
                        self.whole_genome
                    ).genes
                })
                .collect()
        });

        let mut total_genes = 0;
        let mut total_ins = 0;
        let mut total_del = 0;
        for genes in &predictions {
            total_genes += genes.len();
            for g in genes {
                total_ins += g.inserted.len();
                total_del += g.deleted.len();
            }
        }

        let mut sequence_indices = Vec::with_capacity(total_genes);
        let mut starts = Vec::with_capacity(total_genes);
        let mut ends = Vec::with_capacity(total_genes);
        let mut strands = Vec::with_capacity(total_genes);
        let mut frames = Vec::with_capacity(total_genes);
        let mut scores = Vec::with_capacity(total_genes);

        let mut insertions_flat = Vec::with_capacity(total_ins);
        let mut insertions_offsets = Vec::with_capacity(total_genes + 1);
        insertions_offsets.push(0);

        let mut deletions_flat = Vec::with_capacity(total_del);
        let mut deletions_offsets = Vec::with_capacity(total_genes + 1);
        deletions_offsets.push(0);

        for (seq_idx, genes) in predictions.into_iter().enumerate() {
            for g in genes {
                sequence_indices.push(seq_idx as u64);
                starts.push(g.start.saturating_sub(1) as u64);
                ends.push(g.end as u64);
                strands.push(if g.forward_strand { 1 } else { -1 });
                frames.push(g.frame as u64);
                scores.push(g.score);

                for i in g.inserted {
                    insertions_flat.push(i.saturating_sub(1) as u64);
                }
                insertions_offsets.push(insertions_flat.len() as u64);

                for d in g.deleted {
                    deletions_flat.push(d.saturating_sub(1) as u64);
                }
                deletions_offsets.push(deletions_flat.len() as u64);
            }
        }

        Ok(GeneBatch {
            sequence_indices: PyArray1::from_vec(py, sequence_indices).into(),
            starts: PyArray1::from_vec(py, starts).into(),
            ends: PyArray1::from_vec(py, ends).into(),
            strands: PyArray1::from_vec(py, strands).into(),
            frames: PyArray1::from_vec(py, frames).into(),
            scores: PyArray1::from_vec(py, scores).into(),
            insertions_flat: PyArray1::from_vec(py, insertions_flat).into(),
            insertions_offsets: PyArray1::from_vec(py, insertions_offsets).into(),
            deletions_flat: PyArray1::from_vec(py, deletions_flat).into(),
            deletions_offsets: PyArray1::from_vec(py, deletions_offsets).into(),
        })
    }

    /// Run the full CLI pipeline purely in Rust without python loop overhead.
    ///
    /// Args:
    ///     input_path (str): The path to the FASTA/FASTQ file or "-" for stdin.
    ///     is_fastq (bool): True if parsing as FASTQ.
    ///     outputs (Dict[str, str]): Output formats and paths.
    #[pyo3(signature = (input_path, is_fastq, outputs))]
    fn run_file(
        &self,
        py: Python<'_>,
        input_path: &str,
        is_fastq: bool,
        outputs: std::collections::HashMap<String, String>,
    ) -> PyResult<()> {
        crate::cli::run_file(py, self, input_path, is_fastq, outputs)
    }
}

/// Represents a frameshift mutation detected by the ab initio model.
#[gen_stub_pyclass]
#[pyclass(get_all, skip_from_py_object)]
#[derive(Clone, Debug)]
pub struct Mutation {
    /// 1-based VCF anchored position
    pub pos: usize,
    /// "ins" (assembly has extra base) or "del" (assembly missing base)
    pub mut_type: String,
    /// The VCF reference allele
    pub ref_allele: String,
    /// The VCF alternate allele (conceptual FGS fix)
    pub alt_allele: String,
    /// 1-based codon index where the frame breaks
    pub codon_idx: usize,
    /// Snippy-style text annotation
    pub annotation: String,
}

#[gen_stub_pymethods]
#[pymethods]
impl Mutation {
    fn __repr__(&self) -> String {
        format!(
            "<Mutation type='{}' pos={} codon={} mut='{}'>",
            self.mut_type, self.pos, self.codon_idx, self.annotation
        )
    }
}
