use pyo3::prelude::*;
use pyo3::types::PyBytes;
use std::path::PathBuf;
use std::fs::File;

use seq_io::fasta::{Reader, Record as FastaRecord};
use seq_io::fastq::{self, Record as FastqRecord};

// Import directly from the crates.io library (v1.1.0)
use frag_gene_scan_rs::hmm::{self, Global, Local};
use frag_gene_scan_rs::viterbi::viterbi;
use frag_gene_scan_rs::dna::{Nuc, count_cg_content, trinucleotide, CODON_CODE, ANTI_CODON_CODE};

use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

/// The available sequencing error models for FragGeneScanRs.
///
/// Attributes:
///     Illumina1: Illumina sequencing reads with about 0.1% error rate.
///     Illumina5: Illumina sequencing reads with about 0.5% error rate.
///     Illumina10: Illumina sequencing reads with about 1% error rate.
///     Sanger5: Sanger sequencing reads with about 0.5% error rate.
///     Sanger10: Sanger sequencing reads with about 1% error rate.
///     Pyro454_5: 454 pyrosequencing reads with about 0.5% error rate.
///     Pyro454_10: 454 pyrosequencing reads with about 1% error rate.
///     Pyro454_30: 454 pyrosequencing reads with about 3% error rate.
///     Complete: Complete genomic sequences or short sequence reads without sequencing error.
///
/// Example:
///     ```python
///     from pyfgs import Model
///     model = Model.Complete
///     ```
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

/// A memory-efficient FASTA parser yielding (header, sequence).
///
/// Example:
///     ```python
///     from pyfgs import FastaReader
///     
///     reader = FastaReader("genome.fasta")
///     for header, sequence in reader:
///         print(f">{header}\n{sequence}")
///     ```
#[pyclass(unsendable)]
pub struct FastaReader {
    reader: Reader<File>,
}

#[pymethods]
impl FastaReader {
    /// Open a FASTA file for reading.
    ///
    /// Args:
    ///     path (str): The path to the FASTA file.
    #[new]
    fn new(path: String) -> PyResult<Self> {
        let file = File::open(&path)
            .map_err(|e| pyo3::exceptions::PyFileNotFoundError::new_err(format!("Could not open {}: {}", path, e)))?;
        Ok(Self { reader: Reader::new(file) })
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    // Pass `py: Python<'py>` to interact with Python memory
    fn __next__<'py>(mut slf: PyRefMut<'_, Self>, py: Python<'py>) -> PyResult<Option<(Bound<'py, PyBytes>, Bound<'py, PyBytes>)>> {
        match slf.reader.next() {
            Some(Ok(record)) => {
                let head = PyBytes::new(py, record.id_bytes());
                let seq = PyBytes::new(py, record.seq());
                Ok(Some((head, seq)))
            }
            Some(Err(e)) => Err(pyo3::exceptions::PyIOError::new_err(e.to_string())),
            None => Ok(None),
        }
    }
}

/// A memory-efficient FASTQ parser yielding (header, sequence, qualities).
///
/// Example:
///     ```python
///     from pyfgs import FastqReader
///     
///     reader = FastqReader("reads.fastq")
///     for header, sequence, qualities in reader:
///         print(f"@{header}\n{sequence}\n+\n{qualities}")
///     ```
#[pyclass(unsendable)]
pub struct FastqReader {
    reader: fastq::Reader<File>,
}

#[pymethods]
impl FastqReader {
    /// Open a FASTQ file for reading.
    ///
    /// Args:
    ///     path (str): The path to the FASTQ file.
    #[new]
    fn new(path: String) -> PyResult<Self> {
        let file = File::open(&path)
            .map_err(|e| pyo3::exceptions::PyFileNotFoundError::new_err(format!("Could not open {}: {}", path, e)))?;
        Ok(Self { reader: fastq::Reader::new(file) })
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    /// Yields a tuple of (header, sequence, quality_string)
    fn __next__<'py>(mut slf: PyRefMut<'_, Self>, py: Python<'py>) -> PyResult<Option<(Bound<'py, PyBytes>, Bound<'py, PyBytes>, Bound<'py, PyBytes>)>> {
        match slf.reader.next() {
            Some(Ok(record)) => {
                let head = PyBytes::new(py, record.id_bytes());
                let seq = PyBytes::new(py, record.seq());
                let qual = PyBytes::new(py, record.qual());
                Ok(Some((head, seq, qual)))
            }
            Some(Err(e)) => Err(pyo3::exceptions::PyIOError::new_err(e.to_string())),
            None => Ok(None),
        }
    }
}

/// Represents a single predicted Open Reading Frame (ORF).
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
    dna_nucs: Vec<Nuc>,
    forward_strand: bool,
    whole_genome: bool,
}

#[pymethods]
impl Gene {
    fn __repr__(&self) -> String {
        format!("<Gene start={} end={} strand={} frame={} score={:.2}>",
                self.start, self.end, self.strand, self.frame, self.score)
    }

    /// Lazily computes and returns the raw nucleotide sequence as bytes.
    pub fn sequence<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        let dna_bytes: Vec<u8> = if self.forward_strand {
            self.dna_nucs.iter().map(|&n| u8::from(n)).collect()
        } else {
            self.dna_nucs.iter().rev().map(|&n| u8::from(n.rc())).collect()
        };

        PyBytes::new(py, &dna_bytes)
    }

    /// Lazily computes and returns the translated amino acid sequence as bytes.
    pub fn translation<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        let mut protein_bytes: Vec<u8> = if self.forward_strand {
            self.dna_nucs.chunks_exact(3)
                .map(|c| trinucleotide(c).map(|i| CODON_CODE[i]).unwrap_or(b'X'))
                .collect()
        } else {
            self.dna_nucs.rchunks_exact(3)
                .map(|c| trinucleotide(c).map(|i| ANTI_CODON_CODE[i]).unwrap_or(b'X'))
                .collect()
        };

        if protein_bytes.last() == Some(&b'*') {
            protein_bytes.pop();
        }

        if self.whole_genome {
            if self.forward_strand {
                let s = trinucleotide(&self.dna_nucs);
                if s == trinucleotide(&[Nuc::G, Nuc::T, Nuc::G]) || s == trinucleotide(&[Nuc::T, Nuc::T, Nuc::G]) {
                    protein_bytes[0] = b'M';
                }
            } else if self.dna_nucs.len() >= 3 {
                let s = trinucleotide(&self.dna_nucs[self.dna_nucs.len() - 3..]);
                if s == trinucleotide(&[Nuc::C, Nuc::A, Nuc::C]) || s == trinucleotide(&[Nuc::C, Nuc::A, Nuc::A]) {
                    protein_bytes[0] = b'M';
                }
            }
        }

        PyBytes::new(py, &protein_bytes)
    }
}

/// The main engine for finding genes, holding the HMM in memory.
///
/// Example:
///     ```python
///     from pyfgs import GeneFinder, Model
///     
///     finder = GeneFinder(Model.Complete, whole_genome=True)
///     genes = finder.find_genes("seq1", "ATGCGTA...")
///     ```
#[pyclass]
pub struct GeneFinder {
    global: Box<Global>,
    locals: Vec<Local>,
    whole_genome: bool,
}

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

    /// Predict open reading frames in a given DNA sequence.
    ///
    /// This method releases the GIL, allowing for safe multi-threading
    /// across multiple CPU cores.
    ///
    /// Args:
    ///     header (str): The sequence identifier.
    ///     sequence (str): The raw nucleotide string.
    ///
    /// Returns:
    ///     List[Gene]: A list of predicted Gene objects.
    ///     
    /// Example:
    ///     ```python
    ///     genes = finder.find_genes("seq1", "ATGCGTACGTTAG")
    ///     ```
    // Changed String to &str for zero-copy borrowing from Python
    fn find_genes(&self, py: Python, header: Vec<u8>, sequence: Vec<u8>) -> PyResult<Vec<Gene>> {

        // 1. We skip string mapping and iterate the bytes directly
        // Build the Nuc vector directly from Python's memory.
        // We just eliminated an entire heap allocation step!
       let nuc_seq: Vec<Nuc> = sequence
            .iter()
            .filter(|&&b| !b.is_ascii_whitespace())
            .map(|&b| b.to_ascii_uppercase())
            .map(Nuc::from)
            .collect();

        // We still need to own the header to pass it to the background thread
        let head_bytes = header.to_vec();

        // 2. Move the owned data into the background thread (detaching the GIL)
        let results = py.detach(move || {
            let cg = count_cg_content(&nuc_seq);
            let local_model = &self.locals[cg];

            let prediction = viterbi(
                &self.global,
                local_model,
                head_bytes,
                nuc_seq,
                self.whole_genome
            );

            prediction.genes.into_iter().map(|g| {
                Gene {
                    start: g.start.saturating_sub(1),
                    end: g.end,
                    strand: if g.forward_strand { 1 } else { -1 },
                    frame: g.frame,
                    score: g.score,

                    // Make sure the positions are 0-based
                    insertions: g.inserted.into_iter().map(|i| i.saturating_sub(1)).collect(),
                    deletions: g.deleted.into_iter().map(|i| i.saturating_sub(1)).collect(),

                    dna_nucs: g.dna,
                    forward_strand: g.forward_strand,
                    whole_genome: self.whole_genome,
                }
            }).collect()
        });

        Ok(results)
    }
}

/// The Python module definition using the modern PyO3 Bound API
#[pymodule]
fn _pyfgs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Model>()?;
    m.add_class::<FastaReader>()?;
    m.add_class::<FastqReader>()?;
    m.add_class::<Gene>()?;
    m.add_class::<GeneFinder>()?;
    Ok(())
}