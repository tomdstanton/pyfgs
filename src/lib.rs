use pyo3::prelude::*;
use std::path::PathBuf;
use std::fs::File;

use seq_io::fasta::{Reader, Record as FastaRecord};
use seq_io::fastq::{self, Record as FastqRecord};

// Import directly from the crates.io library (v1.1.0)
use frag_gene_scan_rs::hmm::{self, Global, Local};
use frag_gene_scan_rs::viterbi::viterbi;
use frag_gene_scan_rs::dna::{Nuc, count_cg_content};

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

    fn __next__(mut slf: PyRefMut<'_, Self>) -> PyResult<Option<(String, String)>> {
        match slf.reader.next() {
            Some(Ok(record)) => {
                let head = String::from_utf8_lossy(record.id_bytes()).into_owned();
                let seq = String::from_utf8_lossy(record.seq()).into_owned();
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
    fn __next__(mut slf: PyRefMut<'_, Self>) -> PyResult<Option<(String, String, String)>> {
        match slf.reader.next() {
            Some(Ok(record)) => {
                let head = String::from_utf8_lossy(record.id_bytes()).into_owned();
                let seq = String::from_utf8_lossy(record.seq()).into_owned();
                let qual = String::from_utf8_lossy(record.qual()).into_owned();
                Ok(Some((head, seq, qual)))
            }
            Some(Err(e)) => Err(pyo3::exceptions::PyIOError::new_err(e.to_string())),
            None => Ok(None),
        }
    }
}

/// Represents a single predicted Open Reading Frame (ORF).
///
/// Example:
///     ```python
///     for gene in genes:
///         print(f"Start: {gene.start}, End: {gene.end}, Strand: {gene.strand}, Frame: {gene.frame}, Score: {gene.score}")
///     ```
#[pyclass]
pub struct Gene {
    /// int: The 1-based start position of the gene.
    #[pyo3(get)]
    pub start: usize,
    /// int: The 1-based end position of the gene.
    #[pyo3(get)]
    pub end: usize,
    /// str: The strand of the gene ('+' or '-').
    #[pyo3(get)]
    pub strand: String,
    /// int: The reading frame (1, 2, or 3).
    #[pyo3(get)]
    pub frame: usize,
    /// float: The Viterbi score of the gene prediction.
    #[pyo3(get)]
    pub score: f64,
}

#[pymethods]
impl Gene {
    fn __repr__(&self) -> String {
        format!("<Gene start={} end={} strand='{}' frame={} score={:.2}>",
                self.start, self.end, self.strand, self.frame, self.score)
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
    #[pyo3(signature = (model, whole_genome=false))]
    fn new(model: Model, whole_genome: bool) -> PyResult<Self> {
        let model_name = model.as_str();

        let (global, locals) = hmm::get_train_from_file(
            PathBuf::from(""),
            PathBuf::from(model_name)
        ).map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(
            format!("Failed to load internal model '{}'. Error: {}", model_name, e)
        ))?;

        Ok(GeneFinder {
            global,
            locals,
            whole_genome
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
    fn find_genes(&self, py: Python, header: String, sequence: String) -> PyResult<Vec<Gene>> {
        // In PyO3 0.28+, allow_threads is renamed to detach
        let results = py.detach(|| {
            // Because crates.io v1.1.0 doesn't export dna(), we do the string parsing directly
            let nuc_seq: Vec<Nuc> = sequence
                .bytes()
                .map(|b| b.to_ascii_uppercase())
                .map(Nuc::from)
                .collect();

            let head_bytes = header.into_bytes();

            // v1.1.0 viterbi requires us to manually select the local model based on CG content
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
                    start: g.start,
                    end: g.end,
                    strand: if g.forward_strand { "+".to_string() } else { "-".to_string() },
                    frame: g.frame,
                    score: g.score
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