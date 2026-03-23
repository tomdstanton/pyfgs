use std::path::PathBuf;
use std::fs::File;
use std::io::{BufWriter, Write};

use pyo3::prelude::*;
use pyo3::types::PyType;
use pyo3::types::PyBytes;

use seq_io::fasta::{Reader, Record as FastaRecord};
use seq_io::fastq::{self, Record as FastqRecord};

use frag_gene_scan_rs::hmm::{self, Global, Local};
use frag_gene_scan_rs::viterbi::viterbi;
use frag_gene_scan_rs::dna::{Nuc, count_cg_content, trinucleotide, CODON_CODE, ANTI_CODON_CODE};

use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

/// The available sequencing error models for FragGeneScanRs.
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
#[pyclass(unsendable)]
pub struct FastqReader {
    reader: fastq::Reader<File>,
}

#[pymethods]
impl FastqReader {
    /// Open a FASTQ file for reading.
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
// 1. The pure data struct (No methods allowed in here!)
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
                (effective_cds_pos / 3) + 1 // FIXED: Added missing return statement here!
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
        // Re-added the clone/sort because the raw vectors are immutable in this context
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
    pub fn py_mutations(&self, sequence: &[u8]) -> PyResult<Vec<Mutation>> {
        Ok(self.extract_mutations(sequence))
    }
}

/// The main engine for finding genes, holding the HMM in memory.
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
    /// Predict open reading frames in a given DNA sequence.
    ///
    /// This method releases the GIL, allowing for safe multi-threading
    /// across multiple CPU cores.
    ///
    /// Args:
    ///     header (bytes): The sequence identifier.
    ///     sequence (bytes): The raw nucleotide bytes.
    ///
    /// Returns:
    ///     List[Gene]: A list of predicted Gene objects.
    fn find_genes(&self, py: Python, header: &[u8], sequence: &[u8]) -> PyResult<Vec<Gene>> {
        // 1. Build the Nuc vector directly from Python's memory.
        let nuc_seq: Vec<Nuc> = sequence
            .iter()
            .filter(|&&b| !b.is_ascii_whitespace())
            .map(|&b| b.to_ascii_uppercase())
            .map(Nuc::from)
            .collect();

        let head_bytes = header.to_vec();

        // 2. Run the heavy HMM inference without the GIL
        let prediction = py.allow_threads(|| {
            let cg = count_cg_content(&nuc_seq);
            let local_model = &self.locals[cg];

            viterbi(
                &self.global,
                local_model,
                head_bytes,
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
}

/// Represents a frameshift mutation detected by the ab initio model.
#[pyclass(get_all)]
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

#[pymethods]
impl Mutation {
    fn __repr__(&self) -> String {
        format!(
            "<Mutation type='{}' pos={} codon={} mut='{}'>",
            self.mut_type, self.pos, self.codon_idx, self.annotation
        )
    }
}

/// A streaming writer for Extended BED (BED6+1) files.
#[pyclass]
pub struct BedWriter {
    writer: BufWriter<File>,
}

#[pymethods]
impl BedWriter {
    #[new]
    fn new(output_path: &str) -> PyResult<Self> {
        let file = File::create(output_path)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("Could not create BED: {}", e)))?;
        let mut writer = BufWriter::new(file);

        let version = env!("CARGO_PKG_VERSION");
        writeln!(writer, "# source=ab initio prediction: pyfgs v{}", version)?;

        // Added the 7th column for MUTATIONS
        writeln!(writer, "# CHROM\tSTART\tEND\tNAME\tSCORE\tSTRAND\tMUTATIONS")?;

        Ok(Self { writer })
    }

    fn __enter__(slf: Py<Self>) -> Py<Self> {
        slf
    }

    fn __exit__(
        &mut self,
        _exc_type: Option<Bound<'_, PyAny>>,
        _exc_value: Option<Bound<'_, PyAny>>,
        _traceback: Option<Bound<'_, PyAny>>,
    ) -> PyResult<()> {
        self.writer.flush()?;
        Ok(())
    }

    // Now requires the sequence byte slice to extract mutations
    fn write_record(&mut self, genes: Vec<PyRef<Gene>>, header: &str, sequence: &[u8]) -> PyResult<()> {
        for (g_idx, gene) in genes.iter().enumerate() {
            let gene_id = format!("{}_FGS_{}", header, g_idx + 1);
            let strand_char = if gene.strand == 1 { "+" } else { "-" };

            // Query the struct for its mutations
            let mutations = gene.extract_mutations(sequence);

            // Format the 7th column into a parseable key-value list
            let mut_info = if mutations.is_empty() {
                ".".to_string() // Standard VCF/BED placeholder for 'empty'
            } else {
                let formatted_muts: Vec<String> = mutations.iter().map(|m| {
                    format!("pos={};type={};codon={};note={}",
                        m.pos, m.mut_type, m.codon_idx, m.annotation)
                }).collect();

                // Join multiple mutations with a comma so they stay safely in column 7
                formatted_muts.join(",")
            };

            writeln!(
                self.writer,
                "{}\t{}\t{}\t{}\t{:.2}\t{}\t{}",
                header, gene.start, gene.end, gene_id, gene.score, strand_char, mut_info
            )?;
        }
        Ok(())
    }
}

/// A streaming writer for VCF v4.2 files.
#[pyclass]
pub struct VcfWriter {
    writer: BufWriter<File>,
}

#[pymethods]
impl VcfWriter {
    /// Opens the file and writes the VCF headers.
    #[new]
    fn new(output_path: &str) -> PyResult<Self> {
        let file = File::create(output_path)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("Could not create VCF: {}", e)))?;
        let mut writer = BufWriter::new(file);

        writeln!(writer, "##fileformat=VCFv4.2")?;
        writeln!(writer, "##source=pyfgs_ab_initio")?;
        writeln!(writer, "##INFO=<ID=TYPE,Number=1,Type=String,Description=\"Type of sequence discrepancy\">")?;
        writeln!(writer, "##INFO=<ID=GENE,Number=1,Type=String,Description=\"Predicted Gene ID\">")?;
        writeln!(writer, "##INFO=<ID=CODON,Number=1,Type=Integer,Description=\"1-based codon index of the frameshift\">")?;
        writeln!(writer, "##INFO=<ID=FGS_MUT,Number=1,Type=String,Description=\"Snippy-style AA and DNA translation shift\">")?;
        writeln!(writer, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")?;

        Ok(Self { writer })
    }

    /// Context manager entry (`__enter__`)
    fn __enter__(slf: Py<Self>) -> Py<Self> {
        slf
    }

    /// Context manager exit (`__exit__`). Flushes and closes the file.
    fn __exit__(
        &mut self,
        _exc_type: Option<Bound<'_, PyAny>>,
        _exc_value: Option<Bound<'_, PyAny>>,
        _traceback: Option<Bound<'_, PyAny>>,
    ) -> PyResult<()> {
        self.writer.flush()?;
        Ok(())
    }

    /// Writes the frameshift variants for a single contig/sequence.
    fn write_record(&mut self, genes: Vec<PyRef<Gene>>, header: &str, sequence: &[u8]) -> PyResult<()> {
        for (g_idx, gene) in genes.iter().enumerate() {
            let mutations = gene.extract_mutations(sequence);

            for muta in mutations {
                let type_name = if muta.mut_type == "ins" { "insertion" } else { "deletion" };

                // We construct the SnpEff ANN field directly in the writeln! macro.
                // Format: Allele | Annotation | Impact | Gene_Name | Gene_ID | Feature_Type | ...

                writeln!(
                    self.writer,
                    "{chr}\t{pos}\t.\t{ref_all}\t{alt_all}\t.\tPASS\tTYPE=frameshift_{type};GENE={chr}_FGS_{idx};CODON={codon};ANN={alt_all}|frameshift_variant|HIGH|{chr}_FGS_{idx}|{chr}_FGS_{idx}|transcript|{chr}_FGS_{idx}|protein_coding|{codon}/...|{muta_note}",
                    chr = header,
                    pos = muta.pos,
                    ref_all = muta.ref_allele,
                    alt_all = muta.alt_allele,
                    type = type_name,
                    idx = g_idx + 1,
                    codon = muta.codon_idx,
                    muta_note = muta.annotation
                )?;
            }
        }
        Ok(())
    }
}

/// A streaming writer for INSDC-compliant GFF3 files.
#[pyclass]
pub struct Gff3Writer {
    writer: BufWriter<File>,
}

#[pymethods]
impl Gff3Writer {
    #[new]
    fn new(output_path: &str) -> PyResult<Self> {
        let file = File::create(output_path)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("Could not create GFF3: {}", e)))?;
        let mut writer = BufWriter::new(file);

        let version = env!("CARGO_PKG_VERSION");
        writeln!(writer, "##gff-version 3")?;
        writeln!(writer, "#source-ontology=ab initio prediction: pyfgs v{}", version)?;

        Ok(Self { writer })
    }

    fn __enter__(slf: Py<Self>) -> Py<Self> {
        slf
    }

    fn __exit__(
        &mut self,
        _exc_type: Option<Bound<'_, PyAny>>,
        _exc_value: Option<Bound<'_, PyAny>>,
        _traceback: Option<Bound<'_, PyAny>>,
    ) -> PyResult<()> {
        self.writer.flush()?;
        Ok(())
    }

    fn write_record(&mut self, genes: Vec<PyRef<Gene>>, header: &str, sequence: &[u8]) -> PyResult<()> {
        let version = env!("CARGO_PKG_VERSION");
        let source_name = format!("pyfgs_v{}", version);

        for (g_idx, gene) in genes.iter().enumerate() {
            let gene_id = format!("{}_FGS_{}", header, g_idx + 1);

            // GFF3 requires 1-based fully-closed coordinates
            let gff_start = gene.start + 1;
            let gff_end = gene.end;
            let strand_char = if gene.strand == 1 { "+" } else { "-" };

            let mutations = gene.extract_mutations(sequence);

            if mutations.is_empty() {
                // Intact Gene -> Standard CDS
                writeln!(
                    self.writer,
                    "{chr}\t{src}\tCDS\t{start}\t{end}\t{score:.2}\t{strand}\t0\tID={id};inference=ab initio prediction:{src}",
                    chr = header,
                    src = source_name,
                    start = gff_start,
                    end = gff_end,
                    score = gene.score,
                    strand = strand_char,
                    id = gene_id
                )?;
            } else {
                // Frameshifted Gene -> INSDC Pseudogene Flagging
                // We format a human-readable note for the researchers reviewing the file
                let notes: Vec<String> = mutations.iter().map(|m| {
                    let mut_name = if m.mut_type == "ins" { "insertion" } else { "deletion" };
                    format!("Frameshift {} at pos {} (codon {})", mut_name, m.pos, m.codon_idx)
                }).collect();

                let note_str = notes.join(",");

                writeln!(
                    self.writer,
                    "{chr}\t{src}\tpseudogene\t{start}\t{end}\t{score:.2}\t{strand}\t0\tID={id};inference=ab initio prediction:{src};pseudogene=unknown;Note={note}",
                    chr = header,
                    src = source_name,
                    start = gff_start,
                    end = gff_end,
                    score = gene.score,
                    strand = strand_char,
                    id = gene_id,
                    note = note_str
                )?;
            }
        }
        Ok(())
    }
}

/// A streaming writer for non-wrapped nucleotide FASTA (.fna) files.
#[pyclass]
pub struct FnaWriter {
    writer: BufWriter<File>,
}

#[pymethods]
impl FnaWriter {
    #[new]
    fn new(output_path: &str) -> PyResult<Self> {
        let file = File::create(output_path)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("Could not create FNA: {}", e)))?;
        Ok(Self { writer: BufWriter::new(file) })
    }

    fn __enter__(slf: Py<Self>) -> Py<Self> { slf }

    fn __exit__(&mut self, _exc_type: Option<Bound<'_, PyAny>>, _exc_value: Option<Bound<'_, PyAny>>, _traceback: Option<Bound<'_, PyAny>>) -> PyResult<()> {
        self.writer.flush()?;
        Ok(())
    }

    fn write_record(&mut self, genes: Vec<PyRef<Gene>>, header: &str) -> PyResult<()> {
        for (g_idx, gene) in genes.iter().enumerate() {
            let strand_char = if gene.strand == 1 { "+" } else { "-" };

            // Header format: >ID contig:start-end strand=+
            writeln!(
                self.writer,
                ">{}_FGS_{} {}:{}-{} strand={}",
                header, g_idx + 1, header, gene.start + 1, gene.end, strand_char
            )?;

            // Write the raw bytes directly without any string formatting or wrapping!
            self.writer.write_all(&gene.get_dna_bytes())?;
            writeln!(self.writer)?; // Add the final newline
        }
        Ok(())
    }
}

/// A streaming writer for non-wrapped amino acid FASTA (.faa) files.
#[pyclass]
pub struct FaaWriter {
    writer: BufWriter<File>,
}

#[pymethods]
impl FaaWriter {
    #[new]
    fn new(output_path: &str) -> PyResult<Self> {
        let file = File::create(output_path)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("Could not create FAA: {}", e)))?;
        Ok(Self { writer: BufWriter::new(file) })
    }

    fn __enter__(slf: Py<Self>) -> Py<Self> { slf }

    fn __exit__(&mut self, _exc_type: Option<Bound<'_, PyAny>>, _exc_value: Option<Bound<'_, PyAny>>, _traceback: Option<Bound<'_, PyAny>>) -> PyResult<()> {
        self.writer.flush()?;
        Ok(())
    }

    fn write_record(&mut self, genes: Vec<PyRef<Gene>>, header: &str) -> PyResult<()> {
        for (g_idx, gene) in genes.iter().enumerate() {
            let strand_char = if gene.strand == 1 { "+" } else { "-" };

            writeln!(
                self.writer,
                ">{}_FGS_{} {}:{}-{} strand={}",
                header, g_idx + 1, header, gene.start + 1, gene.end, strand_char
            )?;

            self.writer.write_all(&gene.get_aa_bytes())?;
            writeln!(self.writer)?;
        }
        Ok(())
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
    m.add_class::<Mutation>()?;
    m.add_class::<BedWriter>()?;
    m.add_class::<VcfWriter>()?;
    m.add_class::<Gff3Writer>()?;
    m.add_class::<FnaWriter>()?;
    m.add_class::<FaaWriter>()?;
    Ok(())
}