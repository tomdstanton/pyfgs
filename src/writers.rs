use pyo3::prelude::*;
use pyo3::types::PyBytes;
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

use std::fs::File;
use std::io::{BufWriter, Write};

use crate::models::Gene; // Needed to reference Gene in writers

/// A streaming writer for Extended BED (BED6+1) files.
#[gen_stub_pyclass]
#[pyclass]
pub struct BedWriter {
    writer: BufWriter<File>,
}

#[gen_stub_pymethods]
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
    fn write_record(&mut self, genes: Vec<PyRef<Gene>>, header: &str, sequence: &Bound<'_, PyBytes>) -> PyResult<()> {
        let sequence = sequence.as_bytes();
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
#[gen_stub_pyclass]
#[pyclass]
pub struct VcfWriter {
    writer: BufWriter<File>,
}

#[gen_stub_pymethods]
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
    fn write_record(&mut self, genes: Vec<PyRef<Gene>>, header: &str, sequence: &Bound<'_, PyBytes>) -> PyResult<()> {
        let sequence = sequence.as_bytes();
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
#[gen_stub_pyclass]
#[pyclass]
pub struct Gff3Writer {
    writer: BufWriter<File>,
}

#[gen_stub_pymethods]
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

    fn write_record(&mut self, genes: Vec<PyRef<Gene>>, header: &str, sequence: &Bound<'_, PyBytes>) -> PyResult<()> {
        let sequence = sequence.as_bytes();
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
#[gen_stub_pyclass]
#[pyclass]
pub struct FnaWriter {
    writer: BufWriter<File>,
}

#[gen_stub_pymethods]
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
#[gen_stub_pyclass]
#[pyclass]
pub struct FaaWriter {
    writer: BufWriter<File>,
}

#[gen_stub_pymethods]
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
