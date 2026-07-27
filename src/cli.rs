use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Write};
use rayon::prelude::*;
use seq_io::fasta::{Reader as FastaReader, Record as FastaRecordTrait};
use seq_io::fastq::{Reader as FastqReader, Record as FastqRecordTrait};
use pyo3::prelude::*;

use frag_gene_scan_rs::dna::{count_cg_content, Nuc};
use frag_gene_scan_rs::viterbi::viterbi;

use crate::models::{GeneFinder, Gene};

struct CliWriter {
    writer: Box<dyn Write>,
    fmt: String,
}

impl CliWriter {
    fn new(path: &str, fmt: &str) -> io::Result<Self> {
        let writer: Box<dyn Write> = if path == "-" {
            Box::new(BufWriter::new(io::stdout()))
        } else {
            Box::new(BufWriter::new(File::create(path)?))
        };
        
        let mut cli_writer = Self {
            writer,
            fmt: fmt.to_string(),
        };

        // Write headers
        let version = env!("CARGO_PKG_VERSION");
        match fmt {
            "bed" => {
                writeln!(cli_writer.writer, "# source=ab initio prediction: pyfgs v{}", version)?;
                writeln!(cli_writer.writer, "# CHROM\tSTART\tEND\tNAME\tSCORE\tSTRAND\tMUTATIONS")?;
            },
            "gff" => {
                writeln!(cli_writer.writer, "##gff-version 3")?;
                writeln!(cli_writer.writer, "#source-ontology=ab initio prediction: pyfgs v{}", version)?;
            },
            "vcf" => {
                writeln!(cli_writer.writer, "##fileformat=VCFv4.2")?;
                writeln!(cli_writer.writer, "##source=pyfgs_ab_initio")?;
                writeln!(cli_writer.writer, "##INFO=<ID=TYPE,Number=1,Type=String,Description=\"Type of sequence discrepancy\">")?;
                writeln!(cli_writer.writer, "##INFO=<ID=GENE,Number=1,Type=String,Description=\"Predicted Gene ID\">")?;
                writeln!(cli_writer.writer, "##INFO=<ID=CODON,Number=1,Type=Integer,Description=\"1-based codon index of the frameshift\">")?;
                writeln!(cli_writer.writer, "##INFO=<ID=ANN,Number=.,Type=String,Description=\"Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_Biotype | Rank | HGVS.c'\">")?;
                writeln!(cli_writer.writer, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")?;
            },
            _ => {}
        }
        
        Ok(cli_writer)
    }

    fn write_batch(&mut self, headers: &[String], seqs: &[Vec<u8>], all_genes: &[Vec<Gene>]) -> io::Result<()> {
        let version = env!("CARGO_PKG_VERSION");
        let source_name = format!("pyfgs_v{}", version);

        for ((header, sequence), genes) in headers.iter().zip(seqs.iter()).zip(all_genes.iter()) {
            if genes.is_empty() {
                continue;
            }

            for (g_idx, gene) in genes.iter().enumerate() {
                let gene_id = format!("{}_FGS_{}", header, g_idx + 1);
                let strand_char = if gene.strand == 1 { "+" } else { "-" };

                match self.fmt.as_str() {
                    "faa" => {
                        writeln!(
                            self.writer,
                            ">{gene_id} {chr}:{start}-{end} strand={strand}",
                            chr = header, start = gene.start + 1, end = gene.end, strand = strand_char
                        )?;
                        self.writer.write_all(&gene.get_aa_bytes())?;
                        writeln!(self.writer)?;
                    }
                    "fna" => {
                        writeln!(
                            self.writer,
                            ">{gene_id} {chr}:{start}-{end} strand={strand}",
                            chr = header, start = gene.start + 1, end = gene.end, strand = strand_char
                        )?;
                        self.writer.write_all(&gene.get_dna_bytes())?;
                        writeln!(self.writer)?;
                    }
                    "bed" => {
                        let mutations = gene.extract_mutations(sequence);
                        let mut_info = if mutations.is_empty() {
                            ".".to_string()
                        } else {
                            let formatted_muts: Vec<String> = mutations.iter().map(|m| {
                                format!("pos={};type={};codon={};note={}", m.pos, m.mut_type, m.codon_idx, m.annotation)
                            }).collect();
                            formatted_muts.join(",")
                        };

                        writeln!(
                            self.writer,
                            "{chr}\t{start}\t{end}\t{id}\t{score:.2}\t{strand}\t{mut_info}",
                            chr = header, start = gene.start, end = gene.end, id = gene_id,
                            score = gene.score, strand = strand_char, mut_info = mut_info
                        )?;
                    }
                    "gff" => {
                        let mutations = gene.extract_mutations(sequence);
                        let gff_start = gene.start + 1;
                        
                        if mutations.is_empty() {
                            writeln!(
                                self.writer,
                                "{chr}\t{src}\tCDS\t{start}\t{end}\t{score:.2}\t{strand}\t0\tID={id};inference=ab initio prediction:{src}",
                                chr = header, src = source_name, start = gff_start, end = gene.end,
                                score = gene.score, strand = strand_char, id = gene_id
                            )?;
                        } else {
                            let notes: Vec<String> = mutations.iter().map(|m| {
                                let mut_name = if m.mut_type == "ins" { "insertion" } else { "deletion" };
                                format!("Frameshift {} at pos {} (codon {})", mut_name, m.pos, m.codon_idx)
                            }).collect();
                            let note_str = notes.join(",");

                            writeln!(
                                self.writer,
                                "{chr}\t{src}\tpseudogene\t{start}\t{end}\t{score:.2}\t{strand}\t0\tID={id};inference=ab initio prediction:{src};pseudogene=unknown;Note={note}",
                                chr = header, src = source_name, start = gff_start, end = gene.end,
                                score = gene.score, strand = strand_char, id = gene_id, note = note_str
                            )?;
                        }
                    }
                    "vcf" => {
                        let mutations = gene.extract_mutations(sequence);
                        for muta in mutations {
                            let type_name = if muta.mut_type == "ins" { "insertion" } else { "deletion" };
                            writeln!(
                                self.writer,
                                "{chr}\t{pos}\t.\t{ref_all}\t{alt_all}\t.\tPASS\tTYPE=frameshift_{type};GENE={id};CODON={codon};ANN={alt_all}|frameshift_variant|HIGH|{id}|{id}|transcript|{id}|protein_coding|{codon}/...|{muta_note}",
                                chr = header, pos = muta.pos, ref_all = muta.ref_allele, alt_all = muta.alt_allele,
                                type = type_name, id = gene_id, codon = muta.codon_idx, muta_note = muta.annotation
                            )?;
                        }
                    }
                    _ => {}
                }
            }
        }
        self.writer.flush()?;
        Ok(())
    }
}

fn process_chunk(
    finder: &GeneFinder,
    seqs: &[Vec<u8>],
) -> Vec<Vec<Gene>> {
    seqs.par_iter().map(|sequence| {
        let nuc_seq: Vec<Nuc> = sequence
            .iter()
            .filter(|&&b| !b.is_ascii_whitespace())
            .map(|&b| b.to_ascii_uppercase())
            .map(Nuc::from)
            .collect();

        let cg = count_cg_content(&nuc_seq);
        let local_model = &finder.locals[cg];

        let prediction = viterbi(
            &finder.global,
            local_model,
            Vec::new(),
            nuc_seq,
            finder.whole_genome
        );

        prediction.genes.into_iter().map(|g| {
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
                whole_genome: finder.whole_genome,
            }
        }).collect()
    }).collect()
}

pub fn run_cli_pipeline(
    py: Python<'_>,
    finder: &GeneFinder,
    input_path: &str,
    is_fastq: bool,
    outputs: HashMap<String, String>,
) -> PyResult<()> {
    let mut writers = Vec::new();
    for (fmt, path) in outputs {
        writers.push(CliWriter::new(&path, &fmt).map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?);
    }
    
    let mut headers = Vec::with_capacity(10000);
    let mut seqs = Vec::with_capacity(10000);
    
    let mut process_and_write = |headers: &mut Vec<String>, seqs: &mut Vec<Vec<u8>>| -> PyResult<()> {
        if headers.is_empty() { return Ok(()); }
        
        let genes = py.detach(|| process_chunk(finder, seqs));
        
        for writer in &mut writers {
            writer.write_batch(headers, seqs, &genes).map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
        }
        headers.clear();
        seqs.clear();
        Ok(())
    };

    macro_rules! process_reader {
        ($parser_type:ident, $reader:expr) => {
            let mut parser = $parser_type::new($reader);
            while let Some(record) = parser.next() {
                let rec = record.map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
                headers.push(String::from_utf8_lossy(rec.id_bytes()).split_whitespace().next().unwrap_or("").to_string());
                seqs.push(rec.seq().to_vec());
                if headers.len() >= 10000 {
                    process_and_write(&mut headers, &mut seqs)?;
                }
            }
        };
    }

    if input_path == "-" {
        let stdin = io::stdin();
        let reader = BufReader::new(stdin);
        if is_fastq {
            process_reader!(FastqReader, reader);
        } else {
            process_reader!(FastaReader, reader);
        }
    } else {
        let file = File::open(input_path).map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
        let reader = BufReader::new(file);
        if is_fastq {
            process_reader!(FastqReader, reader);
        } else {
            process_reader!(FastaReader, reader);
        }
    }
    
    // Flush remaining
    process_and_write(&mut headers, &mut seqs)?;
    
    Ok(())
}
