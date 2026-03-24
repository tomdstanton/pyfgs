# 🔗🐍⏭️ `pyfgs` [![Stars](https://img.shields.io/github/stars/tomdstanton/pyfgs.svg?style=social&maxAge=3600&label=Star)](https://github.com/tomdstanton/pyfgs/stargazers)

*PyO3 bindings and Python interface to [FragGeneScanRs](https://github.com/unipept/FragGeneScanRs),
a gene prediction model for short and error-prone reads.*

[![Release](https://img.shields.io/github/v/release/tomdstanton/pyfgs?style=flat-square)](https://img.shields.io/github/v/release/tomdstanton/pyfgs)
[![License](https://img.shields.io/github/license/tomdstanton/pyfgs?style=flat-square)](https://img.shields.io/github/license/tomdstanton/pyfgs)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19059429.svg?style=flat-square)](https://doi.org/10.5281/zenodo.19059429)
[![PyPI](https://img.shields.io/pypi/v/pyfgs.svg?style=flat-square&maxAge=3600&logo=PyPI)](https://pypi.org/project/pyfgs)
[![Wheel](https://img.shields.io/pypi/wheel/pyfgs.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/pyfgs/#files)
[![Python Versions](https://img.shields.io/pypi/pyversions/pyfgs.svg?style=flat-square&maxAge=600&logo=python)](https://pypi.org/project/pyfgs/#files)
[![Python Implementations](https://img.shields.io/pypi/implementation/pyfgs.svg?style=flat-square&maxAge=600&label=impl)](https://pypi.org/project/pyfgs/#files)
[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/tomdstanton/pyfgs/)
[![Issues](https://img.shields.io/github/issues/tomdstanton/pyfgs.svg?style=flat-square&maxAge=600)](https://github.com/tomdstanton/pyfgs/issues)
[![Docs](https://img.shields.io/badge/Material_for_MkDocs-526CFE??style=flat-square&maxAge=600)](https://tomdstanton.github.io/pyfgs/)
[![Changelog](https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square)](https://github.com/tomdstanton/pyfgs/blob/main/CHANGELOG.md)
[![Downloads](https://img.shields.io/pypi/dm/pyfgs?style=flat-square&color=303f9f&maxAge=86400&label=downloads)](https://pepy.tech/project/pyfgs)


## 🗺️ Why `pyfgs`?

**Built for noisy data**

Standard ab initio predictors (like Prodigal or Pyrodigal) are fantastic for pristine, fully assembled contigs.
However, they struggle with raw metagenomic reads or error-prone assemblies because they immediately break the open
reading frame at the first sign of an indel. `pyfgs` uses an error-tolerant Hidden Markov Model trained on specific
sequencing profiles (Illumina, 454, Sanger) to power through these sequencing errors, dynamically correct the reading
frame, and salvage the translated protein.

**Native frameshift tracking**

Instead of just silently stitching broken genes together, `pyfgs` exposes the exact coordinates of every hallucinated or
skipped base directly to Python. This allows you to rigorously track structural variants, correctly annotate
INSDC-compliant pseudogenes, or export exact frameshift coordinates for downstream quality control.

**No subprocess I/O tax**

Running standard CLI bioinformatics tools from Python usually requires a heavy I/O penalty: dumping sequences to a
temporary FASTA file, firing a subprocess, and parsing the text outputs back into memory. `pyfgs` binds directly to the
underlying Rust engine. The HMM runs entirely in memory and yields native Python objects ready for immediate analysis.

**True multithreading and zero-copy memory**

`pyfgs` is designed to process massive datasets efficiently:

- GIL-Free Inference: The Rust backend completely releases the Python Global Interpreter Lock (GIL) during the heavy
  HMM math. You can drop the predictor into a standard ThreadPoolExecutor and achieve true parallel processing across
  all your CPU cores.

- Zero-Copy Bytes: The engine borrows raw byte slices `(&[u8])` directly from Python's memory, bypassing the overhead of
  copying strings between languages.

- Lazy Translation: Translating DNA to amino acids is computationally expensive. `pyfgs` evaluates sequence strings
  lazily, meaning you only pay the CPU and memory cost of string allocation if you explicitly request the sequence data.

**A Pythonic API**

Bioinformatics coordinates are notoriously messy. `pyfgs` outputs standard 0-based, half-open intervals ([start, end)),
allowing you to slice sequence arrays immediately without wrestling with 1-based GFF3 coordinate math. When you do need
standardized files, it includes heavily optimized, native-Rust context managers to stream perfectly compliant VCF, BED,
GFF3, and FASTA files directly to disk without bloating your RAM.


## 🔧 Installing

This project is supported on Python 3.10 and later.

`pyfgs` can be installed directly from [PyPI](https://pypi.org/project/pyfgs/):

```console
pip install pyfgs
```

## 💻 Usage

### API Usage
For full API usage, please refer to the [documentation](https://tomdstanton.github.io/pyfgs/api/).

```python
import concurrent.futures
import pyfgs
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO

def process_contig(header_bytes: bytes, seq_bytes: bytes, finder: pyfgs.GeneFinder):
    """
    Worker function to find genes.
    Because find_genes() drops the GIL, this runs in true parallel.
    """
    genes = finder.find_genes(header_bytes, seq_bytes)
    return header_bytes, seq_bytes, genes

def main():
    # 1. Initialize the GeneFinder
    # We use the 'Complete' model for high-quality assemblies, but strictly
    # set whole_genome=False to force the HMM to hunt for frameshifts.
    finder = pyfgs.GeneFinder(pyfgs.Model.Complete, whole_genome=False)

    # 2. Stream the genome using the zero-allocation Rust FastaReader
    reader = pyfgs.FastaReader("bacterial_assembly.fasta")

    results = []

    # 3. Process all contigs concurrently
    with concurrent.futures.ThreadPoolExecutor() as executor:
        # Submit all contigs to the thread pool
        futures = [
            executor.submit(process_contig, header, seq, finder)
            for header, seq in reader
        ]

        # Gather results as they complete
        for future in concurrent.futures.as_completed(futures):
            results.append(future.result())

    # 4. Format into INSDC-compliant GenBank records
    records = []
    for header_bytes, seq_bytes, genes in results:
        header_str = header_bytes.decode('utf-8')
        record = SeqRecord(
            Seq(seq_bytes.decode('utf-8')),
            id=header_str,
            name=header_str,
            description="Annotated by pyfgs"
        )

        for i, gene in enumerate(genes):
            # Query the Rust backend for structural variants
            mutations = gene.mutations(seq_bytes)

            # INSDC Standard: Frameshifted ORFs cannot be 'CDS', must be 'pseudogene'
            feature_type = "pseudogene" if mutations else "CDS"

            qualifiers = {
                "source": "pyfgs",
                "inference": "ab initio prediction:pyfgs",
                "ID": f"{header_str}_FGS_{i+1}"
            }

            if mutations:
                qualifiers["pseudogene"] = ["unknown"]
                notes = []
                for mut in mutations:
                    mut_name = "insertion" if mut.mut_type == "ins" else "deletion"
                    # Include our Snippy-style variant notation
                    notes.append(f"Frameshift {mut_name} at pos {mut.pos} (codon {mut.codon_idx}). {mut.annotation}")
                qualifiers["note"] = notes
            else:
                # Only strictly intact CDS features receive a translation qualifier
                qualifiers["translation"] = [gene.translation().decode('utf-8')]

            # Biopython's FeatureLocation is natively 0-based and half-open,
            # mapping perfectly to our Gene.start and Gene.end!
            location = FeatureLocation(gene.start, gene.end, strand=gene.strand)
            feature = SeqFeature(location=location, type=feature_type, qualifiers=qualifiers)
            record.features.append(feature)

        records.append(record)

    # 5. Export to GenBank
    output_file = "annotated_genome.gbk"
    with open(output_file, "w") as out_handle:
        SeqIO.write(records, out_handle, "genbank")

if __name__ == "__main__":
    main()
```

### CLI Usage

For CLI usage, type `pyfgs --help`

```console
usage: pyfgs <seq> [options]

🔗🐍⏭️	PyO3 bindings and Python interface to FragGeneScanRs,
	a gene prediction model for short and error-prone reads.

Input options 💽:

  seq                 Sequence file (or '-' for stdin)
  -m, --model         Sequence error model (default: complete)
                       - short1: Illumina sequencing reads with about 0.1% error rate
                       - short5: Illumina sequencing reads with about 0.5% error rate
                       - short10: Illumina sequencing reads with about 1% error rate
                       - sanger5: Sanger sequencing reads with about 0.5% error rate
                       - sanger10: Sanger sequencing reads with about 1% error rate
                       - pyro5: 454 pyrosequencing reads with about 0.5% error rate
                       - pyro10: 454 pyrosequencing reads with about 1% error rate
                       - pyro30: 454 pyrosequencing reads with about 3% error rate
                       - complete: Complete genomic sequences or short sequence reads without sequencing error
  -r, --reads         Force FASTQ parsing (Overrides auto-detection)
  -w, --whole-genome  Strict contiguous ORFs. Disables error-tolerant frameshift detection.

Output options ⚙️:
  Provide a PATH to save to a file, or use the flag alone to print to stdout.

  --faa [PATH]        Output protein FASTA
  --fna [PATH]        Output nucleotide FASTA
  --bed [PATH]        Output BED6+1 format
  --gff [PATH]        Output GFF3 format
  --vcf [PATH]        Output VCF v4.2 format

Other options 🚧:

  -t, --threads       Number of threads (default: optimal)
  -v, --version       Print version and exit
  -h, --help          Print help and exit
```


## 🔖 Citation

For now, please cite the original
[FragGeneScanRs paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04736-5):

> Van der Jeugt, F., Dawyndt, P. & Mesuere, B. FragGeneScanRs: faster gene prediction for short reads.
BMC Bioinformatics 23, 198 (2022). https://doi.org/10.1186/s12859-022-04736-5


## 💭 Feedback

### ⚠️ Issue Tracker

Found a bug ? Have an enhancement request ? Head over to the
[GitHub issue tracker](https://github.com/tomdstanton/pyfgs/issues) if you need to report
or ask something. If you are filing in on a bug, please include as much
information as you can about the issue, and try to recreate the same bug
in a simple, easily reproducible situation.

### 🏗️ Contributing

Contributions are more than welcome! See
[`CONTRIBUTING.md`](https://github.com/tomdstanton/pyfgs/blob/main/CONTRIBUTING.md)
for more details.

## 📋 Changelog

This project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html)
and provides a [changelog](https://github.com/tomdstanton/pyfgs/blob/main/CHANGELOG.md)
in the [Keep a Changelog](http://keepachangelog.com/en/1.0.0/) format.


## ⚖️ License

This library is provided under the [GNU General Public License v3.0](https://choosealicense.com/licenses/gpl-3.0/).
The FragGeneScanRs code was written by [Peter Dawyndt](https://github.com/pdawyndt),
[Bart Mesuere](https://github.com/bmesuere) and
[Felix Van der Jeugt](https://github.com/ninewise) and is distributed under the
terms of the GPLv3 as well. See `https://github.com/FragGeneScanRs/LICENSE` for more information.

*This project is in no way affiliated, sponsored, or otherwise endorsed
by the original FragGeneScanRs authors [Peter Dawyndt](https://github.com/pdawyndt),
[Bart Mesuere](https://github.com/bmesuere) and
[Felix Van der Jeugt](https://github.com/ninewise). It was developed
by [Tom Stanton](https://github.com/tomdstanton/) during his Post-doc project
at [Monash University](https://www.monash.edu/medicine/translational/infectious-diseases) in
the [Wryes Lab](https://wyreslab.com/).*
