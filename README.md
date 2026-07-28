# 🔗🐍⏭️ `pyfgs` [![Stars](https://img.shields.io/github/stars/tomdstanton/pyfgs.svg?style=social&maxAge=3600&label=Star)](https://github.com/tomdstanton/pyfgs/stargazers)

*Python bindings to [FragGeneScanRs](https://github.com/unipept/FragGeneScanRs) for fast, error-tolerant bacterial ORF prediction*

[![Release](https://img.shields.io/github/v/release/tomdstanton/pyfgs?style=flat-square)](https://img.shields.io/github/v/release/tomdstanton/pyfgs)
[![License](https://img.shields.io/github/license/tomdstanton/pyfgs?style=flat-square)](https://img.shields.io/github/license/tomdstanton/pyfgs)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19059428.svg?style=flat-square)](https://doi.org/10.5281/zenodo.19059428)
[![PyPI](https://img.shields.io/pypi/v/pyfgs.svg?style=flat-square&maxAge=3600&logo=PyPI)](https://pypi.org/project/pyfgs)
[![Wheel](https://img.shields.io/pypi/wheel/pyfgs.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/pyfgs/#files)
[![Python Versions](https://img.shields.io/pypi/pyversions/pyfgs.svg?style=flat-square&maxAge=600&logo=python)](https://pypi.org/project/pyfgs/#files)
[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/tomdstanton/pyfgs/)
[![Issues](https://img.shields.io/github/issues/tomdstanton/pyfgs.svg?style=flat-square&maxAge=600)](https://github.com/tomdstanton/pyfgs/issues)
[![Downloads](https://img.shields.io/pypi/dm/pyfgs?style=flat-square&color=303f9f&maxAge=86400&label=downloads)](https://pepy.tech/project/pyfgs)

**For full usage, see the [docs](https://tomdstanton.github.io/pyfgs)**

##  Why `pyfgs`?

**Built for noisy data**

Standard ab initio predictors are fantastic for pristine, fully assembled contigs.
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

**Blazing Fast**

Thanks to the memory-safe zero-copy architecture and multi-threading through `Rayon`, `pyfgs` can outperform even the fastest SIMD-optimized Gene Finders (e.g. `pyrodigal`) in Whole Genome mode—especially on highly fragmented assemblies and metagenomes.


## 🚀 Usage

### 1. Vectorized Batch Prediction

If you are a power user integrating gene prediction into a larger Python pipeline, `pyfgs` offers a vectorized, zero-copy NumPy API that sidesteps Python object overhead completely.

```python
import pyfgs

# 1. Initialize the gene finder
# By default, this loads the Illumina HMM model (error-tolerant).
finder = pyfgs.GeneFinder()

# 2. Extract byte strings of sequence
# Notice we pass bytes (b""), not strings, to avoid encoding overhead.
sequences = [
    b"ATGCCCGGGAAATTTTGACCC",
    b"ATGAAAAAA",
]

# 3. Predict genes across all sequences in parallel (releasing the GIL)
batch = finder.find_genes_batch(sequences)

# 4. Access predictions as flat, 1-dimensional NumPy arrays (Structure-of-Arrays)
print(batch.starts)  # array([0, ...])
print(batch.ends)  # array([18, ...])
print(batch.strands)  # array([1, ...])
print(batch.scores)  # array([12.34, ...])

# 5. Fields with variable counts per gene (e.g. insertions/deletions)
# are exposed as flat arrays alongside parallel "offsets" arrays
# defining the slice boundaries (i.e. ragged arrays).
print(batch.insertions_offsets)
print(batch.insertions_flat)
```

### 2. High-Performance I/O Pipeline

If your goal is just to read a FASTA file and output standard format files (GFF3, BED, etc), `pyfgs` provides a pure-Rust IO pipeline that bypasses Python entirely. This is what the `pyfgs` CLI uses under the hood.

```python
import pyfgs

finder = pyfgs.GeneFinder()

# Stream directly from disk to disk using Rust native threading
finder.run_file(
    input_path="input.fasta",
    is_fastq=False,
    outputs={
        "gff3": "output.gff3",
        "faa": "output.faa",
    },
)
```

## 🔧 Installing

This project is supported on Python 3.10 and later.

`pyfgs` can be installed directly from [PyPI](https://pypi.org/project/pyfgs/):

```bash
pip install pyfgs
```

⚡️ Power users ⚡️ can force your local machine to compile the Rust engine specifically for your own CPU by running:

```bash
RUSTFLAGS="-C target-cpu=native" pip install --no-binary pyfgs pyfgs
```

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
