# 🔗🐍⏭️ `pyfgs` [![Stars](https://img.shields.io/github/stars/tomdstanton/pyfgs.svg?style=social&maxAge=3600&label=Star)](https://github.com/tomdstanton/pyfgs/stargazers)

*PyO3 bindings and Python interface to [FragGeneScanRs](https://github.com/unipept/FragGeneScanRs),
a gene prediction model for short and error-prone reads.*

[![Release](https://img.shields.io/github/v/release/tomdstanton/pyfgs)](https://img.shields.io/github/v/release/tomdstanton/pyfgs)
[![Build status](https://img.shields.io/github/actions/workflow/status/tomdstanton/pyfgs/main.yml?branch=main)](https://github.com/tomdstanton/pyfgs/actions/workflows/main.yml?query=branch%3Amain)
[![Commit activity](https://img.shields.io/github/commit-activity/m/tomdstanton/pyfgs)](https://img.shields.io/github/commit-activity/m/tomdstanton/pyfgs)
[![License](https://img.shields.io/github/license/tomdstanton/pyfgs)](https://img.shields.io/github/license/tomdstanton/pyfgs)


## 🗺️ Overview

### 🔬 The Biological Edge: The Metagenomic Short-Read Specialist

- **Reads Through Sequencing Errors:** Unlike Prodigal and Pyrodigal (which are designed for pristine,
  assembled contigs), `pyfgs` uses a Hidden Markov Model trained specifically on sequencing error profiles (Illumina, 454, Sanger).

- **Native Frameshift Correction:** If a raw read contains an indel, standard tools instantly break the open reading
  frame and lose the gene. `pyfgs` detects the error, dynamically corrects the reading frame, and translates the
  protein seamlessly.

- **Granular Indel Tracking:** Every predicted gene exposes native Python lists of exactly where insertions
  and deletions were detected, allowing for rigorous downstream quality control.

### ⚡️ The Engineering Edge: Bare-Metal Rust in Python

- **GIL-Free Multithreading:** The Rust engine completely detaches the Python Global Interpreter Lock (GIL) during
  model inference. You can throw massive FASTQ files at it and watch it perfectly saturate every physical core on
  your machine.

- **True Zero-Copy Memory:** `pyfgs` doesn't waste time copying Python strings into Rust memory. The Rust backend
  borrows raw byte slices (`&[u8]`) directly from the Python interpreter's heap, resulting in a virtually
  non-existent memory footprint.

- **Lazy Byte Evaluation:** Bypasses the massive "UTF-8 Tax" of standard bioinformatics wrappers.
  Translated amino acid sequences and corrected DNA are evaluated lazily—meaning the heavy string math only happens if
  and when you explicitly request it.

- **No FFI Subprocess Tax:** Instead of dumping massive .faa files to your hard drive and parsing them back into
  Python, the HMM runs purely in memory and yields native Python objects ready for immediate downstream analysis.

### 🐍 Pythonic Quality of Life

- **0-Based BED Coordinates:** Say goodbye to wrestling with 1-based, fully closed GFF3 coordinates. `pyfgs` natively
  outputs standard 0-based, half-open intervals ([start, end)), allowing you to slice standard sequence arrays
  immediately.

- **Drop-in CLI Replacement:** Includes a hyper-fast, multithreaded command-line interface that flawlessly mimics the
  original FragGeneScan tool but operates at a fraction of the compute time.


## 🔧 Installing

This project is supported on Python 3.10 and later.

`pyfgs` can be installed directly from [PyPI](https://pypi.org/project/pyfgs/):

```console
$ pip install pyfgs
```

## 💻 CLI Usage

For API usage, please refer to the [documentation](https://tomdstanton.github.io/pyfgs/api/).
For CLI usage, type `pyfgs --help`

```console
usage: pyfgs <seq> [options]

🔗🐍⏭️	PyO3 bindings and Python interface to FragGeneScanRs,
	a gene prediction model for short and error-prone reads.

Input options 💽:

  seq             Sequence file (or '-' for stdin)
  -m, --model     Sequence error model (default: complete):
                   - short1: Illumina sequencing reads with about 0.1% error rate
                   - short5: Illumina sequencing reads with about 0.5% error rate
                   - short10: Illumina sequencing reads with about 1% error rate
                   - sanger5: Sanger sequencing reads with about 0.5% error rate
                   - sanger10: Sanger sequencing reads with about 1% error rate
                   - pyro5: 454 pyrosequencing reads with about 0.5% error rate
                   - pyro10: 454 pyrosequencing reads with about 1% error rate
                   - pyro30: 454 pyrosequencing reads with about 3% error rate
                   - complete: Complete genomic sequences or short sequence reads without sequencing error
  -r, --reads     Force FASTQ parsing (default: False)

Output options ⚙️:

  -o, --out       Output file (default: stdout)
  -f, --format    Output format (default: faa):
                   - faa (protein fasta)
                   - ffn (nucleotide fasta)
                   - bed (BED6 format)

Other options 🚧:

  -t, --threads   Number of threads (default: 8)
  -v, --version   Print version and exit
  -h, --help      Print help and exit
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
