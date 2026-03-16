# ⏭️ pyfgs [![Stars](https://img.shields.io/github/stars/tomdstanton/pyfgs.svg?style=social&maxAge=3600&label=Star)](https://github.com/tomdstanton/pyfgs/stargazers)

*PyO3 bindings and Python interface to [FragGeneScanRs](https://github.com/unipept/FragGeneScanRs),
a gene prediction model for short and error-prone reads.*

[![Release](https://img.shields.io/github/v/release/tomdstanton/pyfgs)](https://img.shields.io/github/v/release/tomdstanton/pyfgs)
[![Build status](https://img.shields.io/github/actions/workflow/status/tomdstanton/pyfgs/main.yml?branch=main)](https://github.com/tomdstanton/pyfgs/actions/workflows/main.yml?query=branch%3Amain)
[![Commit activity](https://img.shields.io/github/commit-activity/m/tomdstanton/pyfgs)](https://img.shields.io/github/commit-activity/m/tomdstanton/pyfgs)
[![License](https://img.shields.io/github/license/tomdstanton/pyfgs)](https://img.shields.io/github/license/tomdstanton/pyfgs)


## 🗺️ Overview

pyfgs is a Python module that provides bindings to FragGeneScanRs using
[Cython](https://cython.org/). It directly interacts with the FragGeneScanRs
internals, which has the following advantages:

- **single dependency**: pyfgs is distributed as a Python package, so you
  can add it as a dependency to your project, and stop worrying about the
  FragGeneScanRs binary being present on the end-user machine.
- **no intermediate files**: Everything happens in memory, in a Python object
  you fully control, so you don't have to invoke the FragGeneScanRs CLI using a
  sub-process and temporary files. Sequences can be passed directly as
  strings or bytes, which avoids the overhead of formatting your input to
  FASTA for FragGeneScanRs.

### 📋 Features

The library features everything from the original FragGeneScanRs CLI:

- **run mode selection**: Choose between *single* mode, using a training
  sequence to count nucleotide hexamers, or *metagenomic* mode, using
  pre-trained data from different organisms (`FragGeneScanRs -p`).
- **region masking**: Prevent genes from being predicted across regions
  containing unknown nucleotides  (`FragGeneScanRs -m`).
- **closed ends**: Genes will be identified as running over edges if they
  are larger than a certain size, but this can be disabled (`FragGeneScanRs -c`).
- **training configuration**: During the training process, a custom
  translation table can be given (`FragGeneScanRs -g`), and the Shine-Dalgarno motif
  search can be forcefully bypassed (`FragGeneScanRs -n`)
- **output files**: Output files can be written in a format mostly
  compatible with the FragGeneScanRs binary, including the protein translations
  in FASTA format (`FragGeneScanRs -a`), the gene sequences in FASTA format
  (`FragGeneScanRs -d`), or the potential gene scores in tabular format
  (`FragGeneScanRs -s`).
- **training data persistence**: Getting training data from a sequence and
  using it for other sequences is supported; in addition, a training data
  file can be saved and loaded transparently (`FragGeneScanRs -t`).

In addition, the **new** features are available:

- **custom gene size threshold**: While FragGeneScanRs uses a minimum gene size
  of 90 nucleotides (60 if on edge), pyfgs allows to customize this
  threshold, allowing for smaller ORFs to be identified if needed.
- **custom metagenomic models**: Since `v3.0.0`, you can use your own 
  metagenomic models to run pyfgs in *meta*-mode. *Check for instance
  [`pyfgs-gv`](https://github.com/tomdstanton/pyfgs-gv) by 
  [Antônio Camargo](https://github.com/apcamargo/), which 
  provides additional models for giant viruses and gut phages,
  or [`pyfgs-rv`](https://github.com/LanderDC/pyfgs-rv)
  by [Lander De Coninck](https://github.com/LanderDC) which
  provides additional models for RNA viruses.*

### 🐏 Memory

pyfgs makes several changes compared to the original FragGeneScanRs binary
regarding memory management:

* Sequences are stored as raw bytes instead of compressed bitmaps. This means
  that the sequence itself takes 3/8th more space, but since the memory used
  for storing the sequence is often negligible compared to the memory used to
  store dynamic programming nodes, this is an acceptable trade-off for better
  performance when extracting said nodes.
* Node fields use smaller data types to fit into 128 bytes, compared to the 
  176 bytes of the original FragGeneScanRs data structure.
* Node arrays are pre-allocated based on the sequence GC% to extrapolate the
  probability to find a start or stop codon.
* Genes are stored in a more compact data structure than in FragGeneScanRs (which
  reserves a buffer to store string data), saving around 1KiB per gene.


### 🧶 Thread-safety

[`pyfgs.GeneFinder`](https://pyfgs.readthedocs.io/en/stable/api/gene_finder.html#pyfgs.GeneFinder)
instances are thread-safe. In addition, the
[`find_genes`](https://pyfgs.readthedocs.io/en/stable/api/gene_finder.html#pyfgs.GeneFinder.find_genes)
method is re-entrant. This means you can train an
[`GeneFinder`](https://pyfgs.readthedocs.io/en/stable/api/gene_finder.html#pyfgs.GeneFinder)
instance once, and then use a pool to process sequences in parallel:
```python
import multiprocessing.pool
import pyfgs

gene_finder = pyfgs.GeneFinder()
gene_finder.train(training_sequence)

with multiprocessing.pool.ThreadPool() as pool:
    predictions = pool.map(gene_finder.find_genes, sequences)
```

## 🔧 Installing

This project is supported on Python 3.7 and later.

pyfgs can be installed directly from [PyPI](https://pypi.org/project/pyfgs/),
which hosts some pre-built wheels for the x86-64 architecture (Linux/MacOS/Windows)
and the Aarch64 architecture (Linux/MacOS), as well as the code required to compile
from source with Cython:
```console
$ pip install pyfgs
```

Otherwise, pyfgs is also available as a [Bioconda](https://bioconda.github.io/)
package:
```console
$ conda install -c bioconda pyfgs
```

Check the [*install* page](https://pyfgs.readthedocs.io/en/stable/guide/install.html)
of the documentation for other ways to install pyfgs on your machine.

## 💡 Example

Let's load a sequence from a
[GenBank](http://www.insdc.org/files/feature_table.html) file, use a `GeneFinder`
to find all the genes it contains, and print the proteins in two-line FASTA
format.

### 🔬 [Biopython](https://github.com/biopython/biopython)

To use the [`GeneFinder`](https://pyfgs.readthedocs.io/en/stable/api/gene_finder.html#pyfgs.GeneFinder)
in single mode (corresponding to `FragGeneScanRs -p single`, the default operation mode of FragGeneScanRs),
you must explicitly call the
[`train`](https://pyfgs.readthedocs.io/en/stable/api/gene_finder.html#pyfgs.GeneFinder.train) method
with the sequence you want to use for training before trying to find genes,
or you will get a [`RuntimeError`](https://docs.python.org/3/library/exceptions.html#RuntimeError):
```python
import Bio.SeqIO
import pyfgs

record = Bio.SeqIO.read("sequence.gbk", "genbank")

gene_finder = pyfgs.GeneFinder()
gene_finder.train(bytes(record.seq))
genes = gene_finder.find_genes(bytes(record.seq))
```

However, in `meta` mode (corresponding to `FragGeneScanRs -p meta`), you can find genes directly:
```python
import Bio.SeqIO
import pyfgs

record = Bio.SeqIO.read("sequence.gbk", "genbank")

gene_finder = pyfgs.GeneFinder(meta=True)
for i, pred in enumerate(gene_finder.find_genes(bytes(record.seq))):
    print(f">{record.id}_{i+1}")
    print(pred.translate())
```

*On older versions of Biopython (before 1.79) you will need to use
`record.seq.encode()` instead of `bytes(record.seq)`*.


## 🔖 Citation

For now, please cite the original [FragGeneScanRs paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04736-5):

> Van der Jeugt, F., Dawyndt, P. & Mesuere, B. FragGeneScanRs: faster gene prediction for short reads. BMC Bioinformatics 23, 198 (2022). https://doi.org/10.1186/s12859-022-04736-5


## 💭 Feedback

### ⚠️ Issue Tracker

Found a bug ? Have an enhancement request ? Head over to the [GitHub issue
tracker](https://github.com/tomdstanton/pyfgs/issues) if you need to report
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

*This project is in no way not affiliated, sponsored, or otherwise endorsed
by the original FragGeneScanRs authors [Peter Dawyndt](https://github.com/pdawyndt),
[Bart Mesuere](https://github.com/bmesuere) and
[Felix Van der Jeugt](https://github.com/ninewise). It was developed
by [Tom Stanton](https://github.com/tomdstanton/) during his Post-doc project
at [Monash University](https://www.monash.edu/medicine/translational/infectious-diseases) in
the [Wryes Lab](https://https://wyreslab.com/).*