# 🔗🐍⏭️ `pyfgs` [![Stars](https://img.shields.io/github/stars/tomdstanton/pyfgs.svg?style=social&maxAge=3600&label=Star)](https://github.com/tomdstanton/pyfgs/stargazers)

*PyO3 bindings and Python interface to [FragGeneScanRs](https://github.com/unipept/FragGeneScanRs),
a gene prediction model for short and error-prone reads.*

[![Release](https://img.shields.io/github/v/release/tomdstanton/pyfgs)](https://img.shields.io/github/v/release/tomdstanton/pyfgs)
[![Build status](https://img.shields.io/github/actions/workflow/status/tomdstanton/pyfgs/main.yml?branch=main)](https://github.com/tomdstanton/pyfgs/actions/workflows/main.yml?query=branch%3Amain)
[![Commit activity](https://img.shields.io/github/commit-activity/m/tomdstanton/pyfgs)](https://img.shields.io/github/commit-activity/m/tomdstanton/pyfgs)
[![License](https://img.shields.io/github/license/tomdstanton/pyfgs)](https://img.shields.io/github/license/tomdstanton/pyfgs)


## 🗺️ Overview

`pyfgs` is a Python module that provides bindings to FragGeneScanRs using
[maturin](https://maturin.rs/). It directly interacts with the FragGeneScanRs
internals, which has the following advantages:

- **single dependency**: `pyfgs` is distributed as a Python package, so you
  can add it as a dependency to your project, and stop worrying about the
  FragGeneScanRs binary being present on the end-user machine.
- **no intermediate files**: Everything happens in memory, in a Python object
  you fully control, so you don't have to invoke the FragGeneScanRs CLI using a
  sub-process and temporary files. Sequences can be passed directly as
  strings or bytes, which avoids the overhead of formatting your input to
  FASTA for FragGeneScanRs.


## 🔧 Installing

This project is supported on Python 3.10 and later.

`pyfgs` can be installed directly from [PyPI](https://pypi.org/project/pyfgs/):

```console
$ pip install pyfgs
```

## Usage

For API usage, please refer to the documentation.
For CLI usage, type `pyfgs --help`

```console
usage: pyfgs <seq> <out> [options]

🔗🐍⏭️	PyO3 bindings and Python interface to FragGeneScanRs,
	a gene prediction model for short and error-prone reads.

Inputs and outputs 💽:

  seq            Sequence file (or '-' for stdin)
  out            Output file, defaults to stdout

Options	⚙️:

  -f, --format   Output format (default: faa):
                  - faa (protein fasta)
                  - ffn (nucleotide fasta)
                  - bed
  -m, --model    Sequence error model (default: complete):
                  - illumina_1
                  - illumina_5
                  - illumina_10
                  - sanger_5
                  - sanger_10
                  - 454_5
                  - 454_10
                  - 454_30
                  - complete

Other 🚧:

  -v, --version  Print version and exit
  -h, --help     Print help and exit
```



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

*This project is in no way affiliated, sponsored, or otherwise endorsed
by the original FragGeneScanRs authors [Peter Dawyndt](https://github.com/pdawyndt),
[Bart Mesuere](https://github.com/bmesuere) and
[Felix Van der Jeugt](https://github.com/ninewise). It was developed
by [Tom Stanton](https://github.com/tomdstanton/) during his Post-doc project
at [Monash University](https://www.monash.edu/medicine/translational/infectious-diseases) in
the [Wryes Lab](https://https://wyreslab.com/).*