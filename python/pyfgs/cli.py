"""Command line interface for `pyfgs` and main entry point"""
import typing
from importlib import resources, metadata
import argparse
from sys import stdin, stdout
from pathlib import Path
import logging

"""
USAGE:
    pyfgs [FLAGS] [OPTIONS] --training-file <train_file_name>

FLAGS:
    -f, --formatted    Format the DNA output.
    -h, --help         Prints help information
    -u, --unordered    Do not preserve record order in output (faster).
    -V, --version      Prints version information

OPTIONS:
    -a, --aa-file <aa_file>                    Output predicted proteins to this file (supersedes -o).
    -w, --complete <complete>                  The input sequence has complete genomic sequences; not short sequence
                                               reads. [default: 0]
    -g, --gff-file <gff_file>                  Output metadata to this gff formatted file (supersedes -o).
    -m, --meta-file <meta_file>                Output metadata to this file (supersedes -o).
    -n, --nucleotide-file <nucleotide_file>    Output predicted genes to this file (supersedes -o).
    -o, --output-prefix <output_prefix>        Output metadata (.out and .gff), proteins (.faa) and genes (.ffn) to
                                               files with this prefix. Use 'stdout' to write the predicted proteins to
                                               standard output.
    -s, --seq-file-name <seq_file_name>        Sequence file name including the full path. Using 'stdin' (or not
                                               suplying this argument) reads from standard input. [default: stdin]
    -p, --thread-num <thread_num>              The number of threads used by FragGeneScan++. [default: 1]
    -t, --training-file <train_file_name>      File name that contains model parameters; this file should be in the -r
                                               directory or one of the following:
                                               [complete] for complete genomic sequences or short sequence reads without
                                               sequencing error
                                               [sanger_5] for Sanger sequencing reads with about 0.5% error rate
                                               [sanger_10] for Sanger sequencing reads with about 1% error rate
                                               [454_5] for 454 pyrosequencing reads with about 0.5% error rate
                                               [454_10] for 454 pyrosequencing reads with about 1% error rate
                                               [454_30] for 454 pyrosequencing reads with about 3% error rate
                                               [illumina_1] for Illumina sequencing reads with about 0.1% error rate
                                               [illumina_5] for Illumina sequencing reads with about 0.5% error rate
                                               [illumina_10] for Illumina sequencing reads with about 1% error rate
"""


def cli(): pass



