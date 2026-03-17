import argparse
import sys
import os
from concurrent.futures import ThreadPoolExecutor
from pyfgs import Model, FastaReader, FastqReader, GeneFinder

# Constants ------------------------------------------------------------------------------------------------------------
_MODEL_MAP = {
    "short1": Model.Illumina1,
    "short5": Model.Illumina5,
    "short10": Model.Illumina10,
    "sanger5": Model.Sanger5,
    "sanger10": Model.Sanger10,
    "pyro5": Model.Pyro454_5,
    "pyro10": Model.Pyro454_10,
    "pyro30": Model.Pyro454_30,
    "complete": Model.Complete,
}


# Stdin streamers ------------------------------------------------------------------------------------------------------
def _FastaStreamer():
    """A lightweight, pure-Python generator for piping FASTA via stdin."""
    buffer = sys.stdin.buffer
    header = None
    seq = bytearray()

    for line in buffer:
        if line.startswith(b">"):
            if header:
                yield header, bytes(seq)
            # rstrip first to ensure \n isn't caught in the header if there are no spaces
            header, _, _ = line[1:].rstrip().partition(b' ')
            seq = bytearray()
        else:
            seq.extend(line.rstrip())

    if header:
        yield header, bytes(seq)


def _FastqStreamer():
    """A lightweight, pure-Python generator for piping FASTQ via stdin."""
    line = sys.stdin.buffer.readline
    while True:
        raw_line = line()
        if not raw_line:
            break
        header, _, _ = raw_line[1:].rstrip().partition(b' ')
        seq = line().rstrip()
        line()  # skip +
        line()  # skip qual
        yield header, seq


# Formatters simply return bytes now (No I/O) --------------------------------------------------------------------------
def _FaaFormatter(header, genes) -> bytes:
    return b''.join(
        b">%s_%d\n%s\n" % (header, n, gene.translation())
        for n, gene in enumerate(genes, start=1)
    )


def _FfnFormatter(header, genes) -> bytes:
    return b''.join(
        b">%s_%d\n%s\n" % (header, n, gene.sequence())
        for n, gene in enumerate(genes, start=1)
    )


def _BedFormatter(header, genes) -> bytes:
    return b''.join(
        b"%s\t%d\t%d\tID=%s_%d\t%.6f\t%s\n" % (
            header, gene.start, gene.end, header, n, gene.score, b'+' if gene.strand == 1 else b'-'
        ) for n, gene in enumerate(genes, start=1)
    )


# Helper functions -----------------------------------------------------------------------------------------------------
def _get_optimal_threads():
    """Safely determine the number of available cores, respecting HPC/Docker limits."""
    # Modern Python (3.13+) handles cgroups and Docker perfectly
    if hasattr(os, "process_cpu_count"):
        return os.process_cpu_count() or 1

    # Fallback for older Pythons on Linux/HPC (respects SLURM/taskset)
    if hasattr(os, "sched_getaffinity"):
        return len(os.sched_getaffinity(0))

    # Absolute fallback to motherboard cores
    return os.cpu_count() or 1


# Main entry point -----------------------------------------------------------------------------------------------------
def main():

    # Define args ------------------------------------------------------------------------------------------------------
    parser = argparse.ArgumentParser(
        prog="pyfgs", usage="%(prog)s <seq> [options]", add_help=False,
        description="🔗🐍⏭️\tPyO3 bindings and Python interface to FragGeneScanRs,\n"
                    "\ta gene prediction model for short and error-prone reads.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    in_group = parser.add_argument_group("Input options 💽", "")
    in_group.add_argument("seq", help="Sequence file (or '-' for stdin)")
    in_group.add_argument("-m", "--model", metavar='', choices=_MODEL_MAP.keys(), default="complete",
                          help="Sequence error model (default: complete):\n"
                               " - short1: Illumina sequencing reads with about 0.1%% error rate\n"
                               " - short5: Illumina sequencing reads with about 0.5%% error rate\n"
                               " - short10: Illumina sequencing reads with about 1%% error rate\n"
                               " - sanger5: Sanger sequencing reads with about 0.5%% error rate\n"
                               " - sanger10: Sanger sequencing reads with about 1%% error rate\n"
                               " - pyro5: 454 pyrosequencing reads with about 0.5%% error rate\n"
                               " - pyro10: 454 pyrosequencing reads with about 1%% error rate\n"
                               " - pyro30: 454 pyrosequencing reads with about 3%% error rate\n"
                               " - complete: Complete genomic sequences or short sequence reads without sequencing error\n"
                          )
    in_group.add_argument("-r", "--reads", action="store_true", help="Force FASTQ parsing (default: False)")

    out_group = parser.add_argument_group("Output options ⚙️", "")
    out_group.add_argument("-o", "--out", help="Output file (default: stdout)", metavar='', default="-")
    out_group.add_argument("-f", "--format", metavar='', choices=['faa', 'ffn', 'bed'], default="faa",
                            help="Output format (default: faa):\n"
                                 " - faa (protein fasta)\n"
                                 " - ffn (nucleotide fasta)\n"
                                 " - bed (BED6 format)\n")

    other_group = parser.add_argument_group("Other options 🚧", "")
    other_group.add_argument("-t", "--threads", type=int, default=_get_optimal_threads(), metavar="",
                            help="Number of threads (default: %(default)s)")
    other_group.add_argument("-v", "--version", action='version', version="pyfgs 0.1.0", help="Print version and exit")
    other_group.add_argument("-h", "--help", action='help', help="Print help and exit")


    # Parse args -------------------------------------------------------------------------------------------------------
    args = parser.parse_args()

    # Init GeneFinder --------------------------------------------------------------------------------------------------
    finder = GeneFinder(_MODEL_MAP[args.model])

    # Init input -------------------------------------------------------------------------------------------------------
    if args.seq == "-":
        reader = _FastqStreamer() if args.reads else _FastaStreamer()
    else:
        reader = (FastqReader if args.reads else FastaReader)(args.seq)

    # Init output ------------------------------------------------------------------------------------------------------
    formatter = {"faa": _FaaFormatter, "ffn": _FfnFormatter, "bed": _BedFormatter}[args.format]
    out_f = sys.stdout.buffer if args.out == "-" else open(args.out, "wb")

    # Begin ------------------------------------------------------------------------------------------------------------
    def _process_record(record):
        """Worker function for the thread pool: pure computation."""
        header, seq = record[0], record[1]
        genes = finder.find_genes(header, seq)
        # We only spend time formatting if genes were actually found
        return formatter(header, genes) if genes else b""

    try:
        with ThreadPoolExecutor(max_workers=args.threads) as executor:
            # executor.map loops in sequence order. The main thread handles the I/O!
            for output_bytes in executor.map(_process_record, reader):
                if output_bytes:
                    out_f.write(output_bytes)
    finally:
        if args.out != "-":
            out_f.close()


if __name__ == "__main__":
    main()