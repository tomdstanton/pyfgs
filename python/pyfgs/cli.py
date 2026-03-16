import argparse
import sys
from concurrent.futures import ThreadPoolExecutor
from pyfgs import Model, FastaReader, FastqReader, GeneFinder

# Constants ------------------------------------------------------------------------------------------------------------
_MODEL_MAP = {
    "illumina_1": Model.Illumina1,
    "illumina_5": Model.Illumina5,
    "illumina_10": Model.Illumina10,
    "sanger_5": Model.Sanger5,
    "sanger_10": Model.Sanger10,
    "454_5": Model.Pyro454_5,
    "454_10": Model.Pyro454_10,
    "454_30": Model.Pyro454_30,
    "complete": Model.Complete,
}


# Functions ------------------------------------------------------------------------------------------------------------
def _stream_fasta():
    """A lightweight, pure-Python generator for piping FASTA via stdin."""
    header = b""
    seq = bytearray()
    
    # Read the raw byte stream to avoid UTF-8 overhead
    for line in sys.stdin.buffer:
        if line.startswith(b">"):
            if header:
                yield header, bytes(seq)
            # Match seq_io behavior: grab ID up to the first space
            header = line[1:].split(b" ", 1)[0].strip()
            seq = bytearray()
        else:
            seq.extend(line.strip())
            
    if header:
        yield header, bytes(seq)


def _stream_fastq():
    pass


def _process_record(finder, record: tuple[bytes, bytes]):
    """Worker function for the thread pool."""
    # record is either (header, seq) or (header, seq, qual)
    header = record[0]
    seq = record[1]
    genes = finder.find_genes(header, seq)
    return header, genes


def main():

    parser = argparse.ArgumentParser(
        prog="pyfgs", usage="%(prog)s <seq> <out> [options]", add_help=False,
        description="🔗🐍⏭️\tPyO3 bindings and Python interface to FragGeneScanRs,\n"
                    "\ta gene prediction model for short and error-prone reads.",
        formatter_class = argparse.RawTextHelpFormatter
    )
    io_group = parser.add_argument_group("Inputs and outputs 💽", "")
    io_group.add_argument("seq", help="Sequence file (or '-' for stdin)")
    io_group.add_argument("out", help="Output file, defaults to stdout", default="-")
    opts_group = parser.add_argument_group("Options\t⚙️", "")
    opts_group.add_argument("-f", "--format", metavar='', choices=['faa', 'ffn', 'bed'],
                            help="Output format (default: faa):\n"
                                 " - faa (protein fasta)\n"
                                 " - ffn (nucleotide fasta)\n"
                                 " - bed",
                            default="faa")
    opts_group.add_argument("-m", "--model", metavar='', choices=_MODEL_MAP.keys(), default="complete",
                        help="Sequence error model (default: complete):" + "".join('\n - ' + i for i in _MODEL_MAP.keys()))
    other_group = parser.add_argument_group("Other 🚧", "")
    other_group.add_argument("-v", "--version", action='version', help="Print version and exit")
    other_group.add_argument("-h", "--help", action='help', help="Print help and exit")

    
    args = parser.parse_args()

    model = _MODEL_MAP[args.model]
    finder = GeneFinder(model)
    if args.seq == "-":
        reader = (_stream_fasta if model == Model.Complete else _stream_fastq)
    else:
        reader = (FastaReader if model == Model.Complete else FastqReader)(args.seq)

    with ThreadPoolExecutor(max_workers=args.thread) as executor:
        # We use executor.map to guarantee the outputs are written in the
        # exact same order as the input file, just like the original CLI!

        # A tiny wrapper lambda to pass the finder along with the record
        futures = executor.map(lambda rec: _process_record(finder, rec), reader)
        pass
