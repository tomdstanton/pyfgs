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


# Formatters -----------------------------------------------------------------------------------------------------------
def _FaaFormatter(header, seq, genes) -> bytes:
    return b''.join(
        b">%s_FGS_%d\n%s\n" % (header, n, gene.translation())
        for n, gene in enumerate(genes, start=1)
    )


def _FnaFormatter(header, seq, genes) -> bytes:
    return b''.join(
        b">%s_FGS_%d\n%s\n" % (header, n, gene.sequence())
        for n, gene in enumerate(genes, start=1)
    )


def _BedFormatter(header, seq, genes) -> bytes:
    lines = []
    for n, gene in enumerate(genes, start=1):
        muts = gene.mutations(seq)
        mut_str = b"."
        if muts:
            mut_info = [f"pos={m.pos};type={m.mut_type};codon={m.codon_idx};note={m.annotation}".encode() for m in muts]
            mut_str = b",".join(mut_info)

        strand = b"+" if gene.strand == 1 else b"-"
        lines.append(b"%s\t%d\t%d\t%s_FGS_%d\t%.2f\t%s\t%s\n" % (
            header, gene.start, gene.end, header, n, gene.score, strand, mut_str
        ))
    return b"".join(lines)


def _GffFormatter(header, seq, genes) -> bytes:
    lines = []
    for n, gene in enumerate(genes, start=1):
        muts = gene.mutations(seq)
        start = gene.start + 1  # GFF3 is 1-based closed
        strand = b"+" if gene.strand == 1 else b"-"
        gene_id = b"%s_FGS_%d" % (header, n)

        if muts:
            notes = b",".join(b"Frameshift %b at pos %d (codon %d)" % (
                b"insertion" if m.mut_type == "ins" else b"deletion", m.pos, m.codon_idx
            ) for m in muts)

            lines.append(
                b"%s\tpyfgs\tpseudogene\t%d\t%d\t%.2f\t%s\t0\tID=%s;inference=ab initio prediction:pyfgs;pseudogene=unknown;Note=%s\n" % (
                    header, start, gene.end, gene.score, strand, gene_id, notes
                ))
        else:
            lines.append(b"%s\tpyfgs\tCDS\t%d\t%d\t%.2f\t%s\t0\tID=%s;inference=ab initio prediction:pyfgs\n" % (
                header, start, gene.end, gene.score, strand, gene_id
            ))
    return b"".join(lines)


def _VcfFormatter(header, seq, genes) -> bytes:
    lines = []
    for n, gene in enumerate(genes, start=1):
        muts = gene.mutations(seq)
        if not muts:
            continue
            
        gene_id = b"%s_FGS_%d" % (header, n)
        
        for m in muts:
            # Convert strings to bytes for fast I/O
            mut_type = b"insertion" if m.mut_type == "ins" else b"deletion"
            ref = m.ref_allele.encode()
            alt = m.alt_allele.encode()
            ann = m.annotation.encode()
            
            # Build the SnpEff compliant ANN string
            info = b"TYPE=frameshift_%s;GENE=%s;CODON=%d;ANN=%s|frameshift_variant|HIGH|%s|%s|transcript|%s|protein_coding|%d/1|%s" % (
                mut_type, gene_id, m.codon_idx, alt, gene_id, gene_id, gene_id, m.codon_idx, ann
            )
            
            lines.append(b"%s\t%d\t.\t%s\t%s\t.\tPASS\t%s\n" % (
                header, m.pos, ref, alt, info
            ))
            
    return b"".join(lines)


# Helper functions -----------------------------------------------------------------------------------------------------
def _get_optimal_threads() -> int:
    if hasattr(os, "process_cpu_count"):
        return os.process_cpu_count() or 1
    if hasattr(os, "sched_getaffinity"):
        return len(os.sched_getaffinity(0))
    return os.cpu_count() or 1


def _is_fastq(seq_arg, force_reads_flag) -> bool:
    """Robust FASTQ detection via explicit flag, file extension, or byte sniffing."""
    if force_reads_flag:
        return True

    if seq_arg == "-":
        return False  # Safely default stdin to FASTA unless forced

    ext = seq_arg.lower().rstrip('.gz')
    if ext.endswith(('.fq', '.fastq')):
        return True
    if ext.endswith(('.fa', '.fna', '.fasta', '.faa')):
        return False

    # If extension is ambiguous, sniff the first byte!
    try:
        with open(seq_arg, 'rb') as f:
            return f.read(1) == b'@'
    except OSError:
        return False


# Main entry point -----------------------------------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        prog="pyfgs", usage="%(prog)s <seq> [options]", add_help=False,
        description="🔗🐍⏭️\tPyO3 bindings and Python interface to FragGeneScanRs,\n"
                    "\ta gene prediction model for short and error-prone reads.",
        formatter_class=argparse.RawTextHelpFormatter
    )

    in_group = parser.add_argument_group("Input options 💽", "")
    in_group.add_argument("seq", help="Sequence file (or '-' for stdin)")
    in_group.add_argument("-m", "--model", metavar='', choices=_MODEL_MAP.keys(), default="complete",
                          help="Sequence error model (default: complete)\n"
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
    # Updated help text to show it overrides the auto-detection
    in_group.add_argument("-r", "--reads", action="store_true",
                          help="Force FASTQ parsing (Overrides auto-detection)")
    in_group.add_argument("-w", "--whole-genome", action="store_true",
                          help="Strict contiguous ORFs. Disables error-tolerant frameshift detection.")

    out_group = parser.add_argument_group("Output options ⚙️",
                                          "Provide a PATH to save to a file, or use the flag alone to print to stdout.")
    out_group.add_argument("--faa", nargs='?', const='-', default=None, metavar="PATH", help="Output protein FASTA")
    out_group.add_argument("--fna", nargs='?', const='-', default=None, metavar="PATH", help="Output nucleotide FASTA")
    out_group.add_argument("--bed", nargs='?', const='-', default=None, metavar="PATH", help="Output BED6+1 format")
    out_group.add_argument("--gff", nargs='?', const='-', default=None, metavar="PATH", help="Output GFF3 format")
    out_group.add_argument("--vcf", nargs='?', const='-', default=None, metavar="PATH", help="Output VCF v4.2 format")

    other_group = parser.add_argument_group("Other options 🚧", "")
    other_group.add_argument("-t", "--threads", type=int, default=_get_optimal_threads(), metavar="",
                             help="Number of threads (default: optimal)")
    other_group.add_argument("-v", "--version", action='version', version="pyfgs 0.1.0", help="Print version and exit")
    other_group.add_argument("-h", "--help", action='help', help="Print help and exit")

    args = parser.parse_args()

    # 1. Output Resolution & Validation --------------------------------------------------------------------------------
    formats_map = {"faa": _FaaFormatter, "fna": _FnaFormatter, "bed": _BedFormatter, "gff": _GffFormatter, "vcf": _VcfFormatter}

    active_outputs = {fmt: getattr(args, fmt) for fmt in formats_map if getattr(args, fmt) is not None}

    if not active_outputs:
        active_outputs["faa"] = "-"

    stdouts = [fmt for fmt, dest in active_outputs.items() if dest == "-"]
    if len(stdouts) > 1:
        sys.stderr.write(
            f"⚠️ WARNING: Multiple outputs ({', '.join(stdouts)}) are writing to stdout. Streams will interleave!\n")

    out_files = {}
    for fmt, path in active_outputs.items():
        handle = sys.stdout.buffer if path == "-" else open(path, "wb")
        out_files[fmt] = handle
        
        # Inject standard headers for tabular formats
        if fmt == "bed":
            handle.write(b"# source=ab initio prediction: pyfgs\n# CHROM\tSTART\tEND\tNAME\tSCORE\tSTRAND\tMUTATIONS\n")
        elif fmt == "gff":
            handle.write(b"##gff-version 3\n#source-ontology=ab initio prediction: pyfgs\n")
        elif fmt == "vcf":
            handle.write(b"##fileformat=VCFv4.2\n##source=pyfgs_ab_initio\n"
                         b"##INFO=<ID=TYPE,Number=1,Type=String,Description=\"Type of sequence discrepancy\">\n"
                         b"##INFO=<ID=GENE,Number=1,Type=String,Description=\"Predicted Gene ID\">\n"
                         b"##INFO=<ID=CODON,Number=1,Type=Integer,Description=\"1-based codon index of the frameshift\">\n"
                         b"##INFO=<ID=ANN,Number=.,Type=String,Description=\"Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_Biotype | Rank | HGVS.c'\">\n"
                         b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    # 2. Init Engine & Input -------------------------------------------------------------------------------------------
    finder = GeneFinder(_MODEL_MAP[args.model], whole_genome=args.whole_genome)

    # Automatically resolve FASTQ parsing!
    is_fastq = _is_fastq(args.seq, args.reads)

    if args.seq == "-":
        reader = _FastqStreamer() if is_fastq else _FastaStreamer()
    else:
        reader = (FastqReader if is_fastq else FastaReader)(args.seq)

    # 3. Execution Pipeline --------------------------------------------------------------------------------------------
    def _process_record(record):
        header, seq = record[0], record[1]
        genes = finder.find_genes(seq) 
        if not genes:
            return None
        return {fmt: formats_map[fmt](header, seq, genes) for fmt in active_outputs}
    

    try:
        with ThreadPoolExecutor(max_workers=args.threads) as executor:
            for formatted_results in executor.map(_process_record, reader):
                if formatted_results:
                    for fmt, byte_string in formatted_results.items():
                        out_files[fmt].write(byte_string)
    finally:
        for path, handle in zip(active_outputs.values(), out_files.values()):
            if path != "-":
                handle.close()


if __name__ == "__main__":
    main()
