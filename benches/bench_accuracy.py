import pytest
import urllib.request
import gzip
from pathlib import Path
import pyfgs
import pyrodigal

_BASE_URL = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000'

def download_and_extract(url: str, name: str, ext: str) -> Path:
    """Downloads a gzipped file from NCBI and extracts it safely by name."""
    file_path = Path(f"{name}_benchmark.{ext}")
    gz_path = file_path.with_suffix(f".{ext}.gz")

    if not file_path.exists():
        urllib.request.urlretrieve(url, gz_path)
        with gzip.open(gz_path, 'rb') as f_in, file_path.open('wb') as f_out:
            f_out.write(f_in.read())
        gz_path.unlink()

    return file_path


def fetch_data(url: str, name: str):
    fna_path = download_and_extract(f'{url}.fna.gz', name, "fna")
    gff_path = download_and_extract(f'{url}.gff.gz', name, "gff")

    # Parse True Stops
    true_stops = set()
    with gff_path.open('r') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.split('\t')
            if len(parts) < 9 or parts[2] != 'CDS': continue

            seq_id = parts[0]
            start = int(parts[3]) - 1
            end = int(parts[4])
            strand = 1 if parts[6] == '+' else -1

            stop = end if strand == 1 else start
            true_stops.add((seq_id, strand, stop))

    # THE FIX: Parse Contigs AND scrub newlines/whitespace so pyrodigal doesn't choke!
    contigs = []
    for header, seq_bytes in pyfgs.FastaReader(str(fna_path)):
        seq_id = header.decode('ascii').split()[0]

        # Strip all whitespace and force uppercase
        clean_seq_str = seq_bytes.decode('ascii').replace('\n', '').replace('\r', '').replace(' ', '').upper()

        # Re-encode to bytes for the benchmark loop
        contigs.append((seq_id, clean_seq_str.encode('ascii')))

    return contigs, true_stops

@pytest.fixture(scope="session")
def genomes() -> dict[str, tuple]:
    """Session fixture to download data and parse GFF only ONCE per test suite run."""
    print("\nFetching Genomes from NCBI...")
    return {
        'E_coli': fetch_data(f"{_BASE_URL}/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic", "E_coli"),
        'S_aureus': fetch_data(f"{_BASE_URL}/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic", "S_aureus"),
        'P_aeruginosa': fetch_data(f"{_BASE_URL}/006/765/GCF_000006765.1_ASM676v1/GCF_000006765.1_ASM676v1_genomic", "P_aeruginosa")
    }


# The parameterize decorator feeds the genome names into the test function one by one
@pytest.mark.parametrize("genome_name", ['E_coli', 'S_aureus', 'P_aeruginosa'])
def test_pyrodigal_accuracy(benchmark, genomes, genome_name):
    """Benchmarks pyrodigal speed and calculates accuracy."""
    contigs, true_stops = genomes[genome_name]
    finder = pyrodigal.GeneFinder(meta=False)

    def run_inference():
        stops = set()
        for seq_id, seq_bytes in contigs:
            # Pyrodigal requires a python string
            seq_str = seq_bytes.decode('ascii')

            # THE FIX: We must train the dynamic programming model on the sequence first!
            finder.train(seq_str)
            genes = finder.find_genes(seq_str)

            for g in genes:
                stop = g.end if g.strand == 1 else (g.begin - 1)
                stops.add((seq_id, g.strand, stop))
        return stops

    pyro_stops = benchmark(run_inference)

    # Calculate Accuracy
    tp = len(pyro_stops & true_stops)
    fn = len(true_stops - pyro_stops)
    fp = len(pyro_stops - true_stops)

    print(f"\n--- pyrodigal Accuracy for {genome_name} ---")
    print(f"Ground Truth: {len(true_stops)} | Perfect Matches: {tp} | Missed: {fn} | Extra: {fp}")


@pytest.mark.parametrize("genome_name", ['E_coli', 'S_aureus', 'P_aeruginosa'])
def test_pyfgs_accuracy(benchmark, genomes, genome_name):
    """Benchmarks pyfgs speed and calculates accuracy."""
    contigs, true_stops = genomes[genome_name]
    finder = pyfgs.GeneFinder(pyfgs.Model.Complete, whole_genome=False)

    def run_inference():
        stops = set()
        for seq_id, seq_bytes in contigs:
            genes = finder.find_genes(seq_bytes)
            for g in genes:
                stop = g.end if g.strand == 1 else g.start
                stops.add((seq_id, g.strand, stop))
        return stops

    pyfgs_stops = benchmark(run_inference)

    # Calculate Accuracy
    tp = len(pyfgs_stops & true_stops)
    fn = len(true_stops - pyfgs_stops)
    fp = len(pyfgs_stops - true_stops)

    print(f"\n--- pyfgs Accuracy for {genome_name} ---")
    print(f"Ground Truth: {len(true_stops)} | Perfect Matches: {tp} | Missed: {fn} | Extra: {fp}")