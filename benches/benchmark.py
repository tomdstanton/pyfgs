import time
from concurrent.futures import ThreadPoolExecutor
import pyfgs


def process_record(finder: pyfgs.GeneFinder, record: tuple[str, str]) -> tuple[str, int]:
    """Worker function to process a single sequence."""
    header, sequence = record

    # This runs in pure C/Rust-land, completely freeing the Python thread!
    genes = finder.find_genes(header, sequence)

    return header, len(genes)


def main():
    fasta_file = "test_reads.fasta"  # Point this to a real file on your machine
    model = pyfgs.Model.Illumina5

    print(f"Loading model {model}...")
    finder = pyfgs.GeneFinder(model, whole_genome=False)

    print(f"Opening {fasta_file}...")
    reader = pyfgs.FastaReader(fasta_file)

    start_time = time.perf_counter()
    total_genes = 0
    total_reads = 0

    # Fire up the thread pool!
    # Python defaults to min(32, os.cpu_count() + 4) workers, which is usually perfect.
    with ThreadPoolExecutor() as executor:
        # We read records sequentially in the main thread,
        # and dispatch the heavy Viterbi math to the worker threads.
        futures = [
            executor.submit(process_record, finder, record)
            for record in reader
        ]

        # Gather the results as they finish
        for future in futures:
            header, gene_count = future.result()
            total_genes += gene_count
            total_reads += 1

    end_time = time.perf_counter()
    duration = end_time - start_time

    print("-" * 30)
    print(f"Processed {total_reads} reads.")
    print(f"Found {total_genes} ORFs.")
    print(f"Total time: {duration:.2f} seconds")
    print(f"Speed: {total_reads / duration:.2f} reads/second")


if __name__ == "__main__":
    main()