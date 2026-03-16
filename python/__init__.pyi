from typing import List, Iterator, Tuple


class Model:
    """The available sequencing error models for FragGeneScanRs."""
    Illumina1: 'Model'
    Illumina5: 'Model'
    Illumina10: 'Model'
    Sanger5: 'Model'
    Sanger10: 'Model'
    Pyro454_5: 'Model'
    Pyro454_10: 'Model'
    Pyro454_30: 'Model'
    Complete: 'Model'


class FastaReader:
    """A memory-efficient FASTA parser yielding (header, sequence)."""

    def __init__(self, path: str) -> None:
        """Open a FASTA file for reading."""
        ...

    def __iter__(self) -> Iterator[Tuple[str, str]]:
        ...

    def __next__(self) -> Tuple[str, str]: ...



class FastqReader:
    """A memory-efficient FASTQ parser yielding (header, sequence, qualities)."""

    def __init__(self, path: str) -> None:
        """Open a FASTQ file for reading."""
        ...

    def __iter__(self) -> Iterator[Tuple[str, str, str]]:
        ...

    def __next__(self) -> Tuple[str, str, str]: ...


class Gene:
    """Represents a single predicted Open Reading Frame (ORF)."""

    @property
    def start(self) -> int: ...

    @property
    def end(self) -> int: ...

    @property
    def strand(self) -> str: ...

    @property
    def frame(self) -> int: ...

    @property
    def score(self) -> float: ...


class GeneFinder:
    """The main engine for finding genes, holding the HMM in memory."""

    def __init__(self, model: Model, whole_genome: bool = False) -> None:
        """
        Initialize the GeneFinder.

        Args:
            model: The sequencing error model to use (e.g., pyfgs.Model.Illumina5).
            whole_genome: Set to True if analyzing complete genomic sequences
                          rather than short reads.
        """
        ...

    def find_genes(self, header: str, sequence: str) -> List[Gene]:
        """
        Predict open reading frames in a given DNA sequence.

        This method releases the GIL, allowing for safe multi-threading
        across multiple CPU cores.

        Args:
            header: The sequence identifier.
            sequence: The raw nucleotide string.

        Returns:
            A list of predicted Gene objects.
        """
        ...