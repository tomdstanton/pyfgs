from typing import List, Iterator, Tuple, Optional, Literal


class Model:
    """The available sequencing error models for FragGeneScanRs.
    
    Attributes:
        Illumina1: Illumina sequencing reads with about 0.1% error rate.
        Illumina5: Illumina sequencing reads with about 0.5% error rate.
        Illumina10: Illumina sequencing reads with about 1% error rate.
        Sanger5: Sanger sequencing reads with about 0.5% error rate.
        Sanger10: Sanger sequencing reads with about 1% error rate.
        Pyro454_5: 454 pyrosequencing reads with about 0.5% error rate.
        Pyro454_10: 454 pyrosequencing reads with about 1% error rate.
        Pyro454_30: 454 pyrosequencing reads with about 3% error rate.
        Complete: Complete genomic sequences or short sequence reads without sequencing error.
        
    Example:
        ```python
        from pyfgs import Model
        model = Model.Complete
        ```
    """
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
    """A memory-efficient FASTA parser yielding (header, sequence).
    
    Example:
        ```python
        from pyfgs import FastaReader
        
        reader = FastaReader("genome.fasta")
        for header, sequence in reader:
            print(f">{header}\\n{sequence}")
        ```
    """

    def __init__(self, path: str) -> None:
        """Open a FASTA file for reading.
        
        Args:
            path (str): The path to the FASTA file.
        """
        ...

    def __iter__(self) -> Iterator[Tuple[bytes, bytes]]:
        ...

    def __next__(self) -> Tuple[bytes, bytes]: ...



class FastqReader:
    """A memory-efficient FASTQ parser yielding (header, sequence, qualities).
    
    Example:
        ```python
        from pyfgs import FastqReader
        
        reader = FastqReader("reads.fastq")
        for header, sequence, qualities in reader:
            print(f"@{header}\\n{sequence}\\n+\\n{qualities}")
        ```
    """

    def __init__(self, path: str) -> None:
        """Open a FASTQ file for reading.
        
        Args:
            path (str): The path to the FASTQ file.
        """
        ...

    def __iter__(self) -> Iterator[Tuple[bytes, bytes, bytes]]:
        ...

    def __next__(self) -> Tuple[bytes, bytes, bytes]: ...


class Gene:
    """Represents a single predicted Open Reading Frame (ORF).
    """

    @property
    def start(self) -> int:
        """The 0-based start coordinate of the ORF."""
        ...

    @property
    def end(self) -> int:
        """The 0-based, half-open end coordinate of the ORF."""
        ...

    @property
    def strand(self) -> int:
        """The strand of the ORF (1 for forward, -1 for reverse)."""
        ...

    @property
    def frame(self) -> int:
        """int: The reading frame (1, 2, or 3)."""
        ...

    @property
    def score(self) -> float:
        """float: The Viterbi score of the gene prediction."""
        ...

    @property
    def insertions(self) -> List[int]:
        """
        A list of 0-based indices representing insertion sequencing errors.
        These are bases that were skipped to maintain the reading frame.
        """
        ...

    @property
    def deletions(self) -> List[int]:
        """
        A list of 0-based indices representing insertion sequencing errors.
        These are bases that were skipped to maintain the reading frame.
        """
        ...

    def sequence(self) -> bytes:
        """Lazily compute the raw nucleotide sequence of the ORF."""
        ...

    def translation(self) -> bytes:
        """Lazily compute the translated amino acid sequence of the ORF."""
        ...


class GeneFinder:
    """The main engine for finding genes, holding the HMM in memory.
    
    Example:
        ```python
        from pyfgs import GeneFinder, Model
        
        finder = GeneFinder(Model.Complete, whole_genome=True)
        genes = finder.find_genes("seq1", "ATGCGTA...")
        ```
    """

    def __init__(self, model: Model, whole_genome: Optional[bool] = None) -> None:
        """
        Initialize the GeneFinder.

        Args:
            model: The sequencing error model to use.
            whole_genome: Set to True if analyzing complete genomic sequences.
                          If None, it defaults to True for Model.Complete,
                          and False for all short-read models.
        """
        ...

    def find_genes(self, header: bytes, sequence: bytes) -> List[Gene]:
        """Predict open reading frames in a given DNA sequence.

        This method releases the GIL, allowing for safe multi-threading
        across multiple CPU cores.

        Args:
            header (bytes): The sequence identifier.
            sequence (bytes): The raw nucleotide buffer.

        Returns:
            List[Gene]: A list of predicted Gene objects.
            
        Example:
            ```python
            genes = finder.find_genes("seq1", "ATGCGTACGTTAG")
            ```
        """
        ...