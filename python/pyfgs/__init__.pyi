from typing import List, Iterator, Tuple


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

    def __iter__(self) -> Iterator[Tuple[str, str]]:
        ...

    def __next__(self) -> Tuple[str, str]: ...



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

    def __iter__(self) -> Iterator[Tuple[str, str, str]]:
        ...

    def __next__(self) -> Tuple[str, str, str]: ...


class Gene:
    """Represents a single predicted Open Reading Frame (ORF).
    
    Example:
        ```python
        for gene in genes:
            print(f"Start: {gene.start}, End: {gene.end}, Strand: {gene.strand}, Frame: {gene.frame}, Score: {gene.score}")
        ```
    """

    @property
    def start(self) -> int:
        """int: The 1-based start position of the gene."""
        ...

    @property
    def end(self) -> int:
        """int: The 1-based end position of the gene."""
        ...

    @property
    def strand(self) -> str:
        """str: The strand of the gene ('+' or '-')."""
        ...

    @property
    def frame(self) -> int:
        """int: The reading frame (1, 2, or 3)."""
        ...

    @property
    def score(self) -> float:
        """float: The Viterbi score of the gene prediction."""
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

    def __init__(self, model: Model, whole_genome: bool = False) -> None:
        """Initialize the GeneFinder.

        Args:
            model (Model): The sequencing error model to use (e.g., pyfgs.Model.Illumina5).
            whole_genome (bool, optional): Set to True if analyzing complete genomic sequences
                rather than short reads. Defaults to False.
        """
        ...

    def find_genes(self, header: str, sequence: str) -> List[Gene]:
        """Predict open reading frames in a given DNA sequence.

        This method releases the GIL, allowing for safe multi-threading
        across multiple CPU cores.

        Args:
            header (str): The sequence identifier.
            sequence (str): The raw nucleotide string.

        Returns:
            List[Gene]: A list of predicted Gene objects.
            
        Example:
            ```python
            genes = finder.find_genes("seq1", "ATGCGTACGTTAG")
            ```
        """
        ...