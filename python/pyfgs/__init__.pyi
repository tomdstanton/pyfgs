# pyfgs.pyi
"""
Python bindings for FragGeneScanRs.

This module provides a high-performance, memory-efficient ab initio gene
predictor and frameshift annotator for short-read and complete genome assemblies.
"""

from typing import List, Tuple, Optional, Iterator, Type, Any
from types import TracebackType


class Model:
    """
    The available sequencing error models for the Hidden Markov Model.

    The model alters the transition probabilities to account for expected
    sequencing error rates, making it more or less forgiving of frameshifts.

    Attributes: 
        Illumina1: Illumina reads with ~0.1% error rate.
        Illumina5: Illumina reads with ~0.5% error rate.
        Illumina10: Illumina reads with ~1% error rate.
        Sanger5: Sanger reads with ~0.5% error rate.
        Sanger10: Sanger reads with ~1% error rate.
        Pyro454_5: 454 pyrosequencing reads with ~0.5% error rate.
        Pyro454_10: 454 pyrosequencing reads with ~1% error rate.
        Pyro454_30: 454 pyrosequencing reads with ~3% error rate.
        Complete: Complete genomic sequences without expected sequencing errors.

    Examples: 
        >>> import pyfgs
        >>> model = pyfgs.Model.Complete

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


class FastaReader(Iterator[Tuple[bytes, bytes]]):
    """
    A fast, memory-efficient FASTA parser implemented in Rust.

    Args: 
        path (str): The file path to the FASTA file.

    Yields: 
        Tuple[bytes, bytes]: A tuple containing the header and sequence as raw bytes.

    Examples: 
        >>> from pyfgs import FastaReader
        >>> reader = FastaReader("genome.fna")
        >>> for header, sequence in reader:
        ...     print(f"ID: {header.decode()} | Length: {len(sequence)}")
    """

    def __init__(self, path: str) -> None: ...

    def __iter__(self) -> 'FastaReader': ...

    def __next__(self) -> Tuple[bytes, bytes]: ...


class FastqReader(Iterator[Tuple[bytes, bytes, bytes]]):
    """
    A fast, memory-efficient FASTQ parser implemented in Rust.

    Args: 
        path (str): The file path to the FASTQ file.

    Yields: 
        Tuple[bytes, bytes, bytes]: A tuple containing the header, sequence,
        and Phred quality string as raw bytes.
    """

    def __init__(self, path: str) -> None: ...

    def __iter__(self) -> 'FastqReader': ...

    def __next__(self) -> Tuple[bytes, bytes, bytes]: ...


class Mutation:
    """
    Represents a frameshift mutation (insertion or deletion) detected by the HMM.

    Attributes: 
        pos (int): The 0-based index of the mutation in the global assembly.
                   (Note: This is mathematically identical to the 1-based VCF anchor position).
        mut_type (str): Either 'ins' (extra base in assembly) or 'del' (missing base).
        ref_allele (str): The reference allele from the raw assembly.
        alt_allele (str): The conceptually corrected allele determined by the model.
        codon_idx (int): The 1-based codon index where the reading frame breaks.
        annotation (str): A Snippy-style text annotation (e.g., 'AA:X->fs|DNA:ATC->AT').
    """

    @property
    def pos(self) -> int: ...

    @property
    def mut_type(self) -> str: ...

    @property
    def ref_allele(self) -> str: ...

    @property
    def alt_allele(self) -> str: ...

    @property
    def codon_idx(self) -> int: ...

    @property
    def annotation(self) -> str: ...


class Gene:
    """
    Represents a single predicted Open Reading Frame (ORF).

    Attributes: 
        start (int): The 0-based, inclusive start coordinate.
        end (int): The 0-based, exclusive end coordinate.
        strand (int): The strand of the feature (1 for forward, -1 for reverse).
        frame (int): The reading frame.
        score (float): The log-probability score of the HMM prediction.
        insertions (List[int]): 1-based global coordinates of predicted insertions.
        deletions (List[int]): 1-based global coordinates of predicted deletions.
    """

    @property
    def start(self) -> int: ...

    @property
    def end(self) -> int: ...

    @property
    def strand(self) -> int: ...

    @property
    def frame(self) -> int: ...

    @property
    def score(self) -> float: ...

    @property
    def insertions(self) -> List[int]: ...

    @property
    def deletions(self) -> List[int]: ...

    def sequence(self) -> bytes: 
        """
        Retrieves the raw nucleotide sequence of the predicted gene.

        If the gene is on the reverse strand, the sequence is automatically
        reverse-complemented. If the gene contains frameshifts, the returned
        sequence represents the conceptual (corrected) HMM path.

        Returns: 
            bytes: The DNA sequence.
        """
        ...

    def translation(self) -> bytes: 
        """
        Translates the predicted gene into an amino acid sequence.

        Alternative start codons (e.g., GTG, TTG) are automatically translated
        to Methionine (M) if the model was initialized with `whole_genome=True`.

        Returns: 
            bytes: The amino acid sequence.
        """
        ...

    def mutations(self, sequence: bytes) -> List[Mutation]:
        """
        Extracts structural variant objects for any predicted frameshifts.

        Args: 
            sequence (bytes): The raw parent contig sequence, used to determine
                VCF anchored alleles.

        Returns: 
            List[Mutation]: A list of structured mutation objects.

        Examples: 
            >>> seq_bytes = str(record.seq).encode()
            >>> for gene in genes: 
            ...     for mut in gene.mutations(seq_bytes):
            ...         print(f"Frameshift at codon {mut.codon_idx}: {mut.annotation}")
        """
        ...


class GeneFinder:
    """
    The core ab initio gene prediction engine.

    Args: 
        model (Model): The sequencing error model to use.
        whole_genome (bool, optional): If False, the HMM permits internal
            frameshifts (insertions/deletions) typical of sequencing errors
            or pseudogenes. If True, strictly enforces contiguous reading frames.
            Defaults to True if `Model.Complete` is used, otherwise False.

    Examples: 
        >>> from Bio import SeqIO
        >>> import pyfgs
        >>> record = SeqIO.read("genome.fasta", "fasta")
        >>>
        >>> # Keep whole_genome=False to allow pseudogene/frameshift detection
        >>> finder = pyfgs.GeneFinder(pyfgs.Model.Complete, whole_genome=False)
        >>> genes = finder.find_genes(record.seq._data)
    """

    def __init__(self, model: Model, whole_genome: Optional[bool] = None) -> None: ...

    def find_genes(self, sequence: bytes) -> List[Gene]:
        """
        Predicts open reading frames in a given DNA sequence.

        This method releases the Python GIL, allowing for safe, lock-free
        multi-threading across multiple CPU cores.

        Args: 
            sequence (bytes): The raw nucleotide sequence.

        Returns: 
            List[Gene]: A list of predicted Gene objects.
        """
        ...


class VcfWriter:
    """
    A high-performance streaming context manager for writing VCF v4.2 files.

    Automatically translates structural frameshifts detected by the HMM into
    anchored VCF variants with SnpEff-compliant `ANN` fields.

    Args: 
        output_path (str): The destination file path.

    Examples: 
        >>> with pyfgs.VcfWriter("variants.vcf") as vcf:
        ...     vcf.write_record(genes, record.id, str(record.seq).encode())
    """

    def __init__(self, output_path: str) -> None: ...

    def __enter__(self) -> 'VcfWriter': ...

    def __exit__(self, exc_type: Optional[Type[BaseException]], exc_val: Optional[BaseException],
                 exc_tb: Optional[TracebackType]) -> None: ...

    def write_record(self, genes: List[Gene], header: str, sequence: bytes) -> None:
        """
        Writes the frameshift variants for a single contig to the VCF buffer.

        Args: 
            genes (List[Gene]): The list of predicted Gene objects.
            header (str): The chromosome/contig ID (used for the #CHROM column).
            sequence (bytes): The raw parent nucleotide sequence.
        """
        ...


class BedWriter:
    """
    A high-performance streaming context manager for writing Extended BED files.

    Outputs a BED6+1 format, where the 7th column contains a VCF-style INFO
    string detailing any frameshifts present in the feature.

    Args: 
        output_path (str): The destination file path.
    """

    def __init__(self, output_path: str) -> None: ...

    def __enter__(self) -> 'BedWriter': ...

    def __exit__(self, exc_type: Optional[Type[BaseException]], exc_val: Optional[BaseException],
                 exc_tb: Optional[TracebackType]) -> None: ...

    def write_record(self, genes: List[Gene], header: str, sequence: bytes) -> None:
        """
        Writes the BED intervals for a single contig to the buffer.

        Args: 
            genes (List[Gene]): The list of predicted Gene objects.
            header (str): The chromosome/contig ID.
            sequence (bytes): The raw parent nucleotide sequence.
        """
        ...


class Gff3Writer:
    """
    A high-performance streaming context manager for writing INSDC-compliant GFF3.

    Automatically shifts 0-based coordinates to 1-based fully-closed coordinates.
    Genes containing frameshifts are flagged as `pseudogene=unknown` to ensure
    compliance with downstream translation tools (like Prokka or Bakta).

    Args: 
        output_path (str): The destination file path.
    """

    def __init__(self, output_path: str) -> None: ...

    def __enter__(self) -> 'Gff3Writer': ...

    def __exit__(self, exc_type: Optional[Type[BaseException]], exc_val: Optional[BaseException],
                 exc_tb: Optional[TracebackType]) -> None: ...

    def write_record(self, genes: List[Gene], header: str, sequence: bytes) -> None:
        """
        Writes the GFF3 annotations for a single contig to the buffer.

        Args: 
            genes (List[Gene]): The list of predicted Gene objects.
            header (str): The chromosome/contig ID.
            sequence (bytes): The raw parent nucleotide sequence.
        """
        ...


class FnaWriter:
    """
    A high-performance streaming context manager for writing nucleotide FASTA files.

    Outputs raw, non-wrapped byte streams for maximum parsing speed by downstream tools.

    Args: 
        output_path (str): The destination file path.
    """

    def __init__(self, output_path: str) -> None: ...

    def __enter__(self) -> 'FnaWriter': ...

    def __exit__(self, exc_type: Optional[Type[BaseException]], exc_val: Optional[BaseException],
                 exc_tb: Optional[TracebackType]) -> None: ...

    def write_record(self, genes: List[Gene], header: str) -> None:
        """
        Writes the conceptual nucleotide sequences for a single contig.

        Args: 
            genes (List[Gene]): The list of predicted Gene objects.
            header (str): The chromosome/contig ID.
        """
        ...


class FaaWriter:
    """
    A high-performance streaming context manager for writing amino acid FASTA files.

    Outputs raw, non-wrapped byte streams of the translated proteins.

    Args: 
        output_path (str): The destination file path.
    """

    def __init__(self, output_path: str) -> None: ...

    def __enter__(self) -> 'FaaWriter': ...

    def __exit__(self, exc_type: Optional[Type[BaseException]], exc_val: Optional[BaseException],
                 exc_tb: Optional[TracebackType]) -> None: ...

    def write_record(self, genes: List[Gene], header: str) -> None:
        """
        Writes the translated amino acid sequences for a single contig.

        Args: 
            genes (List[Gene]): The list of predicted Gene objects.
            header (str): The chromosome/contig ID.
        """
        ...
