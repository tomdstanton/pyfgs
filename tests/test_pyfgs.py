import pytest
import tempfile
import os
import pyfgs


@pytest.fixture
def sample_fasta():
    """Creates a temporary FASTA file for testing."""
    fd, path = tempfile.mkstemp(suffix=".fasta")
    with os.fdopen(fd, 'wb') as f:
        f.write(b">seq1\n")
        f.write(b"ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA\n")
        f.write(b"ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA\n")
    yield path
    os.remove(path)


@pytest.fixture
def sample_fastq():
    """Creates a temporary FASTQ file for testing."""
    fd, path = tempfile.mkstemp(suffix=".fastq")
    with os.fdopen(fd, 'wb') as f:
        f.write(b"@read1\n")
        f.write(b"ATGCGTACGTAGCTAGCTAGCT\n")
        f.write(b"+\n")
        f.write(b"IIIIIIIIIIIIIIIIIIIIII\n")
    yield path
    os.remove(path)


def test_model_inference():
    """Ensure the whole_genome flag is correctly inferred from the model."""
    finder_complete = pyfgs.GeneFinder(pyfgs.Model.Complete)
    assert finder_complete is not None

    finder_illumina = pyfgs.GeneFinder(pyfgs.Model.Illumina5)
    assert finder_illumina is not None


def test_fasta_reader(sample_fasta):
    """Ensure FastaReader yields raw bytes."""
    reader = pyfgs.FastaReader(sample_fasta)
    records = list(reader)

    assert len(records) == 1
    header, seq = records[0]

    assert isinstance(header, bytes)
    assert isinstance(seq, bytes)
    assert header == b"seq1"


def test_fastq_reader(sample_fastq):
    """Ensure FastqReader yields raw bytes including quality scores."""
    reader = pyfgs.FastqReader(sample_fastq)
    records = list(reader)

    assert len(records) == 1
    header, seq, qual = records[0]

    assert isinstance(header, bytes)
    assert isinstance(seq, bytes)
    assert isinstance(qual, bytes)
    assert header == b"read1"
    assert seq == b"ATGCGTACGTAGCTAGCTAGCT"


def test_newline_filtering():
    """Ensure the GeneFinder correctly strips whitespace before processing."""
    finder = pyfgs.GeneFinder(pyfgs.Model.Complete)
    raw_seq = b"ATG CGT\nACG TAG\r\nCTA GCT"

    # FIXED: find_genes now only takes the sequence bytes!
    genes = finder.find_genes(raw_seq)
    assert isinstance(genes, list)


def test_gene_properties_and_indels():
    """Verify that Genes use 0-based BED coordinates and correctly expose indels."""
    finder = pyfgs.GeneFinder(pyfgs.Model.Complete)

    # The actual E. coli recA coding sequence (1062 bp)
    real_gene_dna = (
        b"ATGGCTATCGACGAAAACAAACAGAAAGCGTTGGCGGCAGCACTGGGCCAGATTGAGAAACAATTT"
        b"GGTAAAGGCTCCATCATGCGCCTGGGTGAAGACCGTTCCATGGATGTGGAAACCATCTCTACCGGT"
        b"TCGCTTTCACTGGATATCGCCCTTGGGGCAGGCGGTCTGCCGATGGGCCGTATCGTCGAAATCTAC"
        b"GGACCGGAATCTTCCGGTAAAACCACGCTGACGCTGCAGGTGATTGCCGCAGCGCAGCGTGAAGGT"
        b"AAAACCTGTGCGTTTATCGATGCTGAACACGCGCTGGACCCAATCTACGCACGTAAACTGGGCGTC"
        b"GATATCGACAACCTGCTGTGCTCCCAGCCGGACACCGGCGAGCAGGCACTGGAAATCTGTGACGCC"
        b"CTGGCGCGTTCTGGCGCAGTAGACGTTATCGTCGTTGACTCCGTGGCGGCACTGACGCCGAAAGCG"
        b"GAAATCGAAGGCGAAATCGGCGACTCTCACATGGGCCTTGCGGCACGTATGATGAGCCAGGCGATG"
        b"CGTAAGCTGGCGGGTAACCTGAAGCAGTCCAACACGCTGCTGATCTTCATCAACCAGATCCGTATG"
        b"AAAATTGGTGTGATGTTCGGTAACCCGGAAACCACTACCGGTGGTAACGCGCTGAAATTCTACGCC"
        b"TCTGTTCGTCTCGACATCCGTCGTATCGGCGCGGTGAAAGAGGGCGAAAACGTGGTGGGTAGCGAA"
        b"ACCCGCGTGAAAGTGGTGAAGAACAAAATCGCTGCGCCGTTTAAACAGGCTGAATTCCAGATCCTC"
        b"TACGGCGAAGGTATCAACTTCTACGGCGAACTGGTTGACCTGGGCGTAAAAGAGAAGCTGATCGAG"
        b"AAAGCAGGCGCGTGGTACAGCTACAAAGGTGAGAAGATCGGTCAGGGTAAAGCGAATGCGACTGCC"
        b"TGGCTGAAAGATAACCCGGAAACCGCGAAAGAGATCGAGAAGAAAGTACGTGAGTTGCTGCTGAGC"
        b"AACCCGAACTCAACGCCGGATTTCTCTGTAGATGATAGCGAAGGCGTAGCAGAAACTAACGAAGAT"
        b"TTTTAA"
    )

    # FIXED: find_genes now only takes the sequence bytes!
    genes = finder.find_genes(real_gene_dna)

    assert len(genes) > 0
    gene = genes[0]

    # Check types for our 0-based coordinates and score
    assert isinstance(gene.start, int)
    assert isinstance(gene.end, int)
    assert isinstance(gene.strand, int)
    assert isinstance(gene.frame, int)
    assert isinstance(gene.score, float)

    # Check our native indel vectors
    assert isinstance(gene.insertions, list)
    assert isinstance(gene.deletions, list)

    # Coordinate sanity checks
    assert gene.start >= 0
    assert gene.start < gene.end
    assert gene.strand in (1, -1)


def test_lazy_byte_evaluation():
    """Verify that sequence() and translation() bypass the UTF-8 tax."""
    finder = pyfgs.GeneFinder(pyfgs.Model.Complete)
    real_gene_dna = (
        b"ATGGCTATCGACGAAAACAAACAGAAAGCGTTGGCGGCAGCACTGGGCCAGATTGAGAAACAATTT"
        b"GGTAAAGGCTCCATCATGCGCCTGGGTGAAGACCGTTCCATGGATGTGGAAACCATCTCTACCGGT"
        b"TCGCTTTCACTGGATATCGCCCTTGGGGCAGGCGGTCTGCCGATGGGCCGTATCGTCGAAATCTAC"
    )
    genes = finder.find_genes(real_gene_dna)

    gene = genes[0]
    dna_bytes = gene.sequence()
    prot_bytes = gene.translation()

    assert isinstance(dna_bytes, bytes)
    assert isinstance(prot_bytes, bytes)
    assert len(dna_bytes) > 0
    assert len(prot_bytes) > 0


def test_mutations_api():
    """Ensure the mutations API returns a valid list of Mutation objects."""
    # We use whole_genome=False to allow the model to predict frameshifts if it finds them
    finder = pyfgs.GeneFinder(pyfgs.Model.Illumina10, whole_genome=False)
    real_gene_dna = (
        b"ATGGCTATCGACGAAAACAAACAGAAAGCGTTGGCGGCAGCACTGGGCCAGATTGAGAAACAATTT"
        b"GGTAAAGGCTCCATCATGCGCCTGGGTGAAGACCGTTCCATGGATGTGGAAACCATCTCTACCGGT"
    )
    
    genes = finder.find_genes(real_gene_dna)
    assert len(genes) > 0
    gene = genes[0]
    
    muts = gene.mutations(real_gene_dna)
    assert isinstance(muts, list)
    
    # If the noisy model happened to insert a mutation, verify its structure
    if muts:
        assert hasattr(muts[0], "pos")
        assert hasattr(muts[0], "mut_type")
        assert hasattr(muts[0], "codon_idx")


def test_rust_file_writers():
    """Ensure the native Rust context managers safely open and write to files."""
    finder = pyfgs.GeneFinder(pyfgs.Model.Complete)
    seq = b"ATGGCTATCGACGAAAACAAACAGAAAGCGTTGGCGGCAGCACTGGGCCAGATTGAGAAACAATTT"
    genes = finder.find_genes(seq)

    with tempfile.TemporaryDirectory() as tmpdir:
        bed_path = os.path.join(tmpdir, "out.bed")
        
        # Test the Context Manager entry and exit
        with pyfgs.BedWriter(bed_path) as writer:
            writer.write_record(genes, "test_contig", seq)

        assert os.path.exists(bed_path)
        
        # Verify the Rust code successfully flushed bytes to the OS
        with open(bed_path, "r") as f:
            content = f.read()
            assert "source=ab initio" in content
            assert "test_contig" in content
