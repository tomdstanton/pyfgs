# import python.pyfgs as pyfgs
#
#
# def test_read_fasta():
#     header, chromosome = next(pyfgs.FastaReader('tests/data/GCF_001457455.1_NCTC11397_genomic.fna'))
#     finder = pyfgs.GeneFinder(pyfgs.Model.Complete)
#     genes = finder.find_genes(header, chromosome)
