from pyfgs import FastaReader, GeneFinder, Model


def main():

    finder = GeneFinder(pyfgs.Model.Complete)

    reader = FastaReader("reads.fastq")

    for header, seq, qual in reader:
        genes = finder.find_genes(header, seq)