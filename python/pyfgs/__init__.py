# Import the compiled Rust classes from the hidden internal module
from ._pyfgs import (
    BedWriter,
    FaaWriter,
    FastaReader,
    FastqReader,
    FnaWriter,
    Gene,
    GeneFinder,
    Gff3Writer,
    Model,
    Mutation,
    VcfWriter,
)

# Define what gets imported when a user types `from pyfgs import *`
__all__ = [
    "Model",
    "FastaReader",
    "FastqReader",
    "Gene",
    "GeneFinder",
    "Mutation",
    "BedWriter",
    "VcfWriter",
    "Gff3Writer",
    "FnaWriter",
    "FaaWriter",
]
