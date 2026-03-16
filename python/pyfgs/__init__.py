# Import the compiled Rust classes from the hidden internal module
from ._pyfgs import Model, FastaReader, FastqReader, Gene, GeneFinder

# Define what gets imported when a user types `from pyfgs import *`
__all__ = ["Model", "FastaReader", "FastqReader", "Gene", "GeneFinder"]