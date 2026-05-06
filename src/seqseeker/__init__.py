"""
seqseeker — fast k-mer–based sequence detection in paired FASTQ files.

Designed for single-cell RNA-seq (10x Chromium) but works with any paired
FASTQ data where reads from R1 carry a cell barcode + UMI and R2 carries
the insert.
"""

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("seqseeker")
except PackageNotFoundError:
    __version__ = "0.0.0"

from ._core import search
from ._sequences import BUILTIN_SEQUENCES, list_sequences

__all__ = ["search", "BUILTIN_SEQUENCES", "list_sequences", "__version__"]
