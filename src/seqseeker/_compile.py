"""Compile the kmer_search C binary, caching next to the source file."""

import shutil
import subprocess
import sys
from pathlib import Path

_C_SRC = Path(__file__).parent / "_c_src" / "kmer_search.c"
_DEFAULT_BIN = Path(__file__).parent / "_c_src" / "kmer_search"


def get_binary(binary_path: Path | None = None, *, force: bool = False) -> Path:
    """Return path to a compiled kmer_search binary, compiling if necessary.

    Parameters
    ----------
    binary_path:
        Where to write the binary. Defaults to the _c_src/ directory inside
        the installed package.
    force:
        Recompile even if the binary is newer than the source.
    """
    src = _C_SRC
    bin_ = Path(binary_path) if binary_path else _DEFAULT_BIN

    if not src.exists():
        raise FileNotFoundError(f"C source not found: {src}")

    if not force and bin_.exists():
        if bin_.stat().st_mtime >= src.stat().st_mtime:
            return bin_

    compiler = shutil.which("gcc") or shutil.which("cc")
    if not compiler:
        raise RuntimeError(
            "No C compiler found (gcc/cc). "
            "Install build-essential (Linux) or Xcode Command Line Tools (macOS)."
        )

    print(f"[seqseeker] Compiling {src.name} → {bin_} ...", file=sys.stderr, flush=True)
    result = subprocess.run(
        [compiler, "-O3", "-march=native", str(src), "-o", str(bin_), "-lz"],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"Compilation failed:\n{result.stderr}")

    print(f"[seqseeker] Compiled OK ({compiler})", file=sys.stderr, flush=True)
    return bin_
