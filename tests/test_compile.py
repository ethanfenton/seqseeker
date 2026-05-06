"""Tests that the C binary compiles cleanly."""

import subprocess
from pathlib import Path
import pytest
from seqseeker._compile import get_binary, _C_SRC


def test_c_source_exists():
    assert _C_SRC.exists(), f"C source not found: {_C_SRC}"


def test_binary_compiles(tmp_path):
    bin_path = tmp_path / "kmer_search"
    binary = get_binary(bin_path, force=True)
    assert binary.exists()
    assert binary.stat().st_size > 0


def test_binary_runs_help(tmp_path):
    bin_path = tmp_path / "kmer_search"
    binary = get_binary(bin_path, force=True)
    result = subprocess.run([str(binary)], capture_output=True, text=True)
    # exits with error code (no args) but should print usage
    assert "Usage" in result.stderr or result.returncode == 1


def test_binary_list_sequences(tmp_path):
    bin_path = tmp_path / "kmer_search"
    binary = get_binary(bin_path, force=True)
    result = subprocess.run([str(binary), "--list-sequences"], capture_output=True, text=True)
    assert result.returncode == 0
    assert "egfp" in result.stderr.lower()
