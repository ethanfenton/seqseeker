"""Tests for the aggregation layer."""

import csv
from pathlib import Path
import pytest
from seqseeker._aggregate import aggregate_results


def _write_pair_csv(path: Path, rows: list[dict]) -> None:
    fields = ["sample","lane","barcode","umi","n_kmer_hits",
              "target_pos","strand","read_name","r2_seq"]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def test_aggregate_basic(tmp_path):
    pair_dir = tmp_path / "_per_pair"
    pair_dir.mkdir()
    csv1 = pair_dir / "A1_L001.csv"
    _write_pair_csv(csv1, [
        dict(sample="A1", lane="L001", barcode="AAACCCTTGGAA",
             umi="AACCTTGGAACC", n_kmer_hits=3, target_pos=42,
             strand="+", read_name="read1", r2_seq="ACGT"),
        dict(sample="A1", lane="L001", barcode="AAACCCTTGGAA",
             umi="TTGACCTTGGAA", n_kmer_hits=2, target_pos=100,
             strand="+", read_name="read2", r2_seq="GCTA"),
    ])

    total, n_cells = aggregate_results([str(csv1)], tmp_path)
    assert total == 2
    assert n_cells == 1

    per_cell = list(csv.DictReader(open(tmp_path / "hits_per_cell.csv")))
    assert len(per_cell) == 1
    assert int(per_cell[0]["n_reads"]) == 2
    assert int(per_cell[0]["n_umi"]) == 2


def test_aggregate_empty(tmp_path):
    total, n_cells = aggregate_results([], tmp_path)
    assert total == 0
    assert n_cells == 0


def test_aggregate_missing_file(tmp_path):
    total, n_cells = aggregate_results([str(tmp_path / "nonexistent.csv")], tmp_path)
    assert total == 0
