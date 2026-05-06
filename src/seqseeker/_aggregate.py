"""Aggregate per-pair CSV outputs into summary tables."""

import csv
from collections import defaultdict
from pathlib import Path


def aggregate_results(csv_paths: list[str | Path], out_dir: Path) -> tuple[int, int]:
    """Concatenate per-pair CSVs and write summary tables.

    Returns
    -------
    (total_rows, n_unique_barcodes)
    """
    out_dir = Path(out_dir)
    raw_path      = out_dir / "hits_raw.csv"
    per_cell_path = out_dir / "hits_per_cell.csv"
    pos_path      = out_dir / "kmer_positions.csv"

    raw_fields = [
        "sample", "lane", "barcode", "umi",
        "n_kmer_hits", "target_pos", "strand", "read_name", "r2_seq",
    ]

    cell_data  = defaultdict(lambda: dict(n_reads=0, umis=set(), samples=set()))
    pos_counts = defaultdict(int)
    total_rows = 0

    with open(raw_path, "w", newline="") as raw_fh:
        writer = csv.DictWriter(raw_fh, fieldnames=raw_fields)
        writer.writeheader()

        for csv_path in csv_paths:
            p = Path(csv_path)
            if not p.exists() or p.stat().st_size == 0:
                continue
            with open(p, newline="") as fh:
                for row in csv.DictReader(fh):
                    writer.writerow({f: row.get(f, "") for f in raw_fields})
                    total_rows += 1
                    bc = row.get("barcode", "")
                    cell_data[bc]["n_reads"] += 1
                    cell_data[bc]["umis"].add(row.get("umi", ""))
                    cell_data[bc]["samples"].add(row.get("sample", ""))
                    try:
                        pos_counts[int(row.get("target_pos", -1))] += 1
                    except ValueError:
                        pass

    per_cell = sorted(
        [
            dict(
                barcode=bc,
                n_reads=d["n_reads"],
                n_umi=len(d["umis"]),
                samples=";".join(sorted(d["samples"])),
            )
            for bc, d in cell_data.items()
        ],
        key=lambda x: -x["n_reads"],
    )

    with open(per_cell_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["barcode", "n_reads", "n_umi", "samples"])
        w.writeheader()
        w.writerows(per_cell)

    with open(pos_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["target_pos", "n_reads"])
        for pos in sorted(pos_counts):
            w.writerow([pos, pos_counts[pos]])

    return total_rows, len(per_cell)
