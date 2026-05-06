"""Core Python orchestration layer for seqseeker."""

import multiprocessing
import os
import re
import subprocess
import time
from datetime import datetime, timedelta
from pathlib import Path

from ._aggregate import aggregate_results
from ._compile import get_binary
from ._sequences import get_sequence, BUILTIN_SEQUENCES


def _find_fastq_pairs(
    input_dirs: list[str | Path],
    sample_filter: list[str] | None = None,
) -> list[dict]:
    """Discover R1/R2 paired FASTQ files under input_dirs."""
    r1_files: dict = {}
    r2_files: dict = {}

    for base in input_dirs:
        base = Path(base)
        if not base.exists():
            print(f"[seqseeker WARN] Directory not found: {base}", flush=True)
            continue
        for f in sorted(base.rglob("*.fastq.gz")):
            if re.search(r"_I[12]_", f.name):
                continue
            key = (str(f.parent), re.sub(r"_R[12]_", "_RX_", f.name))
            if re.search(r"_R1_", f.name):
                r1_files[key] = f
            elif re.search(r"_R2_", f.name):
                r2_files[key] = f

    pairs = []
    for key, r1 in r1_files.items():
        r2 = r2_files.get(key)
        if r2 is None:
            print(f"[seqseeker WARN] No R2 for {r1.name}", flush=True)
            continue
        sample = r1.parent.name
        if sample_filter and sample not in sample_filter:
            continue
        m = re.search(r"_([A-Z0-9]+_S\d+_L\d+)_R1_", r1.name)
        lane = m.group(1) if m else r1.stem
        pairs.append(dict(sample=sample, lane=lane, r1=str(r1), r2=str(r2)))

    pairs.sort(key=lambda x: (x["sample"], x["lane"]))
    return pairs


def _run_pair(args: dict) -> dict:
    """Run kmer_search binary for one R1/R2 pair (worker function)."""
    cmd = [
        args["binary"],
        "-1", args["pair"]["r1"],
        "-2", args["pair"]["r2"],
        "-s", args["pair"]["sample"],
        "-l", args["pair"]["lane"],
        "-o", str(args["out_csv"]),
        "-k", str(args["kmer_size"]),
        "-b", str(args["barcode_len"]),
        "-u", str(args["umi_len"]),
        "-m", str(args["min_kmers"]),
        "-r", str(args["report_interval"]),
        "-R", args["region"],
    ]
    if args["region"] != "full":
        cmd += ["-L", str(args["region_len"])]

    # Sequence source
    if args.get("sequence_raw"):
        cmd += ["-q", args["sequence_raw"]]
    elif args.get("sequence_file"):
        cmd += ["-f", args["sequence_file"]]
    elif args.get("sequence_name"):
        cmd += ["-n", args["sequence_name"]]

    proc = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    for line in proc.stderr:
        print(line, end="", flush=True)
    proc.wait()

    return dict(
        sample=args["pair"]["sample"],
        lane=args["pair"]["lane"],
        out_csv=str(args["out_csv"]),
        returncode=proc.returncode,
    )


def search(
    input_dirs: str | list[str | Path],
    output_dir: str | Path = "seqseeker_output",
    *,
    sequence: str | None = None,
    sequence_file: str | Path | None = None,
    sequence_name: str | None = None,
    region: str = "full",
    region_len: int = 200,
    kmer_size: int = 20,
    barcode_len: int = 16,
    umi_len: int = 12,
    min_kmers: int = 1,
    workers: int = 1,
    sample_filter: list[str] | None = None,
    report_interval: int = 1_000_000,
    binary_path: str | Path | None = None,
) -> dict:
    """Search FASTQ files for reads containing a target sequence.

    Parameters
    ----------
    input_dirs:
        One or more directories to search recursively for *.fastq.gz pairs.
    output_dir:
        Directory where output CSVs are written.
    sequence:
        Raw nucleotide string to search for.
    sequence_file:
        Path to a FASTA file containing the target sequence.
    sequence_name:
        Name of a built-in sequence (e.g. ``"egfp"``, ``"mcherry"``).
        Exactly one of sequence/sequence_file/sequence_name must be given.
    region:
        Which part of the target to index: ``"full"``, ``"5prime"``, or
        ``"3prime"``.
    region_len:
        Length in bp for 5prime/3prime region slicing (default 200).
    kmer_size:
        k-mer length for matching (default 20; lower = more sensitive,
        higher = more specific).
    barcode_len:
        Cell barcode length in R1 (default 16 for 10x v3).
    umi_len:
        UMI length in R1 (default 12 for 10x v3).
    min_kmers:
        Minimum k-mer hits to call a read as positive (default 1).
    workers:
        Number of parallel processes. 0 = auto (cpu_count − 1), 1 = serial.
    sample_filter:
        Only process samples whose directory name is in this list.
    report_interval:
        Print progress every N reads.
    binary_path:
        Override path to a pre-compiled kmer_search binary.

    Returns
    -------
    dict with keys: output_dir, raw_csv, per_cell_csv, positions_csv,
    n_hits, n_unique_barcodes.
    """
    # Validate sequence source
    sources = sum([
        sequence is not None,
        sequence_file is not None,
        sequence_name is not None,
    ])
    if sources != 1:
        raise ValueError(
            "Exactly one of sequence=, sequence_file=, or sequence_name= must be provided."
        )

    if region not in ("full", "5prime", "3prime"):
        raise ValueError("region must be 'full', '5prime', or '3prime'.")

    if isinstance(input_dirs, (str, Path)):
        input_dirs = [input_dirs]

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir = out_dir / "_per_pair"
    tmp_dir.mkdir(exist_ok=True)

    binary = str(get_binary(binary_path))

    pairs = _find_fastq_pairs(input_dirs, sample_filter)
    if not pairs:
        raise RuntimeError("No FASTQ pairs found in the specified directories.")

    n_cpu = os.cpu_count() or 1
    if workers == 0:
        n_workers = min(len(pairs), max(1, n_cpu - 1))
    else:
        n_workers = workers

    print("=" * 70, flush=True)
    print(f"  seqseeker", flush=True)
    print(f"  Started   : {datetime.now()}", flush=True)
    print(f"  CPUs      : {n_cpu}   workers: {n_workers}", flush=True)
    print(f"  Pairs     : {len(pairs)} across "
          f"{len({p['sample'] for p in pairs})} sample(s)", flush=True)
    print(f"  k-mer     : {kmer_size}  region: {region}"
          + (f"({region_len} bp)" if region != "full" else ""), flush=True)
    print(f"  Min hits  : {min_kmers}", flush=True)
    print("=" * 70, flush=True)

    worker_args = []
    csv_paths   = []
    for p in pairs:
        out_csv = tmp_dir / f"{p['sample']}_{p['lane']}.csv"
        csv_paths.append(str(out_csv))
        worker_args.append(dict(
            binary=binary,
            pair=p,
            out_csv=out_csv,
            kmer_size=kmer_size,
            barcode_len=barcode_len,
            umi_len=umi_len,
            min_kmers=min_kmers,
            report_interval=report_interval,
            region=region,
            region_len=region_len,
            sequence_raw=sequence,
            sequence_file=str(sequence_file) if sequence_file else None,
            sequence_name=sequence_name,
        ))

    t0 = time.time()
    if n_workers == 1:
        results = [_run_pair(wa) for wa in worker_args]
    else:
        with multiprocessing.Pool(processes=n_workers) as pool:
            results = list(pool.imap_unordered(_run_pair, worker_args))
            for r in results:
                if r["returncode"] != 0:
                    print(f"[WARN] {r['sample']}/{r['lane']} exited with "
                          f"code {r['returncode']}", flush=True)

    elapsed = time.time() - t0

    print(f"\n[seqseeker] Aggregating results ...", flush=True)
    total_hits, n_cells = aggregate_results(csv_paths, out_dir)

    print("=" * 70, flush=True)
    print(f"  Wall time    : {str(timedelta(seconds=int(elapsed)))}", flush=True)
    print(f"  Hits         : {total_hits:,}", flush=True)
    print(f"  Unique cells : {n_cells:,}", flush=True)
    print(f"  Raw CSV      : {out_dir}/hits_raw.csv", flush=True)
    print(f"  Per-cell     : {out_dir}/hits_per_cell.csv", flush=True)
    print(f"  k-mer pos    : {out_dir}/kmer_positions.csv", flush=True)
    print("=" * 70, flush=True)
    print(f"""
  R snippet (Seurat):
    hits <- read.csv("{out_dir}/hits_per_cell.csv")
    hits$bc <- paste0(hits$barcode, "-1")
    seurat_obj$seq_hits <- hits$n_reads[match(colnames(seurat_obj), hits$bc)]
    seurat_obj$seq_hits[is.na(seurat_obj$seq_hits)] <- 0
""", flush=True)

    return dict(
        output_dir=str(out_dir),
        raw_csv=str(out_dir / "hits_raw.csv"),
        per_cell_csv=str(out_dir / "hits_per_cell.csv"),
        positions_csv=str(out_dir / "kmer_positions.csv"),
        n_hits=total_hits,
        n_unique_barcodes=n_cells,
    )
