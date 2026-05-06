"""Command-line interface for seqseeker."""

import argparse
import sys

from . import __version__
from ._compile import get_binary
from ._core import search
from ._sequences import list_sequences


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="seqseeker",
        description=(
            "seqseeker: fast k-mer–based detection of any target sequence "
            "in paired single-cell FASTQ files.\n\n"
            "Sub-commands:\n"
            "  search              Search FASTQ files for a target sequence\n"
            "  compile             Pre-compile the C binary\n"
            "  list-sequences      Show built-in sequences\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--version", action="version", version=f"seqseeker {__version__}")
    sub = p.add_subparsers(dest="command")

    # ── search ────────────────────────────────────────────────────────────
    sp = sub.add_parser("search", help="Search FASTQ files for a target sequence")

    # Input / output
    sp.add_argument(
        "--input", "-i", nargs="+", required=True, metavar="DIR",
        help="One or more directories containing *.fastq.gz pairs",
    )
    sp.add_argument(
        "--output", "-o", default="seqseeker_output", metavar="DIR",
        help="Output directory (default: seqseeker_output)",
    )

    # Target sequence — mutually exclusive group
    seq_grp = sp.add_mutually_exclusive_group(required=True)
    seq_grp.add_argument(
        "--sequence", "-q", metavar="ACGT…",
        help="Raw nucleotide string to search for",
    )
    seq_grp.add_argument(
        "--sequence-file", "-f", metavar="FILE",
        help="FASTA file containing the target sequence",
    )
    seq_grp.add_argument(
        "--sequence-name", "-n", metavar="NAME",
        help="Built-in sequence name (e.g. egfp, mcherry). "
             "Run `seqseeker list-sequences` to see all options.",
    )

    # Region
    sp.add_argument(
        "--region", "-R", choices=["full", "5prime", "3prime"], default="full",
        help="Which part of the target to index (default: full)",
    )
    sp.add_argument(
        "--region-len", "-L", type=int, default=200, metavar="N",
        help="Length in bp when --region is 5prime or 3prime (default: 200)",
    )

    # k-mer / scoring
    sp.add_argument("--kmer-size",   "-k", type=int, default=20, metavar="K",
                    help="k-mer length (default 20; range 8–31)")
    sp.add_argument("--min-kmers",   "-m", type=int, default=1,  metavar="N",
                    help="Minimum k-mer hits to call a read positive (default 1)")

    # 10x read layout
    sp.add_argument("--barcode-len", "-b", type=int, default=16, metavar="N",
                    help="Cell barcode length in R1 (default 16, 10x v3)")
    sp.add_argument("--umi-len",     "-u", type=int, default=12, metavar="N",
                    help="UMI length in R1 (default 12, 10x v3)")

    # Parallelism
    sp.add_argument(
        "--workers", "-w", type=int, default=1, metavar="N",
        help="Parallel worker processes (default 1; 0 = auto-detect)",
    )

    # Filtering / misc
    sp.add_argument(
        "--sample-filter", nargs="*", metavar="SAMPLE",
        help="Process only these sample IDs (directory names)",
    )
    sp.add_argument(
        "--report-interval", type=int, default=1_000_000, metavar="N",
        help="Print progress every N reads (default 1 000 000)",
    )
    sp.add_argument(
        "--binary", metavar="PATH",
        help="Path to a pre-compiled kmer_search binary",
    )

    # ── compile ───────────────────────────────────────────────────────────
    cp = sub.add_parser("compile", help="Pre-compile the C search binary")
    cp.add_argument("--binary", metavar="PATH",
                    help="Destination path for the compiled binary")
    cp.add_argument("--force", action="store_true", help="Recompile even if up to date")

    # ── list-sequences ────────────────────────────────────────────────────
    sub.add_parser("list-sequences", help="Show built-in sequences")

    return p


def main(argv: list[str] | None = None) -> None:
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        sys.exit(0)

    if args.command == "list-sequences":
        list_sequences()
        return

    if args.command == "compile":
        binary = get_binary(
            getattr(args, "binary", None),
            force=getattr(args, "force", False),
        )
        print(f"Binary ready: {binary}")
        return

    if args.command == "search":
        search(
            input_dirs=args.input,
            output_dir=args.output,
            sequence=args.sequence,
            sequence_file=args.sequence_file,
            sequence_name=args.sequence_name,
            region=args.region,
            region_len=args.region_len,
            kmer_size=args.kmer_size,
            barcode_len=args.barcode_len,
            umi_len=args.umi_len,
            min_kmers=args.min_kmers,
            workers=args.workers,
            sample_filter=args.sample_filter,
            report_interval=args.report_interval,
            binary_path=args.binary,
        )
        return

    parser.print_help()


if __name__ == "__main__":
    main()
