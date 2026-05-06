# seqseeker

**Fast k-mer–based detection of any target sequence in paired single-cell FASTQ files**

`seqseeker` finds reads containing a user-specified nucleotide sequence (GFP,
mCherry, Cre, any transgene) directly in raw gzipped FASTQ files from 10x
Chromium experiments, attributes each hit to a cell barcode + UMI, and writes
output tables ready for Seurat.

No custom genome build. No Cell Ranger re-run. Per-sample runtimes of
**5–15 minutes** vs. hours for alignment.

[![CI](https://github.com/ethanfenton/seqseeker/actions/workflows/ci.yml/badge.svg)](https://github.com/ethanfenton/seqseeker/actions)
[![PyPI](https://img.shields.io/pypi/v/seqseeker.svg)](https://pypi.org/project/seqseeker)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

---

## Table of contents

- [Installation](#installation)
- [Quick start](#quick-start)
- [CLI reference](#cli-reference)
- [Python API](#python-api)
- [R interface](#r-interface)
- [Bash wrapper](#bash-wrapper)
- [Built-in sequences](#built-in-sequences)
- [Output files](#output-files)
- [Loading into Seurat](#loading-into-seurat)
- [Why not just build a custom genome?](#why-not-just-build-a-custom-genome)
- [How it works](#how-it-works)
- [Tuning guide](#tuning-guide)
- [Citation](#citation)

---

## Installation

**Requirements:** Python ≥ 3.10, gcc or cc (for compiling the C core), zlib.

```bash
pip install seqseeker
```

On first use `seqseeker` automatically compiles its C backend. You can also
pre-compile manually:

```bash
seqseeker compile
```

### macOS
```bash
# Install Xcode Command Line Tools if you don't have gcc:
xcode-select --install
pip install seqseeker
```

### Linux
```bash
sudo apt-get install build-essential zlib1g-dev   # Debian/Ubuntu
# or
sudo yum install gcc zlib-devel                    # RHEL/CentOS
pip install seqseeker
```

### R package

```r
# Requires the Python package to be installed first (above).
# Install from GitHub:
if (!requireNamespace("remotes")) install.packages("remotes")
remotes::install_github("ethanfenton/seqseeker", subdir = "R/seqseeker")
```

---

## Quick start

### Python / CLI

```bash
# Search for eGFP using a built-in sequence
seqseeker search \
  --sequence-name egfp \
  --input /data/novogene/01.RawData \
  --output results/gfp

# Any gene — pass a FASTA file
seqseeker search \
  --sequence-file my_transgene.fasta \
  --input /data/fastq/ \
  --output results/transgene

# Raw sequence string
seqseeker search \
  --sequence ATGGTGAGCAAGGGCGAGGAG... \
  --input /data/fastq/ \
  --output results/custom

# Only the 3' end (faster, specific)
seqseeker search \
  --sequence-name egfp \
  --region 3prime --region-len 200 \
  --input /data/fastq/ \
  --output results/gfp_3p

# Parallel processing
seqseeker search \
  --sequence-name egfp \
  --input /data/fastq/ \
  --workers 0            # 0 = auto-detect CPUs
```

### Python API

```python
import seqseeker

results = seqseeker.search(
    input_dirs=["/data/novogene/01.RawData",
                "/data/novogene2/01.RawData"],
    sequence_name="egfp",
    output_dir="results/gfp",
    workers=8,
)

print(f"Total hits:     {results['n_hits']:,}")
print(f"Unique barcodes:{results['n_unique_barcodes']:,}")
```

### R

```r
library(seqseeker)

results <- seqseeker(
  input         = c("/data/novogene/01.RawData",
                    "/data/novogene2/01.RawData"),
  sequence_name = "egfp",
  output        = "results/gfp",
  workers       = 8L
)

# Add to Seurat
seurat_obj <- seqseeker_add_to_seurat(seurat_obj, results)
FeaturePlot(seurat_obj, features = "seq_hits")
```

---

## CLI reference

```
seqseeker <command> [options]

Commands:
  search            Search FASTQ files for a target sequence
  compile           Pre-compile the C search binary
  list-sequences    Show built-in sequences

seqseeker search options:

  Input/output:
    --input/-i DIR [DIR ...]    FASTQ directories (required)
    --output/-o DIR             Output directory [default: seqseeker_output]

  Target sequence (exactly one required):
    --sequence/-q ACGT...       Raw nucleotide string
    --sequence-file/-f FILE     FASTA file
    --sequence-name/-n NAME     Built-in name (egfp, mcherry, venus, tdtomato)

  Region selection:
    --region/-R full|5prime|3prime   Part of target to index [default: full]
    --region-len/-L N               Length for 5'/3' slice [default: 200]

  k-mer tuning:
    --kmer-size/-k K            k-mer length [default: 20, range: 8–31]
    --min-kmers/-m N            Min hits to call read positive [default: 1]

  10x read layout:
    --barcode-len/-b N          Barcode length in R1 [default: 16]
    --umi-len/-u N              UMI length in R1 [default: 12]

  Parallelism:
    --workers/-w N              Parallel processes [default: 1, 0=auto]

  Filtering:
    --sample-filter SAMPLE ...  Process only these sample IDs
    --report-interval N         Progress every N reads [default: 1000000]

  Advanced:
    --binary PATH               Path to pre-compiled kmer_search binary
```

---

## Python API

```python
seqseeker.search(
    input_dirs,            # str | list[str | Path]
    output_dir = "seqseeker_output",
    *,
    sequence      = None,  # raw nucleotide string
    sequence_file = None,  # FASTA file path
    sequence_name = None,  # built-in name
    region        = "full",     # "full" | "5prime" | "3prime"
    region_len    = 200,
    kmer_size     = 20,
    barcode_len   = 16,
    umi_len       = 12,
    min_kmers     = 1,
    workers       = 1,          # 0 = auto
    sample_filter = None,
    report_interval = 1_000_000,
    binary_path   = None,
) -> dict
```

Returns a dict with keys:
`output_dir`, `raw_csv`, `per_cell_csv`, `positions_csv`,
`n_hits`, `n_unique_barcodes`.

---

## R interface

```r
# Search
results <- seqseeker(input, output, sequence_name, ...)

# Load existing results
results <- seqseeker_load_results("results/gfp")

# Add to Seurat
seurat_obj <- seqseeker_add_to_seurat(seurat_obj, results,
                                      barcode_suffix = "-1",
                                      meta_col       = "gfp_reads")

# List built-ins
seqseeker_list_sequences()
```

---

## Bash wrapper

The `bash/seqseeker` script is a thin wrapper that delegates to the Python
package. After `pip install seqseeker` the `seqseeker` entry point is already
on your PATH. The bash script is useful in cluster environments where you want
an explicit shell command without activating a Python environment:

```bash
# Set a custom Python path
SEQSEEKER_PYTHON=/path/to/env/bin/python3 seqseeker search ...
```

---

## Built-in sequences

```
Name           Length  Description
--------------------------------------------------------------
egfp              720 bp  Enhanced GFP (canonical, non-codon-optimised)
venus             720 bp  Venus YFP (T203Y/S65G/V68L/S72A mutations)
mcherry           711 bp  mCherry red fluorescent protein
tdtomato         1386 bp  tdTomato (tandem dimer Tomato)
```

Run `seqseeker list-sequences` to see the current list.

For a **codon-optimised** variant of eGFP (common in AAV constructs), use
`--sequence-file` pointing to your construct FASTA rather than the built-in.

---

## Output files

All outputs land in `--output` (default: `seqseeker_output/`):

| File | Description |
|---|---|
| `hits_raw.csv` | One row per positive read; sample, lane, barcode, UMI, hit count, position, strand, read name, R2 sequence |
| `hits_per_cell.csv` | One row per unique barcode; n_reads, n_umi, samples seen |
| `kmer_positions.csv` | Histogram of first k-mer match position along the target (useful for diagnosing partial integrations) |
| `_per_pair/` | Intermediate per-FASTQ-pair CSVs (safe to delete after aggregation) |

---

## Loading into Seurat

```r
# Method 1: via seqseeker R package
library(seqseeker)
results    <- seqseeker_load_results("results/gfp")
seurat_obj <- seqseeker_add_to_seurat(seurat_obj, results)

# Method 2: base R
hits         <- read.csv("results/gfp/hits_per_cell.csv")
hits$bc      <- paste0(hits$barcode, "-1")
seurat_obj$gfp_reads <- hits$n_reads[match(colnames(seurat_obj), hits$bc)]
seurat_obj$gfp_reads[is.na(seurat_obj$gfp_reads)] <- 0
```

---

## Why not just build a custom genome?

See [docs/why_kmer_search.md](docs/why_kmer_search.md) for the full argument.
In short:

1. **~100× faster** — no STAR index build, no full alignment pass.
2. **Higher sensitivity** — searches *all* raw reads, not just those that
   survive MAPQ filtering; catches partial integrations and junction reads.
3. **No false negatives from multi-mapping** — transgene k-mers that share
   similarity with endogenous loci are still detected if any portion is unique.
4. **No genome required** — works from the sequence string alone.
5. **Direct barcode attribution** — output is immediately Seurat-ready.

---

## How it works

1. The target sequence (and its reverse complement) is encoded into a
   compact open-addressing hash table of 2-bit k-mers at startup (~1 ms).
2. R2 reads are scanned with a rolling window; any k-mer hit increments a
   counter. Once `--min-kmers` threshold is met the read is emitted.
3. The cell barcode and UMI are extracted from the corresponding R1 read
   (first `barcode_len` + next `umi_len` bases).
4. Results are aggregated per barcode across all FASTQ pairs.

The C core processes ~50 M reads/min/core on a modern server; total time
for a 200 M read 10x library is typically **<5 min** on 4 cores.

---

## Tuning guide

| Goal | Recommended settings |
|---|---|
| Initial survey (maximum sensitivity) | k=20, min=1, region=full |
| Reduce noise after confirming signal | k=25, min=3 |
| Fast specific detection of 3' end | region=3prime, region-len=200 |
| Partial / truncated integration | k=15, min=1, region=full |
| Codon-optimised variant | --sequence-file construct.fasta |

Inspect `kmer_positions.csv` to see where hits cluster along the target — a
peak at position 600/720 (3' end of eGFP) is expected for normal expression;
a flat distribution suggests unspecific noise.

---

## Citation

If you use `seqseeker` in published work, please cite:

> Fenton E (2026). *seqseeker: fast k-mer detection of target sequences in
> single-cell FASTQ files.* GitHub: https://github.com/ethanfenton/seqseeker

---

## License

MIT © 2026 Ethan Fenton
