# Why k-mer search outperforms custom-genome alignment

## The alignment approach and its hidden costs

A common strategy for detecting exogenous sequences (GFP, mCherry, Cre, etc.)
in scRNA-seq data is to append the transgene CDS to the reference genome and
re-run Cell Ranger. This works, but it introduces several failure modes that
`seqseeker` avoids entirely.

---

## 1. Speed

| Task | Time (10x snRNA-seq, ~200M reads) |
|---|---|
| Build custom STAR index (once) | 45–90 min, 32 GB RAM |
| Cell Ranger count (per sample) | 6–24 h, 64 GB RAM |
| **seqseeker search** (per sample) | **5–15 min, <1 GB RAM** |

The k-mer scan runs at ~50 M reads / min / core on commodity hardware and
requires **no pre-built index** beyond what is generated in <1 second from
the target sequence itself.

---

## 2. No reads are thrown away before you search

Cell Ranger aligns R2 reads to the reference, then discards anything that
doesn't map above a minimum MAPQ. Reads that:

- span the transgene–genome junction,
- come from a partial/truncated integration,
- contain sequencing errors near the only mapping region, or
- map equally well to the endogenous locus (multi-mappers suppressed by MAPQ)

…are all silently dropped. `seqseeker` searches **every R2 read in the raw
FASTQ**, regardless of mappability.

---

## 3. Partial integration detection

If the viral or transposon insertion is truncated (common), the transgene may
produce reads covering only a portion of the CDS. Custom genome alignment
still requires a seed of ~50 bp to anchor; k-mer scanning at k=20 detects
signal even in reads with just one k-mer hit overlapping the target.

Use `--region 3prime --region-len 200` to index only the last 200 bp of the
CDS — the region most likely to be intact and least likely to be present in
any endogenous transcript — maximising specificity for partial integrations.

---

## 4. Strand-independent detection

Both the forward and reverse-complement of the target sequence are indexed.
This matters when:

- the construct drives antisense transcription,
- the sequencing library is not strand-specific (many multiome preps), or
- reads span the polyA signal on the reverse strand.

---

## 5. No genome assembly required

Some target sequences (codon-optimised variants, synthetic constructs) may not
exist in any public assembly. With `seqseeker` you paste in the sequence
directly (`--sequence`) or point to a FASTA file (`--sequence-file`). There is
no "reference build" step.

---

## 6. Transparent barcode attribution

Standard Cell Ranger output gives you a BAM aligned to the transgene, but
attributing those reads back to individual cell barcodes requires additional
parsing (samtools, pysam). `seqseeker` reads the barcode and UMI directly from
R1 — the same way Cell Ranger does — and outputs a per-cell table ready for
`Seurat::AddMetaData` in one step.

---

## When custom-genome alignment is still preferable

- You need **precise read counts for isoform/exon-level analysis** of the
  transgene (requires a proper GTF annotation).
- You are doing **multi-omic alignment** (ATAC + RNA in the same Cell Ranger
  run) and want unified output.
- The transgene is **highly expressed** (hundreds of reads per cell) and you
  need accurate UMI deduplication across cells.

For the use case of **binary or low-count classification** (is this cell
GFP-positive?) `seqseeker` is faster, more sensitive to partial signal, and
requires no genome rebuilding.

---

## Choosing parameters

| Parameter | Conservative | Sensitive |
|---|---|---|
| `--kmer-size` | 25 | 15 |
| `--min-kmers` | 3 | 1 |
| `--region` | `3prime` | `full` |

Start with the defaults (k=20, min=1, full) for an initial survey, then tighten
to reduce noise once you have a sense of the signal strength.
