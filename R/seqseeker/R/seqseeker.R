#' Search paired FASTQ files for a target sequence
#'
#' Calls the \code{seqseeker} Python CLI (which must be installed and on PATH,
#' or supplied via \code{python} / \code{cli_path}). Results are written to
#' \code{output_dir} and optionally returned as a list of data frames.
#'
#' @param input Character vector of directories containing \code{*.fastq.gz}
#'   pairs.
#' @param output Output directory (created if needed). Default:
#'   \code{"seqseeker_output"}.
#' @param sequence Raw nucleotide string to search for. Exactly one of
#'   \code{sequence}, \code{sequence_file}, or \code{sequence_name} must be
#'   supplied.
#' @param sequence_file Path to a FASTA file with the target sequence.
#' @param sequence_name Name of a built-in sequence, e.g. \code{"egfp"},
#'   \code{"mcherry"}. Run \code{\link{seqseeker_list_sequences}} for options.
#' @param region Which part of the target to index: \code{"full"} (default),
#'   \code{"5prime"}, or \code{"3prime"}.
#' @param region_len Length in bp for the 5'/3' slice (default 200).
#' @param kmer_size k-mer length (default 20).
#' @param barcode_len Cell barcode length in R1 (default 16 for 10x v3).
#' @param umi_len UMI length in R1 (default 12 for 10x v3).
#' @param min_kmers Minimum k-mer hits to call a read positive (default 1).
#' @param workers Number of parallel processes (default 1; 0 = auto-detect).
#' @param sample_filter Character vector of sample IDs to process; \code{NULL}
#'   processes all.
#' @param report_interval Print progress every N reads (default 1e6).
#' @param return_results If \code{TRUE} (default), load and return the output
#'   CSVs as a named list.
#' @param python Path to Python executable (default \code{"python3"}).
#'
#' @return If \code{return_results = TRUE}, a named list with elements
#'   \code{per_cell} (data frame) and \code{raw} (data frame). Otherwise
#'   \code{invisible(output)}.
#'
#' @examples
#' \dontrun{
#' # Search for eGFP in a 10x FASTQ directory
#' results <- seqseeker(
#'   input         = "/data/novogene/01.RawData",
#'   sequence_name = "egfp",
#'   output        = "gfp_results",
#'   workers       = 4
#' )
#'
#' # 3' end only (faster, still highly specific)
#' results <- seqseeker(
#'   input         = "/data/novogene/01.RawData",
#'   sequence_name = "egfp",
#'   region        = "3prime",
#'   region_len    = 200
#' )
#' }
#'
#' @export
seqseeker <- function(
    input,
    output          = "seqseeker_output",
    sequence        = NULL,
    sequence_file   = NULL,
    sequence_name   = NULL,
    region          = "full",
    region_len      = 200L,
    kmer_size       = 20L,
    barcode_len     = 16L,
    umi_len         = 12L,
    min_kmers       = 1L,
    workers         = 1L,
    sample_filter   = NULL,
    report_interval = 1000000L,
    return_results  = TRUE,
    python          = "python3"
) {
  # Validate sequence source
  n_sources <- sum(!is.null(sequence), !is.null(sequence_file),
                   !is.null(sequence_name))
  if (n_sources != 1L) {
    stop("Exactly one of `sequence`, `sequence_file`, or `sequence_name` must be supplied.")
  }
  if (!region %in% c("full", "5prime", "3prime")) {
    stop("`region` must be 'full', '5prime', or '3prime'.")
  }

  cmd <- c(
    python, "-m", "seqseeker", "search",
    "--input",           input,
    "--output",          output,
    "--region",          region,
    "--region-len",      as.integer(region_len),
    "--kmer-size",       as.integer(kmer_size),
    "--barcode-len",     as.integer(barcode_len),
    "--umi-len",         as.integer(umi_len),
    "--min-kmers",       as.integer(min_kmers),
    "--workers",         as.integer(workers),
    "--report-interval", as.integer(report_interval)
  )

  if (!is.null(sequence))      cmd <- c(cmd, "--sequence",      sequence)
  if (!is.null(sequence_file)) cmd <- c(cmd, "--sequence-file", sequence_file)
  if (!is.null(sequence_name)) cmd <- c(cmd, "--sequence-name", sequence_name)
  if (!is.null(sample_filter)) cmd <- c(cmd, "--sample-filter", sample_filter)

  message("Running seqseeker ...")
  status <- system(paste(shQuote(cmd), collapse = " "))
  if (status != 0L) stop("seqseeker exited with status ", status)

  if (return_results) {
    return(seqseeker_load_results(output))
  }
  invisible(output)
}


#' Load seqseeker output CSVs
#'
#' @param output_dir Path to the seqseeker output directory.
#' @return Named list: \code{per_cell}, \code{raw}, \code{positions}.
#' @export
seqseeker_load_results <- function(output_dir = "seqseeker_output") {
  per_cell_path <- file.path(output_dir, "hits_per_cell.csv")
  raw_path      <- file.path(output_dir, "hits_raw.csv")
  pos_path      <- file.path(output_dir, "kmer_positions.csv")

  out <- list()
  if (file.exists(per_cell_path)) out$per_cell  <- utils::read.csv(per_cell_path)
  if (file.exists(raw_path))      out$raw        <- utils::read.csv(raw_path)
  if (file.exists(pos_path))      out$positions  <- utils::read.csv(pos_path)
  out
}


#' Add seqseeker hit counts to a Seurat object
#'
#' Matches barcodes from \code{hits_per_cell.csv} to cell names in a Seurat
#' object and writes hit counts into the object metadata.
#'
#' @param seurat_obj A Seurat object.
#' @param results Output from \code{\link{seqseeker}} or
#'   \code{\link{seqseeker_load_results}}.
#' @param barcode_suffix Suffix appended to bare barcodes to match Seurat
#'   column names (default \code{"-1"}).
#' @param meta_col Name of the new metadata column (default
#'   \code{"seq_hits"}).
#'
#' @return The Seurat object with the new metadata column added.
#'
#' @examples
#' \dontrun{
#' results    <- seqseeker_load_results("gfp_results")
#' seurat_obj <- seqseeker_add_to_seurat(seurat_obj, results)
#' FeaturePlot(seurat_obj, features = "seq_hits")
#' }
#'
#' @export
seqseeker_add_to_seurat <- function(
    seurat_obj,
    results,
    barcode_suffix = "-1",
    meta_col       = "seq_hits"
) {
  if (is.null(results$per_cell)) {
    stop("`results` must contain a `per_cell` data frame.")
  }
  hits     <- results$per_cell
  hits$bc  <- paste0(hits$barcode, barcode_suffix)
  counts   <- hits$n_reads[match(colnames(seurat_obj), hits$bc)]
  counts[is.na(counts)] <- 0L
  seurat_obj[[meta_col]] <- counts
  seurat_obj
}


#' List built-in target sequences
#'
#' @param python Path to Python executable.
#' @return Invisibly, the character output from the CLI.
#' @export
seqseeker_list_sequences <- function(python = "python3") {
  out <- system(paste(shQuote(python), "-m seqseeker list-sequences"),
                intern = TRUE)
  cat(out, sep = "\n")
  invisible(out)
}
