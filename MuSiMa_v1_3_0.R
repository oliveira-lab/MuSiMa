#!/usr/bin/env Rscript

# --- Install and Load Required Packages ---
required_packages <- c("seqinr", "circlize", "dplyr", "Biostrings", "ComplexHeatmap", "parallel", "optparse")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(paste("Installing package:", pkg, "\n"))
    if (pkg %in% c("Biostrings", "ComplexHeatmap")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg, update = FALSE)
    } else {
      install.packages(pkg, dependencies = TRUE)
    }
  }
}

# Load packages silently
suppressPackageStartupMessages({
  library(seqinr)
  library(circlize)
  library(dplyr)
  library(Biostrings)
  library(ComplexHeatmap)
  library(parallel)
  library(optparse)
})

# --- Define Command-Line Options ---
option_list <- list(
  make_option(c("-f", "--fasta"), type = "character", help = "Comma-separated list of FASTA files"),
  make_option(c("-m", "--motif"), type = "character", help = "Motif string or file with motifs (one per line)"),
  make_option(c("-w", "--windows"), type = "character", help = "Comma-separated list of window sizes"),
  make_option(c("-s", "--step"), type = "integer", help = "Step size for sliding windows"),
  make_option(c("-t", "--method"), type = "character", help = "Method: 'uniform' or 'markov'"),
  make_option(c("-o", "--order"), type = "character", help = "Order for Markov method, or 'NA' for uniform"),
  make_option(c("-c", "--cores"), type = "integer", default = NULL, help = "Number of cores [default: all - 1]")
)

# Parse arguments
parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

# --- Extract and Validate Arguments ---
fasta_files <- strsplit(args$fasta, ",")[[1]]
motif_input <- args$motif
window_sizes <- as.numeric(strsplit(args$windows, ",")[[1]])
step_size <- args$step
method <- args$method
order_str <- args$order
cores <- args$cores

# Input validation
if (any(is.na(window_sizes)) || any(window_sizes <= 0)) stop("Error: Window sizes must be positive integers.")
if (is.na(step_size) || step_size <= 0) stop("Error: Step size must be a positive integer.")
if (!method %in% c("uniform", "markov")) stop("Error: Method must be 'uniform' or 'markov'.")
if (method == "uniform" && order_str != "NA") stop("Error: For uniform method, order must be 'NA'.")
if (method == "markov" && (is.na(as.integer(order_str)) || as.integer(order_str) < 0)) stop("Error: For Markov method, order must be a non-negative integer.")
order <- if (method == "uniform") NA else as.integer(order_str)

# Set cores for parallel processing
num_cores <- if (is.null(cores)) detectCores() - 1 else cores
if (num_cores < 1) stop("Error: Number of cores must be at least 1.")

# Load motifs (from string or file)
motifs <- if (file.exists(motif_input)) readLines(motif_input) else strsplit(motif_input, ",")[[1]]

# --- Function Definitions ---

# Load FASTA sequence
load_fasta_sequence <- function(fasta_file) {
  seq <- read.fasta(fasta_file, as.string = TRUE, forceDNAtolower = FALSE)[[1]]
  DNAString(toupper(as.character(seq)))
}

# Find motif positions (supports IUPAC degenerate codes)
find_motif_positions <- function(dna_seq, motif) {
  matches <- matchPattern(motif, dna_seq, fixed = FALSE)
  if (length(matches) == 0) return(data.frame(start = integer(), end = integer()))
  data.frame(start = start(matches), end = end(matches))
}

# Precompute motif counts per step
precompute_motif_counts <- function(motif_positions, seq_length, step_size) {
  bins <- floor((seq_length - 1) / step_size) + 1
  counts <- numeric(bins)
  for (i in seq_along(motif_positions$start)) {
    bin <- floor((motif_positions$start[i] - 1) / step_size) + 1
    if (bin <= bins) counts[bin] <- counts[bin] + 1
  }
  counts
}

# Compute transition probabilities for Markov model up to order k
get_all_transitions <- function(dna_seq, k) {
  transitions <- list()
  for (m in 0:k) {
    if (m == 0) {
      freq <- table(strsplit(as.character(dna_seq), "")[[1]]) / length(dna_seq)
      transitions[[1]] <- freq
    } else {
      mmer_counts <- oligonucleotideFrequency(dna_seq, width = m, step = 1)
      m1mer_counts <- oligonucleotideFrequency(dna_seq, width = m + 1, step = 1)
      transitions[[m + 1]] <- list()
      for (x in names(mmer_counts)) {
        if (mmer_counts[x] > 0) {
          possible_b <- c("A", "C", "G", "T")
          m1mers_x <- paste0(x, possible_b)
          counts_xb <- m1mer_counts[m1mers_x]
          counts_xb[is.na(counts_xb)] <- 0
          probs <- counts_xb / mmer_counts[x]
          transitions[[m + 1]][[x]] <- setNames(probs, possible_b)
        }
      }
    }
  }
  transitions
}

# Calculate expected motif occurrences
calc_expected_motif <- function(dna_seq, motif, window_size, method, order) {
  motif_length <- nchar(motif)
  if (method == "uniform") {
    total_positions <- length(dna_seq) - motif_length + 1
    observed_total <- length(matchPattern(motif, dna_seq, fixed = FALSE))
    p_motif <- observed_total / total_positions
  } else { # Markov
    transitions <- get_all_transitions(dna_seq, order)
    iupac_map <- list(
      "A" = "A", "C" = "C", "G" = "G", "T" = "T", "N" = c("A", "C", "G", "T"),
      "R" = c("A", "G"), "Y" = c("C", "T"), "M" = c("A", "C"), "K" = c("G", "T"),
      "S" = c("C", "G"), "W" = c("A", "T"), "B" = c("C", "G", "T"), "D" = c("A", "G", "T"),
      "H" = c("A", "C", "T"), "V" = c("A", "C", "G")
    )
    sequences <- list(strsplit(motif, "")[[1]])
    for (i in 1:length(sequences[[1]])) {
      bases <- iupac_map[[sequences[[1]][i]]]
      if (length(bases) > 1) {
        new_seqs <- list()
        for (seq in sequences) {
          for (b in bases) {
            new_seq <- seq
            new_seq[i] <- b
            new_seqs <- c(new_seqs, list(new_seq))
          }
        }
        sequences <- new_seqs
      }
    }
    p_motif <- 0
    for (seq in sequences) {
      p_seq <- 1
      for (i in 1:length(seq)) {
        m <- min(order, i - 1)
        if (m == 0) {
          p <- transitions[[1]][seq[i]]
        } else {
          context <- paste(seq[(i - m):(i - 1)], collapse = "")
          p <- transitions[[m + 1]][[context]][seq[i]]
        }
        if (is.na(p) || p == 0) {
          p_seq <- 0
          break
        }
        p_seq <- p_seq * p
      }
      p_motif <- p_motif + p_seq
    }
  }
  p_motif * (window_size - motif_length + 1)
}

# Process window size
process_window_size <- function(window_size, dna_sequences, motif, step_size, method, order, motif_counts_list, fasta_names) {
  cat(sprintf("Processing window size: %d for motif: %s\n", window_size, motif))
  results <- data.frame()
  offset <- 0
  for (chr_index in seq_along(dna_sequences)) {
    dna_seq <- dna_sequences[[chr_index]]
    seq_length <- length(dna_seq)
    motif_counts <- motif_counts_list[[chr_index]]
    num_windows <- floor((seq_length - window_size) / step_size) + 1
    expected <- calc_expected_motif(dna_seq, motif, window_size, method, order)
    steps_per_window <- ceiling(window_size / step_size)
    for (i in 1:num_windows) {
      window_start <- (i - 1) * step_size + 1 + offset
      window_end <- window_start + window_size - 1
      if (window_end > seq_length + offset) next
      start_bin <- (i - 1) + 1
      end_bin <- min(start_bin + steps_per_window - 1, length(motif_counts))
      observed <- sum(motif_counts[start_bin:end_bin])
      n <- window_size - nchar(motif) + 1
      p_motif <- expected / n
      sd <- sqrt(n * p_motif * (1 - p_motif))
      z_score <- if (sd == 0) 0 else (observed - expected) / sd
      results <- rbind(results, data.frame(
        chrom = fasta_names[chr_index], start = window_start, end = window_end,
        observed = observed, expected = expected, sd = sd, z_score = z_score
      ))
    }
    offset <- offset + seq_length
  }
  results
}

# --- Main Processing Loop ---
summary_df <- data.frame()
for (motif in motifs) {
  motif <- toupper(motif)
  L <- nchar(motif)
  if (method == "markov" && (order > L - 2)) stop(sprintf("Error: For motif %s (length %d), order must be <= %d", motif, L, L - 2))
  
  # Load sequences
  dna_sequences <- lapply(fasta_files, load_fasta_sequence)
  motif_positions_list <- lapply(dna_sequences, find_motif_positions, motif)
  motif_counts_list <- mapply(precompute_motif_counts, motif_positions_list, sapply(dna_sequences, length), MoreArgs = list(step_size), SIMPLIFY = FALSE)
  
  # Extract FASTA file base names without extension
  fasta_names <- tools::file_path_sans_ext(basename(fasta_files))
  
  # Chromosome data for plotting using FASTA file names
  chr_lengths <- sapply(dna_sequences, length)
  df <- data.frame(
    name = fasta_names,
    start = cumsum(c(1, head(chr_lengths, -1))),
    end = cumsum(chr_lengths)
  )
  
  # Generate plot
  pdf_file <- paste0("musima_plot_", gsub("[^A-Za-z0-9]", "_", motif), ".pdf")
  pdf(pdf_file, width = 15, height = 10)
  circos.par("track.height" = 0.05, "gap.degree" = 5, "start.degree" = 90)
  circos.genomicInitialize(df)
  col_fun <- colorRamp2(c(-10, -5, 0, 5, 10), c("blue", "navy", "green", "red", "darkred"))
  
  # Process windows in parallel
  results_list <- mclapply(window_sizes, process_window_size, dna_sequences, motif, step_size, method, order, motif_counts_list, fasta_names, mc.cores = num_cores)
  
  # Plot tracks
  for (results in results_list) {
    plot_data <- results[, c("chrom", "start", "end", "z_score")]
    circos.genomicTrack(plot_data, stack = TRUE, panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = col_fun(value[[1]]), border = NA, ...)
    })
  }
  draw(Legend(col_fun = col_fun, title = "Motif Z-Score", direction = "vertical"), x = unit(0.9, "npc"), y = unit(0.5, "npc"))
  circos.clear()
  dev.off()
  
  # Save results
  output_file <- paste0("results_list_", gsub("[^A-Za-z0-9]", "_", motif), ".txt")
  file_conn <- file(output_file, "w")
  for (i in seq_along(window_sizes)) {
    cat(sprintf("\n=== Window Size: %d ===\n", window_sizes[i]), file = file_conn)
    write.table(results_list[[i]], file = file_conn, sep = "\t", quote = FALSE, row.names = FALSE)
    cat("\n", file = file_conn)
  }
  close(file_conn)
  
  # Summary statistics per FASTA file
  for (i in seq_along(window_sizes)) {
    results <- results_list[[i]]
    summary_per_file <- results %>%
      group_by(chrom) %>%
      summarize(mean_z = mean(z_score),
                num_sig = sum(abs(z_score) > 2))
    summary_df <- rbind(summary_df, data.frame(
      motif = motif,
      window_size = window_sizes[i],
      fasta_file = summary_per_file$chrom,
      mean_z = summary_per_file$mean_z,
      num_sig = summary_per_file$num_sig
    ))
  }
}

# Save summary statistics
write.table(summary_df, "summary_stats.txt", sep = "\t", quote = FALSE, row.names = FALSE)
cat("Analysis complete. Plots and results saved for each motif.\n")
