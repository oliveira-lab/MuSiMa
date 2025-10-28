#!/usr/bin/env Rscript

# ---- Install and Load Required Packages ----
required_packages <- c("seqinr", "circlize", "dplyr", "Biostrings", "ComplexHeatmap", "parallel", "optparse", "logger")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    tryCatch({
      if (pkg %in% c("Biostrings", "ComplexHeatmap")) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager", dependencies = TRUE, repos = "https://cloud.r-project.org")
        }
        BiocManager::install(pkg, update = FALSE)
      } else {
        install.packages(pkg, dependencies = TRUE, repos = "https://cloud.r-project.org")
      }
      if (!requireNamespace(pkg, quietly = TRUE)) stop(sprintf("Failed to install package %s", pkg))
    }, error = function(e) stop(sprintf("Error installing %s: %s", pkg, e$message)))
  }
}

# Load packages
suppressPackageStartupMessages({
  suppressWarnings({
    library(seqinr)
    library(circlize)
    library(dplyr)
    library(Biostrings)
    library(ComplexHeatmap)
    library(parallel)
    library(optparse)
    library(logger)
  })
})

# ---- Define Command-Line Options ----
option_list <- list(
  make_option(c("-f", "--fasta"), type = "character", 
              help = "For motif mode: Comma-separated list of FASTA files (each can be a multifasta file). For Signal mode: Single FASTA file containing chromosome sequences."),
  make_option(c("-m", "--motif"), type = "character", 
              help = "Motif string or file with motifs (one per line)"),
  make_option(c("-w", "--windows"), type = "character", 
              help = "Comma-separated list of window sizes"),
  make_option(c("-s", "--step"), type = "integer", 
              help = "Step size for sliding windows"),
  make_option(c("-t", "--method"), type = "character", 
              help = "Method: 'uniform' or 'markov'"),
  make_option(c("-o", "--order"), type = "character", 
              help = "Order for Markov method, or 'NA' for uniform"),
  make_option(c("-c", "--cores"), type = "integer", default = NULL, 
              help = "Number of cores [default: all - 1]"),
  make_option(c("-b", "--bed"), type = "character", default = NULL,
              help = "BED file (4 columns: ID, start, end, strand). If provided, script runs BED mode."),
  make_option(c("-r", "--randomization"), type = "integer", default = NULL,
              help = "Number of randomization iterations for Signal mode (recommended >100)."),
  make_option(c("--plotrange"), type = "character", action = "store", default = "-5,5", 
              help = "Range for plotting z-scores in both modes (e.g. --plotrange=-5,5)")
)

# Parse arguments
parser <- OptionParser(option_list = option_list, description = "MuSiMa: Motif and region enrichment analysis.\nMotif mode: Scan for motif occurrences in sliding windows.\nSignal mode: Analyze region enrichment with randomization. In Signal mode, p-values are empirical two-tailed.")
args <- parse_args(parser)

# --- Extract and Validate Arguments ----
fasta_files <- if (!is.null(args$fasta)) strsplit(args$fasta, ",")[[1]] else NULL
motif_input <- args$motif
window_sizes <- if (!is.null(args$windows)) as.numeric(strsplit(args$windows, ",")[[1]]) else NULL
step_size <- args$step
method <- args$method
order_str <- args$order
cores <- args$cores
bed_file <- args$bed
randomization_n <- args$randomization
zrange <- as.numeric(strsplit(args$plotrange, ",")[[1]])

# Determine mode
is_signal_mode <- !is.null(bed_file)

# Initialize logging after mode determination
ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_file <- if (is_signal_mode) paste0("musima_signal_", ts, ".log") else paste0("musima_motif_", ts, ".log")
log_appender(appender_file(log_file))

# Input validation for non-Signal mode (motif mode)
if (!is_signal_mode) {
  if (is.null(fasta_files) || is.null(motif_input) || is.null(window_sizes) || is.null(step_size) || is.null(method) || is.null(order_str)) {
    stop("Error: For motif mode you must provide --fasta, --motif, --windows, --step, --method and --order (or provide --bed to run Signal mode).")
  }
  if (!is.null(randomization_n)) {
    stop("Error: Option --randomization is not allowed in motif mode.")
  }
  if (!all(file.exists(fasta_files))) stop("Error: One or more FASTA files not found.")
  if (any(is.na(window_sizes)) || any(window_sizes <= 0)) stop("Error: Window sizes must be positive integers.")
  if (is.na(step_size) || step_size <= 0) stop("Error: Step size must be a positive integer.")
  if (!method %in% c("uniform", "markov")) stop("Error: Method must be 'uniform' or 'markov'.")
  if (method == "uniform" && order_str != "NA") stop("Error: For uniform method, order must be 'NA'.")
  if (method == "markov" && (is.na(as.integer(order_str)) || as.integer(order_str) < 0)) 
    stop("Error: For Markov method, order must be a non-negative integer.")
  if (any(step_size / window_sizes < 0.1)) 
    warning("Step size is small relative to window size, which may lead to correlated results.")
  order <- if (method == "uniform") NA else as.integer(order_str)
}

# Signal mode validation
if (is_signal_mode) {
  if (is.null(fasta_files)) stop("Error: --fasta (single genome FASTA) is required for Signal mode.")
  if (length(fasta_files) != 1) stop("Error: For Signal mode, --fasta must specify exactly one FASTA file.")
  if (!file.exists(bed_file)) stop("Error: Signal file not found.")
  if (!file.exists(fasta_files)) stop("Error: Genome FASTA file not found.")
  if (is.null(step_size) || step_size <= 0) stop("Error: Step size must be a positive integer.")
  if (is.null(window_sizes) || any(is.na(window_sizes)) || any(window_sizes <= 0)) stop("Error: --windows must be provided and contain positive integer(s).")
  if (is.null(randomization_n)) stop("Error: --randomization must be provided and must be a positive integer in Signal mode.")
  if (randomization_n < 1) stop("Error: --randomization must be a positive integer.")
  if (!is.null(motif_input) || !is.null(method) || !is.null(order_str)) {
    stop("Error: Options --motif, --method, and --order are not allowed in Signal mode.")
  }
}

# Set cores
num_cores <- if (Sys.info()["sysname"] == "Windows") 1 else min(detectCores(), cores %||% detectCores() - 1)
if (num_cores < 1) stop("Error: Number of cores must be at least 1.")

# Timestamped filename function
timestamped_filename <- function(base, ext = ".pdf") {
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  paste0(base, "_", ts, ext)
}

# Load motifs (only if motif mode)
if (!is_signal_mode) {
  motifs <- if (file.exists(motif_input)) readLines(motif_input) else strsplit(motif_input, ",")[[1]]
  motifs <- trimws(motifs)
  valid_iupac <- c("A","C","G","T","N","R","Y","M","K","S","W","B","D","H","V")
  for (m in motifs) {
    if (!all(strsplit(m, "")[[1]] %in% valid_iupac)) 
      stop(sprintf("Error: Motif '%s' contains invalid IUPAC codes.", m))
  }
}

# ---- Function Definitions ----

# Load sequences from FASTA files
load_all_sequences <- function(fasta_files) {
  all_sequences <- lapply(fasta_files, read.fasta, as.string = TRUE, forceDNAtolower = FALSE)
  combined_sequences <- do.call(c, all_sequences)
  dna_sequences <- lapply(combined_sequences, function(seq) DNAString(toupper(as.vector(seq))))
  sequence_names <- names(combined_sequences)
  list(dna_sequences = dna_sequences, sequence_names = sequence_names)
}

# Load single genome fasta (for Signal mode)
load_genome_fasta <- function(fasta_file) {
  seqs <- read.fasta(fasta_file, as.string = TRUE, forceDNAtolower = FALSE)
  dna_sequences <- lapply(seqs, function(seq) DNAString(toupper(as.vector(seq))))
  sequence_names <- names(seqs)
  list(dna_sequences = dna_sequences, sequence_names = sequence_names)
}

# Find motif positions (supports IUPAC degenerate codes)
find_motif_positions <- function(dna_seq, motif) {
  has_ambiguities <- any(strsplit(motif, "")[[1]] %in% c("N","R","Y","M","K","S","W","B","D","H","V"))
  fixed_val <- !has_ambiguities  # TRUE for exact motifs, FALSE for ambiguous
  matches <- matchPattern(motif, dna_seq, fixed = fixed_val)
  if (length(matches) == 0) return(data.frame(start = integer(), end = integer()))
  data.frame(start = start(matches), end = end(matches))
}

# Compute transition probabilities for Markov model up to order k
get_all_transitions <- function(dna_seq, k) {
  transitions <- list()
  for (m in 0:k) {
    if (m == 0) {
      freq <- (table(strsplit(as.character(dna_seq), "")[[1]]) + 1) / (length(dna_seq) + 4)
      transitions[[1]] <- freq
    } else {
      m1mer_counts <- oligonucleotideFrequency(dna_seq, width = m + 1, step = 1) + 1
      transitions[[m + 1]] <- list()
      possible_b <- c("A", "C", "G", "T")
      all_mmers <- names(oligonucleotideFrequency(dna_seq, width = m, step = 1))
      for (x in all_mmers) {
        m1mers_x <- paste0(x, possible_b)
        counts_xb <- m1mer_counts[m1mers_x]
        probs <- counts_xb / sum(counts_xb)
        transitions[[m + 1]][[x]] <- setNames(probs, possible_b)
      }
    }
  }
  transitions
}

# Calculate expected motif occurrences
calc_expected_motif <- function(dna_seq, motif, window_size, method, order) {
  motif_length <- nchar(motif)
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
  if (method == "uniform") {
    base_counts <- letterFrequency(dna_seq, c("A","C","G","T")) + 1
    base_probs <- base_counts / sum(base_counts)
    p_motif <- 0
    for (seq in sequences) {
      p_seq <- 1
      for (b in seq) {
        p <- base_probs[b]
        p_seq <- p_seq * p
      }
      p_motif <- p_motif + p_seq
    }
  } else { # Markov
    transitions <- get_all_transitions(dna_seq, order)
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

# Process window size (motif mode)
process_window_size <- function(window_size, dna_sequences, motif, step_size, method, order, prefix_sums_list, sequence_names, strand) {
  log_info("Processing window size {window_size} for motif {motif} on {strand} strand")
  results <- data.frame()
  L <- nchar(motif)
  for (chr_index in seq_along(dna_sequences)) {
    dna_seq <- dna_sequences[[chr_index]]
    seq_length <- length(dna_seq)
    prefix <- prefix_sums_list[[chr_index]]
    num_windows <- floor((seq_length - window_size) / step_size) + 1
    expected <- calc_expected_motif(dna_seq, motif, window_size, method, order)
    n <- window_size - L + 1
    if (n <= 0) next
    for (i in 1:num_windows) {
      window_start <- (i - 1) * step_size + 1
      window_end <- window_start + window_size - 1
      if (window_end > seq_length) next
      start_pos <- window_start
      end_pos <- start_pos + n - 1
      observed <- prefix[end_pos + 1] - prefix[start_pos]
      p_motif <- expected / n
      sd <- sqrt(n * p_motif * (1 - p_motif))
      lower <- observed < expected
      if (expected < 5) {
        p_value <- if (lower) {
          ppois(observed, expected, lower.tail = TRUE)
        } else {
          ppois(observed - 1, expected, lower.tail = FALSE)
        }
      } else {
        cont_obs <- observed + ifelse(lower, 0.5, -0.5)
        p_value <- pnorm(cont_obs, expected, sd, lower.tail = lower)
      }
      z_score <- ifelse(sd == 0, 0, (observed - expected) / sd)
      results <- rbind(results, data.frame(
        chrom = sequence_names[chr_index], 
        start = window_start, 
        end = window_end,
        strand = strand,
        observed = observed, 
        expected = expected, 
        sd = sd, 
        z_score = z_score,
        p_value = p_value
      ))
    }
  }
  results$p_value <- p.adjust(results$p_value, method = "fdr")
  results
}

# Check if a motif is palindromic
is_palindromic <- function(motif) {
  motif <- toupper(motif)
  rc_motif <- as.character(reverseComplement(DNAString(motif)))
  return(motif == rc_motif)
}

# --- Helper functions for Signal mode ----

# Read BED with 4 columns (ID, start, end, strand). Start assumed 1-based.
read_bed4 <- function(bed_path) {
  bed <- read.table(bed_path, header = FALSE, stringsAsFactors = FALSE)
  if (ncol(bed) < 4) stop("BED file must have at least 4 columns: ID, start, end, strand.")
  bed <- bed[, 1:4]
  colnames(bed) <- c("id", "start", "end", "strand")
  bed$start <- as.integer(bed$start)
  bed$end <- as.integer(bed$end)
  if (any(is.na(bed$start)) || any(is.na(bed$end)) || any(bed$start > bed$end)) stop("Invalid start/end coordinates in BED file.")
  bed
}

# Convert per-chromosome coordinates to global coordinates using offsets
build_chr_offsets <- function(sequence_lengths, sequence_names) {
  starts <- cumsum(c(1, head(sequence_lengths, -1)))
  ends <- cumsum(sequence_lengths)
  data.frame(name = sequence_names, start = starts, end = ends, length = sequence_lengths, stringsAsFactors = FALSE)
}

# Map bed to global coordinates (chrom name matching required)
bed_to_global <- function(bed_df, chr_offsets) {
  bed_df$chrom <- bed_df$id
  merged <- merge(bed_df, chr_offsets, by.x = "chrom", by.y = "name", all.x = TRUE, sort = FALSE)
  if (any(is.na(merged$start.y))) stop("Some BED chromosomes not found in genome FASTA.")
  merged$global_start <- merged$start.x + merged$start.y - 1
  merged$global_end <- merged$end.x + merged$start.y - 1
  merged[, c("chrom", "start.x", "end.x", "strand", "global_start", "global_end")]
}

# Count overlaps of regions against windows (both in global coordinates)
count_overlaps_windows <- function(regions_global, windows_global_start, windows_global_end) {
  obs <- integer(length(windows_global_start))
  for (i in seq_along(windows_global_start)) {
    s <- windows_global_start[i]
    e <- windows_global_end[i]
    obs[i] <- sum(regions_global$start <= e & regions_global$end >= s)
  }
  obs
}

# Randomize regions: sample random global start positions for each region (preserve length)
randomize_regions_global <- function(regions, chr_offsets) {
  randomized <- data.frame()
  for (i in seq_len(nrow(chr_offsets))) {
    chr <- chr_offsets$name[i]
    chr_start <- chr_offsets$start[i]
    chr_len <- chr_offsets$length[i]
    chr_regions <- regions[regions$chrom == chr, ]
    if (nrow(chr_regions) == 0) next
    widths <- chr_regions$global_end - chr_regions$global_start + 1
    new_starts <- integer(nrow(chr_regions))
    for (j in seq_len(nrow(chr_regions))) {
      max_start <- chr_len - widths[j] + 1
      if (max_start <= 0) next
      new_starts[j] <- sample.int(max_start, 1, replace = TRUE)
    }
    if (all(new_starts == 0)) next
    new_ends <- new_starts + widths - 1
    randomized <- rbind(randomized, data.frame(
      chrom  = chr,
      start  = chr_start + new_starts - 1,
      end    = chr_start + new_ends - 1,
      strand = chr_regions$strand
    ))
  }
  if (nrow(randomized) == 0) warning("No regions generated after randomization")
  randomized
}

# Map a global coordinate to chromosome (binary search like)
global_to_chr <- function(pos, chr_offsets) {
  idx <- which(chr_offsets$start <= pos & chr_offsets$end >= pos)
  if (length(idx) == 0) return(NA)
  chr_offsets$name[idx]
}

# ---- Signal Mode Implementation --------
if (is_signal_mode) {
  log_info("Starting Signal mode analysis")
  cat("Computing. Please wait...\n")
  # Load genome FASTA (single file)
  genome <- load_genome_fasta(fasta_files[1])
  dna_sequences <- genome$dna_sequences
  sequence_names <- genome$sequence_names
  chr_lengths <- sapply(dna_sequences, length)
  chr_offsets <- build_chr_offsets(chr_lengths, sequence_names)
  total_genome_length <- sum(chr_lengths)
  
  # Read BED and convert to global coords
  bed_df <- read_bed4(bed_file)
  regions_global <- bed_to_global(bed_df, chr_offsets)
  
  # Build windows in global coordinates
  windows_info_per_window_size <- list()
  for (w in window_sizes) {
    windows <- data.frame()
    for (i in seq_along(sequence_names)) {
      chr <- sequence_names[i]
      chr_len <- chr_lengths[i]
      offset <- chr_offsets$start[i] - 1
      num_windows <- floor((chr_len - w) / step_size) + 1
      if (num_windows < 1) next
      for (j in 1:num_windows) {
        local_start <- (j - 1) * step_size + 1
        local_end <- local_start + w - 1
        global_start <- local_start + offset
        global_end <- local_end + offset
        windows <- rbind(windows, data.frame(
          chrom = chr,
          start = global_start,
          end = global_end,
          local_start = local_start,
          local_end = local_end
        ))
      }
    }
    windows_info_per_window_size[[as.character(w)]] <- windows
  }
  
  # Observed counts per strand
  obs_counts_list <- list()
  for (w in window_sizes) {
    log_info("Processing window size {w} in Signal mode")
    windows <- windows_info_per_window_size[[as.character(w)]]
    if (nrow(windows) == 0) {
      obs_counts_list[[as.character(w)]] <- list(plus = integer(0), minus = integer(0))
      next
    }
    obs_plus  <- count_overlaps_windows(regions_global[regions_global$strand == "+", ], windows$start, windows$end)
    obs_minus <- count_overlaps_windows(regions_global[regions_global$strand == "-", ], windows$start, windows$end)
    obs_counts_list[[as.character(w)]] <- list(plus = obs_plus, minus = obs_minus)
  }
  
  # Helper: empirical p-value (two-tailed)
  calc_empirical_p <- function(all_counts, observed) {
    n_rand <- nrow(all_counts)
    pvals <- numeric(ncol(all_counts))
    for (i in seq_along(observed)) {
      left_tail <- sum(all_counts[, i] <= observed[i]) / n_rand
      right_tail <- sum(all_counts[, i] >= observed[i]) / n_rand
      pvals[i] <- min(1, 2 * min(left_tail, right_tail))
    }
    pvals
  }
  
  # Signal mode randomization
  random_stats_per_window_size <- list()
  for (w in window_sizes) {
    windows <- windows_info_per_window_size[[as.character(w)]]
    n_windows <- nrow(windows)
    if (n_windows == 0) {
      random_stats_per_window_size[[as.character(w)]] <- list(
        plus = list(mean = numeric(0), sd = numeric(0), all_counts = matrix(0,0,0)),
        minus = list(mean = numeric(0), sd = numeric(0), all_counts = matrix(0,0,0))
      )
      next
    }
    counts_list <- mclapply(seq_len(randomization_n), function(iter_idx) {
      tryCatch({
        rand_regs <- randomize_regions_global(regions_global, chr_offsets)
        if (nrow(rand_regs) == 0) {
          return(list(plus = rep(0, n_windows), minus = rep(0, n_windows)))
        }
        plus_counts <- count_overlaps_windows(rand_regs[rand_regs$strand == "+", ], windows$start, windows$end)
        minus_counts <- count_overlaps_windows(rand_regs[rand_regs$strand == "-", ], windows$start, windows$end)
        list(plus = plus_counts, minus = minus_counts)
      }, error = function(e) {
        message(sprintf("Randomization iteration %d failed: %s", iter_idx, e$message))
        return(list(plus = rep(0, n_windows), minus = rep(0, n_windows)))
      })
    }, mc.cores = num_cores)
    if (length(counts_list) == 0) {
      stop(sprintf("No randomization results for window size %d", w))
    }
    all_plus <- do.call(rbind, lapply(counts_list, `[[`, "plus"))
    all_minus <- do.call(rbind, lapply(counts_list, `[[`, "minus"))
    mean_plus <- colMeans(all_plus)
    mean_minus <- colMeans(all_minus)
    sd_plus <- apply(all_plus, 2, sd)
    sd_minus <- apply(all_minus, 2, sd)
    random_stats_per_window_size[[as.character(w)]] <- list(
      plus = list(mean = mean_plus, sd = sd_plus, all_counts = all_plus),
      minus = list(mean = mean_minus, sd = sd_minus, all_counts = all_minus)
    )
  }
  
  # Assemble results with strand + p-values
  results_list_bed <- list()
  for (w in window_sizes) {
    windows <- windows_info_per_window_size[[as.character(w)]]
    obs_plus <- obs_counts_list[[as.character(w)]]$plus
    obs_minus <- obs_counts_list[[as.character(w)]]$minus
    stats_plus <- random_stats_per_window_size[[as.character(w)]]$plus
    stats_minus <- random_stats_per_window_size[[as.character(w)]]$minus
    res_plus <- data.frame(
      chrom = windows$chrom,
      start = windows$local_start,
      end = windows$local_end,
      strand = "+",
      observed = obs_plus,
      expected = stats_plus$mean,
      sd = stats_plus$sd,
      z_score = ifelse(stats_plus$sd == 0, 0, (obs_plus - stats_plus$mean) / stats_plus$sd),
      p_value = calc_empirical_p(stats_plus$all_counts, obs_plus)
    )
    res_minus <- data.frame(
      chrom = windows$chrom,
      start = windows$local_start,
      end = windows$local_end,
      strand = "-",
      observed = obs_minus,
      expected = stats_minus$mean,
      sd = stats_minus$sd,
      z_score = ifelse(stats_minus$sd == 0, 0, (obs_minus - stats_minus$mean) / stats_minus$sd),
      p_value = calc_empirical_p(stats_minus$all_counts, obs_minus)
    )
    res <- rbind(res_plus, res_minus)
    res$p_value <- p.adjust(res$p_value, method = "fdr")
    results_list_bed[[as.character(w)]] <- res
  }
  
  # Compute summary statistics for signal mode
  summary_df <- data.frame()
  for (w in window_sizes) {
    results <- results_list_bed[[as.character(w)]]
    summary_per_sequence <- results %>%
      group_by(chrom, strand) %>%
      summarize(mean_p = mean(-log10(p_value)),
                num_sig = sum(p_value < 0.05),
                .groups = "drop")
    summary_df <- rbind(summary_df, data.frame(
      window_size = w,
      sequence = summary_per_sequence$chrom,
      strand = summary_per_sequence$strand,
      mean_p = summary_per_sequence$mean_p,
      num_sig = summary_per_sequence$num_sig
    ))
  }
  summary_file <- timestamped_filename("signal_summary_stats", ".txt")
  write.table(summary_df, summary_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Plotting
  df_plot <- data.frame(
    name = sequence_names,
    start = rep(1, length(sequence_names)),
    end = chr_lengths
  )
  
  pdf_file <- timestamped_filename("signal_musima_plot", ".pdf")
  pdf(pdf_file, width = 15, height = 10)
  
  num_sectors <- nrow(df_plot)
  gap_degrees <- c(15, rep(5, num_sectors - 1))
  sum_gap_degrees <- sum(gap_degrees)
  total_angle_for_sectors <- 360 - sum_gap_degrees
  prop_first <- chr_lengths[1] / sum(chr_lengths)
  theta1 <- total_angle_for_sectors * prop_first
  start_degree <- (90 - theta1 - (gap_degrees[1] / 2)) %% 360
  
  circos.par(
    track.height = 0.045,
    gap.degree = gap_degrees,
    track.margin = c(0.005, 0.005),
    start.degree = start_degree
  )
  
  circos.genomicInitialize(df_plot)
  col_fun <- colorRamp2(seq(zrange[1], zrange[2], length.out = 5), c("blue", "navy", "green", "red", "darkred"))
  
  # Outer (+) tracks
  for (i in seq_along(window_sizes)) {
    w <- window_sizes[i]
    plot_data <- results_list_bed[[as.character(w)]]
    if (nrow(plot_data) == 0) next
    sense_data <- plot_data[plot_data$strand == "+", ]
    if (nrow(sense_data) > 0) {
      track_df <- sense_data[, c("chrom", "start", "end", "z_score")]
      track_margin <- if (i == length(window_sizes)) c(0.04, 0.005) else c(0.005, 0.005)
      circos.genomicTrack(track_df, stack = TRUE, track.margin = track_margin,
        panel.fun = function(region, value, ...) {
          circos.genomicRect(region, value, col = col_fun(value[[1]]), border = NA, ...)
          if (get.current.sector.index() == sequence_names[1]) {
            xlim <- get.cell.meta.data("xlim")
            circos.text(x = xlim[2], y = 1.05, labels = paste(w, "(+)"),
                        facing = "inside", niceFacing = TRUE, adj = c(0,0.5), cex = 0.6)
          }
        })
    }
  }
  
  # Inner (-) tracks
  for (i in seq_along(window_sizes)) {
    w <- window_sizes[i]
    plot_data <- results_list_bed[[as.character(w)]]
    if (nrow(plot_data) == 0) next
    antisense_data <- plot_data[plot_data$strand == "-", ]
    if (nrow(antisense_data) > 0) {
      track_df <- antisense_data[, c("chrom", "start", "end", "z_score")]
      track_margin <- if (i == 1) c(0.005, 0.04) else c(0.005, 0.005)
      circos.genomicTrack(track_df, stack = TRUE, track.margin = track_margin,
        panel.fun = function(region, value, ...) {
          circos.genomicRect(region, value, col = col_fun(value[[1]]), border = NA, ...)
          if (get.current.sector.index() == sequence_names[1]) {
            xlim <- get.cell.meta.data("xlim")
            circos.text(x = xlim[2], y = 1.05, labels = paste(w, "(-)"),
                        facing = "inside", niceFacing = TRUE, adj = c(0,0.5), cex = 0.6)
          }
        })
    }
  }
  
  draw(Legend(col_fun = col_fun, title = "BED region Z-Score",
              at = seq(zrange[1], zrange[2], length.out = 5),
              direction = "vertical"),
       x = unit(0.9, "npc"), y = unit(0.5, "npc"))
  circos.clear()
  dev.off()
  
  # Concatenate all results into one file with header
  all_res <- data.frame()
  for (w in window_sizes) {
    res <- results_list_bed[[as.character(w)]]
    res$window_size <- w
    all_res <- rbind(all_res, res)
  }
  outname <- timestamped_filename("signal_results_all", ".bed")
  write.table(
    all_res[, c("chrom", "start", "end", "strand", "observed", "expected", "sd", "z_score", "p_value", "window_size")],
    outname,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    append = FALSE
  )
  
  log_info("Signal mode analysis complete. Plot saved to: {pdf_file}")
  cat("Signal mode analysis complete. Plot saved to:", pdf_file, "\n")
  quit(save = "no", status = 0)
}

# ---- Main processing loop for motif mode ----
log_info("Starting motif mode analysis")
cat("Computing. Please wait...\n")
summary_df <- data.frame()

# Load sequences and their names
seq_data <- load_all_sequences(fasta_files)
dna_sequences <- seq_data$dna_sequences
sequence_names <- seq_data$sequence_names

for (motif in motifs) {
  log_info("Processing motif {motif}")
  motif <- toupper(motif)
  L <- nchar(motif)
  if (L > 0.1 * min(window_sizes)) warning("Motif length is significant relative to window size; overlaps may bias variance.")
  if (method == "markov" && (order >= L)) 
    stop(sprintf("Error: For motif %s (length %d), order must be < %d", motif, L, L))
  
  # Determine if motif is palindromic and set motifs to process
  palindromic <- is_palindromic(motif)
  if (palindromic) {
    motifs_to_process <- list(list(motif = motif, strand = "sense"))
  } else {
    rc_motif <- as.character(reverseComplement(DNAString(motif)))
    motifs_to_process <- list(list(motif = motif, strand = "sense"), list(motif = rc_motif, strand = "antisense"))
  }
  
  # Process each motif (sense and/or antisense)
  motif_results_list <- list()
  for (mt in motifs_to_process) {
    current_motif <- mt$motif
    strand <- mt$strand
    # Find motif positions and precompute prefix sums
    motif_positions_list <- lapply(dna_sequences, find_motif_positions, current_motif)
    prefix_sums_list <- mclapply(seq_along(dna_sequences), function(i) {
      starts <- motif_positions_list[[i]]$start
      seq_length <- length(dna_sequences[[i]])
      max_start <- seq_length - L + 1
      if (max_start < 1) return(c(0))
      counts <- tabulate(starts, nbins = max_start)
      c(0, cumsum(counts))
    }, mc.cores = num_cores)
    # Process windows in parallel
    results_list <- mclapply(window_sizes, process_window_size, dna_sequences, current_motif, step_size, 
                             method, order, prefix_sums_list, sequence_names, strand, mc.cores = num_cores)
    motif_results_list[[strand]] <- results_list
  }
  
  # Chromosome data for plotting using sequence names
  chr_lengths <- sapply(dna_sequences, length)
  df <- data.frame(
    name = sequence_names,
    start = rep(1, length(sequence_names)),
    end = chr_lengths
  )
  
  # Generate plot
  pdf_base <- paste0("motif_musima_plot_", gsub("[^A-Za-z0-9]", "_", motif))
  pdf_file <- timestamped_filename(pdf_base, ".pdf")
  pdf(pdf_file, width = 15, height = 10)
  
  num_sectors <- nrow(df)
  gap_degrees <- c(15, rep(5, num_sectors - 1))
  sum_gap_degrees <- sum(gap_degrees)
  total_angle_for_sectors <- 360 - sum_gap_degrees
  prop_first <- chr_lengths[1] / sum(chr_lengths)
  theta1 <- total_angle_for_sectors * prop_first
  start_degree <- (90 - theta1 - (gap_degrees[1] / 2)) %% 360
  
  circos.par(
    track.height = 0.045,
    gap.degree = gap_degrees,
    track.margin = c(0.005, 0.005),
    start.degree = start_degree
  )
  
  circos.genomicInitialize(df)
  col_fun <- colorRamp2(seq(zrange[1], zrange[2], length.out = 5), c("blue", "navy", "green", "red", "darkred"))
  
  if (palindromic) {
    results_list <- motif_results_list[["sense"]]
    for (i in seq_along(window_sizes)) {
      plot_data <- results_list[[i]][, c("chrom", "start", "end", "z_score")]
      if (nrow(plot_data) == 0) next
      window_size_label <- window_sizes[i]
      circos.genomicTrack(plot_data, stack = TRUE, track.margin = c(0.005, 0.005), panel.fun = function(region, value, ...) {
        circos.genomicRect(region, value, col = col_fun(value[[1]]), border = NA, ...)
        if (get.current.sector.index() == sequence_names[1]) {
          xlim <- get.cell.meta.data("xlim")
          circos.text(x = xlim[2], y = 1.05, labels = paste(window_size_label), 
                      facing = "inside", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.6)
        }
      })
    }
  } else {
    # Plot sense tracks with legends
    results_list_sense <- motif_results_list[["sense"]]
    for (i in seq_along(window_sizes)) {
      plot_data <- results_list_sense[[i]][, c("chrom", "start", "end", "z_score")]
      if (nrow(plot_data) == 0) next
      window_size_label <- window_sizes[i]
      track_margin <- if (i == length(window_sizes)) c(0.04, 0.005) else c(0.005, 0.005)
      circos.genomicTrack(plot_data, stack = TRUE, track.margin = track_margin, panel.fun = function(region, value, ...) {
        circos.genomicRect(region, value, col = col_fun(value[[1]]), border = NA, ...)
        if (get.current.sector.index() == sequence_names[1]) {
          xlim <- get.cell.meta.data("xlim")
          circos.text(x = xlim[2], y = 1.05, labels = paste(window_size_label, "(+)"), 
                      facing = "inside", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.6)
        }
      })
    }
    # Plot antisense tracks with legends
    results_list_antisense <- motif_results_list[["antisense"]]
    for (i in seq_along(window_sizes)) {
      plot_data <- results_list_antisense[[i]][, c("chrom", "start", "end", "z_score")]
      if (nrow(plot_data) == 0) next
      window_size_label <- window_sizes[i]
      track_margin <- if (i == 1) c(0.005, 0.04) else c(0.005, 0.005)
      circos.genomicTrack(plot_data, stack = TRUE, track.margin = track_margin, panel.fun = function(region, value, ...) {
        circos.genomicRect(region, value, col = col_fun(value[[1]]), border = NA, ...)
        if (get.current.sector.index() == sequence_names[1]) {
          xlim <- get.cell.meta.data("xlim")
          circos.text(x = xlim[2], y = 1.05, labels = paste(window_size_label, "(-)"), 
                      facing = "inside", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.6)
        }
      })
    }
  }
  
  # Add legend
  draw(Legend(col_fun = col_fun, title = "Motif Z-Score", at = seq(zrange[1], zrange[2], length.out = 5), direction = "vertical"), 
       x = unit(0.9, "npc"), y = unit(0.5, "npc"))
  circos.clear()
  dev.off()
  
  # Concatenate all results into one file with header
  all_res <- data.frame()
  for (mt in motifs_to_process) {
    strand <- mt$strand
    results_list <- motif_results_list[[strand]]
    for (i in seq_along(window_sizes)) {
      res <- results_list[[i]]
      res$window_size <- window_sizes[i]
      all_res <- rbind(all_res, res)
    }
  }
  output_base <- paste0("motif_results_all_", gsub("[^A-Za-z0-9]", "_", motif))
  output_file <- timestamped_filename(output_base, ".bed")
  write.table(
    all_res[, c("chrom", "start", "end", "strand", "observed", "expected", "sd", "z_score", "p_value", "window_size")],
    output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    append = FALSE
  )
  
  # Summary statistics per sequence
  for (mt in motifs_to_process) {
    strand <- mt$strand
    results_list <- motif_results_list[[strand]]
    for (i in seq_along(window_sizes)) {
      results <- results_list[[i]]
      summary_per_sequence <- results %>%
        group_by(chrom) %>%
        summarize(mean_p = mean(-log10(p_value)),
                  num_sig = sum(p_value < 0.05),
                  .groups = "drop")
      summary_df <- rbind(summary_df, data.frame(
        motif = motif,
        strand = strand,
        window_size = window_sizes[i],
        sequence = summary_per_sequence$chrom,
        mean_p = summary_per_sequence$mean_p,
        num_sig = summary_per_sequence$num_sig
      ))
    }
  }
}

# Save summary statistics
summary_file <- timestamped_filename("motif_summary_stats", ".txt")
write.table(summary_df, summary_file, sep = "\t", quote = FALSE, row.names = FALSE)
log_info("Analysis complete. Plots and results saved for each motif")
cat("Analysis complete. Plots and results saved for each motif.\n")
                          
