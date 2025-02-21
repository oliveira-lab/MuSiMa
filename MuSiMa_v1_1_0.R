# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript scriptname.R FASTA1 FASTA2 ... FASTAN Motif")
}

# Extract FASTA file names and motif from arguments
fasta_files <- args[-length(args)]  # FASTA files
motif <- toupper(args[length(args)])  # DNA motif

# Check if files exist
for (file in fasta_files) {
  if (!file.exists(file)) {
    stop(paste("Error: The file", file, "does not exist. Please provide valid file paths."))
  }
}

# Check and install required packages
required_packages <- c("seqinr", "circlize", "dplyr", "Biostrings", "ComplexHeatmap", "parallel")
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
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# Parameters
step_size <- 10000
window_sizes <- c(500000, 400000, 300000, 200000, 100000, 90000, 80000, 70000, 60000, 50000)

# Function to load DNA sequence from FASTA file
load_fasta_sequence <- function(fasta_file) {
  seq <- read.fasta(fasta_file, as.string = TRUE, forceDNAtolower = FALSE)[[1]]
  seq_string <- toupper(as.character(seq))
  return(DNAString(seq_string))  # Convert to DNAString
}

# Function to find degenerate motif positions
find_motif_positions <- function(dna_seq, motif) {
  matches <- matchPattern(motif, dna_seq, fixed = FALSE)
  if (length(matches) == 0) return(data.frame(start = integer(), end = integer()))
  starts <- start(matches)
  ends <- end(matches)
  return(data.frame(start = starts, end = ends))
}

# Function to precompute motif counts per step
precompute_motif_counts <- function(motif_positions, sequence_length, step_size) {
  bins <- floor((sequence_length - 1) / step_size) + 1
  counts <- numeric(bins)
  for (i in seq_along(motif_positions$start)) {
    bin <- floor((motif_positions$start[i] - 1) / step_size) + 1
    if (bin <= bins) counts[bin] <- counts[bin] + 1
  }
  return(counts)
}

# Function to calculate expected motif occurrences for degenerate motifs
calc_expected_motif <- function(dna_seq, motif, window_size) {
  motif_length <- nchar(motif)
  total_possible_positions <- length(dna_seq) - motif_length + 1
  observed_total <- length(matchPattern(motif, dna_seq, fixed = FALSE))
  p_motif <- observed_total / total_possible_positions
  expected <- p_motif * (window_size - motif_length + 1)
  return(expected)
}

# Function to calculate z-score
calc_z_score <- function(observed, expected, window_size, motif_length) {
  p <- expected / (window_size - motif_length + 1)
  variance <- (window_size - motif_length + 1) * p * (1 - p)
  sd <- sqrt(variance)
  z_score <- if (sd == 0) 0 else (observed - expected) / sd
  return(z_score)
}

# Load sequences
dna_sequences <- lapply(fasta_files, load_fasta_sequence)

# Precompute motif positions and counts
motif_positions_list <- lapply(dna_sequences, find_motif_positions, motif)
motif_counts_list <- mapply(precompute_motif_counts, motif_positions_list, sapply(dna_sequences, length), MoreArgs = list(step_size), SIMPLIFY = FALSE)

# Create chromosome data frame
chr_lengths <- sapply(dna_sequences, length)
df <- data.frame(
  name = paste0("chr", seq_along(fasta_files)),
  start = cumsum(c(1, head(chr_lengths, -1))),
  end = cumsum(chr_lengths)
)

# Open PDF output
pdf("musima_plot.pdf", width = 15, height = 10)
circos.par("track.height" = 0.05, "gap.degree" = 5, "start.degree" = 90)
circos.genomicInitialize(df)

# Color scale
col_fun <- colorRamp2(c(-10, -5, 0, 5, 10), c("blue", "navy", "green", "red", "darkred"))

# Process window sizes in parallel
process_window_size <- function(window_size) {
  cat(sprintf("Processing window size: %d\n", window_size))
  results <- data.frame()
  offset <- 0
  
  for (chr_index in seq_along(fasta_files)) {
    dna_seq <- dna_sequences[[chr_index]]
    sequence_length <- length(dna_seq)
    motif_counts <- motif_counts_list[[chr_index]]
    num_windows <- floor((sequence_length - window_size) / step_size) + 1
    
    expected <- calc_expected_motif(dna_seq, motif, window_size)
    steps_per_window <- ceiling(window_size / step_size)
    
    for (i in 1:num_windows) {
      window_start <- (i - 1) * step_size + 1 + offset
      window_end <- window_start + window_size - 1
      if (window_end > sequence_length + offset) next
      
      start_bin <- (i - 1) + 1
      end_bin <- min(start_bin + steps_per_window - 1, length(motif_counts))
      observed <- sum(motif_counts[start_bin:end_bin])
      
      z_score <- calc_z_score(observed, expected, window_size, nchar(motif))
      results <- rbind(results, data.frame(
        chrom = df$name[chr_index],
        start = window_start,
        end = window_end,
        value1 = z_score
      ))
    }
    offset <- offset + sequence_length
  }
  return(results)
}

# Run in parallel
num_cores <- detectCores() - 1
results_list <- mclapply(window_sizes, process_window_size, mc.cores = num_cores)

# Plot tracks
for (results in results_list) {
  circos.genomicTrack(results, stack = TRUE, panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = col_fun(value[[1]]), border = NA, ...)
  })
}

# Add legend
draw(Legend(col_fun = col_fun, title = "Motif Z-Score", direction = "vertical"), x = unit(0.9, "npc"), y = unit(0.5, "npc"))

circos.clear()
dev.off()

cat("Analysis complete. Plot saved as 'musima_plot.pdf'.\n")
