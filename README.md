# MuSiMa: Multiscale Signal Mapping

This repository contains an R script (`musima.R`) designed to analyze DNA sequences from FASTA files for the occurrence of a specified motif (including degenerate motifs) and visualize the results as a circular plot of z-scores. MuSiMa stands for **Multiscale Signal Mapping**, reflecting its ability to map motif signals across multiple scales (window sizes).

## Purpose
MuSiMa:
- Identifies positions of user-specified motifs across multiple DNA sequences provided in FASTA format or a single multi-sequence FASTA file.
- Calculates observed versus expected motif occurrences in sliding windows of varying sizes (user defined) with a given step (also user defined). Innermost layer corresponds to the smallest window size. The expected motif occurrence is computed either as i) the product of the genome-wide motif frequency and the number of potential motif sites within each window, under a null model of uniform random distribution across the sequence; or ii) using a Markov chain of a specified order through sequence transition probabilities (maximum order allowed is L-2, where L is the motif's length). If motif is palindromic, only one set of layers is produced. Otherwise two sets of layers will be produced (sense + antisense strands).
- Computes z-scores to assess statistical significance of motif enrichment or depletion.
- Generates circular visualizations saved as a PDF file (`musima_plot_MOTIF.pdf`), prints the raw z-score values for each window (`results_list_MOTIF.txt`), and a summary file (`summary_stats.txt`) providing mean z-scores and significant (|z|>2) window counts per motif and FASTA file.

This tool is particularly suited for genomic analyses where understanding motif distribution across chromosomes or contigs is of interest.

## Assumptions
- **Input Files**: MuSiMa assumes that all provided FASTA files are valid, single-sequence DNA files. Multi-sequence FASTA files are also supported.
- **Motif**: The motif is a valid DNA string, potentially including IUPAC degenerate codes (e.g., "ATGC" or "AT[G|C]C"). It is case-insensitive (converted to uppercase internally).
- **Environment**: An R installation with internet access is required to install missing packages. MuSiMa uses Bioconductor packages (`Biostrings`, `ComplexHeatmap`) and CRAN packages (`seqinr`, `circlize`, etc.).
- **Sequence Length**: Sequences should be long enough to accommodate the largest window size (10 times smaller as a rule of thumb) for meaningful analysis.
- **Hardware**: Parallel processing is utilized, assuming a multi-core system (uses all but one core).

## Dependencies
MuSiMa automatically checks for and installs the following R packages if they are not already present:
- **CRAN**: `seqinr`, `circlize`, `dplyr`, `parallel`, `optparse`
- **Bioconductor**: `Biostrings`, `ComplexHeatmap`, `BiocManager`

## Installation and running
1. Ensure R is installed on your system (version 4.0 or higher recommended).
2. Clone this repository:
   ```bash
   git clone https://github.com/oliveira-lab/musima.git
   cd musima
3. Command-line options
   ```TXT
   -f, --fasta: Comma-separated list of FASTA files (e.g., file1.fa,file2.fa) or a single multi-sequence FASTA file.
   -m, --motif: Motif string (e.g., GATC,CCWGG) or file path with motifs (one per line).
   -w, --windows: Comma-separated list of window sizes (e.g., 500000,400000).
   -s, --step: Step size for sliding windows (e.g., 10000).
   -t, --method: Method for expected occurrences: 'uniform' or 'markov'.
   -o, --order: Markov method order (integer, ) or 'NA' for uniform method.
   -c, --cores: Number of CPU cores for parallel processing (default: all available minus one).
5. Run it as:
   ```bash
   Rscript musima.R -f file1.fa,file2.fa -m Motif1,Motif2 -w window_size1,window_size2,...,window_sizeN -s step_size -t uniform -o NA -c 4

## Output
A PDF file named musima_plot_MOTIF.pdf containing a circular plot of z-scores across chromosomes for multiple window sizes and a given step size (user defined). Auxiliary results_list_MOTIF.txt and summary_stats.txt files are also produced. If motif is palindromic, only one set of layers is produced. Otherwise two sets of layers will be produced (sense + antisense strands).

![Output Musima](/test/musima_plot_CAAAAA.jpg "CAAAAA distribution across E. coli MG1655 and R. solanacearum GMI1000")

In this example, we produced a circular plot of z-scores for CAAAAA over / under abundance in the main chromosomes of <em>E. coli</em> MG1655 and <em>R. solanacearum</em> GMI1000 using window sizes 500000, 250000, 100000, 90000, 80000, 70000, 50000, a step size of 10000, and a 'uniform' background model.

## License and citing
This project is licensed under a GPL-3.0 License. See the LICENSE file for details. Please cite MuSiMa by including the link to https://github.com/oliveira-lab/musima.git.
