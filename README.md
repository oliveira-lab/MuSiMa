# MuSiMa: Multiscale Signal Mapping

This repository contains an R script (`musima.R`) designed to analyze DNA sequences from FASTA files for the occurrence of a specified motif (including degenerate motifs) and visualize the results as a circular plot of z-scores. MuSiMa stands for **Multiscale Signal Mapping**, reflecting its ability to map motif signals across multiple scales (window sizes).

## Purpose
MuSiMa:
- Identifies positions of a user-specified motif across multiple DNA sequences provided in FASTA format.
- Calculates observed versus expected motif occurrences in sliding windows of varying sizes (50 to 500 kb) with a 10 kb step.
- Computes z-scores to assess statistical significance of motif enrichment or depletion.
- Generates a circular visualization saved as a PDF file (`musima_plot.pdf`).

This tool is particularly suited for genomic analyses where understanding motif distribution across chromosomes or contigs is of interest.

## Assumptions
- **Input Files**: MuSiMa assumes that all provided FASTA files are valid, single-sequence DNA files. Multi-sequence FASTA files are not supported (only the first sequence is used).
- **Motif**: The motif is a valid DNA string, potentially including IUPAC degenerate codes (e.g., "ATGC" or "AT[G|C]C"). It is case-insensitive (converted to uppercase internally).
- **Environment**: An R installation with internet access is required to install missing packages. MuSiMa uses Bioconductor packages (`Biostrings`, `ComplexHeatmap`) and CRAN packages (`seqinr`, `circlize`, etc.).
- **Sequence Length**: Sequences should be long enough to accommodate the largest window size (500 kb by default) for meaningful analysis. Subsequent versions will allow the user to finetune these parameters.
- **Hardware**: Parallel processing is utilized, assuming a multi-core system (uses all but one core).

## Dependencies
MuSiMa automatically checks for and installs the following R packages if they are not already present:
- **CRAN**: `seqinr`, `circlize`, `dplyr`, `parallel`
- **Bioconductor**: `Biostrings`, `ComplexHeatmap`, `BiocManager` (for installation)

## Installation and running
1. Ensure R is installed on your system (version 4.0 or higher recommended).
2. Clone this repository:
   ```bash
   git clone https://github.com/oliveira-lab/musima.git
   cd musima
   Rscript musima.R FASTA1 FASTA2 ... FASTAN Motif

## Output
A PDF file named musima_plot.pdf containing a circular plot of z-scores across chromosomes for multiple window sizes (50–500 kb).

## License and citing
This project is licensed under the MIT License. See the LICENSE file for details. Please cite MuSiMa by including the link to https://github.com/oliveira-lab/musima.git.
