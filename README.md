# MuSiMa: Multiscale Signal Mapping
This repository contains an R script (`musima.R`) designed to detect enrichment or depletion of DNA motifs or genomic signals in sliding windows across DNA sequences in FASTA format, and visualize the results as a circular plot of z-scores. This tool is ideal for genomic motif scanning (e.g., TF binding sites, methylation target sites) or region enrichment (e.g., ChIP-seq peaks, variants) to reveal non-random distributions.

## Features
MuSiMa supports two modes:

**'Motif' mode**: Scans for motif occurrences (exact or degenerate) in sliding windows and computes statistical significance using probabilistic models (uniform or Markov chain-based backgrounds).
- Identifies positions of user-specified motifs across multiple DNA sequences provided in single or multi-sequence FASTA file.
- Calculates observed versus expected motif occurrences in sliding windows of varying sizes (user defined) with a given step (also user defined). Expected values use a probabilistic null model: i) "uniform" (independent bases, adjusted for sequence composition with +1 pseudocounts); or ii) "markov" (smoothed transition probabilities up to user-specified order, max < motif length).
- Computes z-scores and FDR-adjusted two-tailed p-values (normal/Poisson approximation with continuity correction). Palindromic motifs produce one set of tracks; non-palindromic produce separate sense (+) and antisense (-) tracks.

**'Signal' mode**: Analyzes enrichment or depletion of genomic signals via randomization-based null distributions.
- Analyzes user-defined genomic signals (from a 4-column BED file: chrom/ID, start, end, strand) for enrichment in sliding windows across a single genome FASTA. Computes observed signals per window/strand, expected values via randomization (user-specified iterations, >100 recommended) as well as Z-scores and FDR-adjusted two-tailed empirical p-values.

## Assumptions
- **Input Files**: MuSiMa assumes that all provided FASTA files are valid, single- or multi-sequence DNA files. BED files must contain exactly 4 columns (chrom/ID, start, end, strand as +/−).
- **Motif**: The motif is a valid DNA string, potentially including IUPAC degenerate codes (e.g., "ATGC" or "AT[G|C]C"). It is case-insensitive (converted to uppercase internally).
- **Environment**: An R installation with internet access is required to install missing packages. MuSiMa uses Bioconductor packages (`Biostrings`, `ComplexHeatmap`) and CRAN packages (`seqinr`, `circlize`, etc.).
- **Sequence Length**: Sequences should be long enough to accommodate the largest window size (10 times smaller as a rule of thumb) for meaningful analysis.
- **Hardware**: Parallel processing is utilized, assuming a multi-core system (uses all but one core).

## Dependencies
MuSiMa automatically checks for and installs the following R packages if they are not already present:
- **CRAN**: `seqinr`, `circlize`, `dplyr`, `parallel`, `optparse`, `logger`
- **Bioconductor**: `Biostrings`, `ComplexHeatmap`, `BiocManager`

## Installation and running
1. Ensure R is installed on your system (version 4.0 or higher recommended).
2. Clone this repository:
   ```bash
   git clone https://github.com/oliveira-lab/musima.git
   cd musima
   ```
3. Command-line options:
   - `-f, --fasta`: Comma-separated list of FASTA files (e.g., file1.fa,file2.fa) or a single multi-sequence FASTA file. Single file for 'Signal' mode.
   - `-m, --motif`: Motif string (e.g., GATC,CCWGG) or file path with motifs (one per line).
   - `-w, --windows`: Comma-separated list of window sizes (e.g., 500000,400000).
   - `-s, --step`: Step size for sliding windows (e.g., 10000) (default 1).
   - `-t, --method`: Method for expected occurrences: 'uniform' or 'markov'.
   - `-o, --order`: Markov method order (integer or 'NA' for uniform method).
   - `-c, --cores`: Number of CPU cores for parallel processing (default: all available minus one).
   - `-b, --bed`: BED file (4 cols: ID, start, end, strand). Triggers Signal mode.
   - `-r, --randomization`: Randomization iterations (>100 recommended).
   - `--plotrange`: Z-score plot range (e.g., "-5,5"; default: -5,5).

    Notes:
- Signal mode ignores --motif, --method, --order
- Motif mode ignores --bed, --randomization
- Only one set of tracks is printed for palindromic motifs.

4. **'Motif'** mode example:
   - Distribution of "GATC" in *E. coli* MG1655 and *R. solanacearum* GMI1000 using as background a markov model of order 2:
     ```bash
     Rscript musima.R -f coli.fa,solanacearum.fa -m GATC -w 500000,250000,100000,90000,80000,70000,50000 -s 10000 -t markov -o 2 -c 4 --plotrange=-5,5
     ```

     ![MuSiMa Output](/test/musima_plot_GATC.jpg "GATC distribution across *E. coli* MG1655 and *R. solanacearum* GMI1000")

   - Distribution of "GTWWAC" and "RAATTY" in *E. coli* MG1655 and *R. solanacearum* GMI1000 using as background a markov model of order 5:
     ```bash
     Rscript musima.R -f coli.fa,solanacearum.fa -m GTWWAC,RAATTY -w 500000,250000,100000,90000,80000,70000,50000 -s 10000 -t markov -o 5 -c 4 --plotrange=-5,5
     ```

     ![MuSiMa Output](/test/musima_plot_multimotif.jpg "GTWWAC and RAATTY distribution across *E. coli* MG1655 and *R. solanacearum* GMI1000")


5. **'Signal'** mode example:
   - TFBS in *C. difficile* 630 (from [Oliveira *et al* (2020) Nat. Microbiol](https://www.nature.com/articles/s41564-019-0613-4)).
     ```bash
     Rscript musima.R -f difficile.fasta --bed TFBS.bed -w 500000,250000,100000,90000,80000,70000,50000 -s 10000 --randomization 200 -c 4 --plotrange=-5,5
     ```

     ![MuSiMa Output](/test/musima_plot_TFBS.jpg "TFBS distribution across *C. difficile* 630")

## Output
Multiple files are produced:
- Circular visualizations (`musima_plot_*.pdf` with stacked tracks by window size/strand and color-coded by z-score).
- Per-window stats (`results_all_*.bed`: ID, start, end, strand, observed, expected, sd, z_score, p_value, window_size). 
- Summary (`summary_stats.txt`: window_size, ID, strand, mean -log10(p), significant windows per chrom/motif/region).
- Run logs (`musima.log`).

## License and citing
This project is licensed under a GPL-3.0 License. See the LICENSE file for details. Please cite MuSiMa by including the link to https://github.com/oliveira-lab/musima.git.
