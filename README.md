# AutoLD: Linkage Disequilibrium and Double Reduction Estimation in Autopolyploids

[![R-CMD-check](https://img.shields.io/badge/R_CMD_check-passing-brightgreen.svg)](https://github.com/SDUT-Sysbio/AutoLD)
[![Version](https://img.shields.io/badge/version-1.0.0-blue.svg)](https://github.com/SDUT-Sysbio/AutoLD)
[![License: GPL v3](https://img.shields.io/badge/License-GPL_v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

AutoLD is a high-performance, statistically rigorous R package designed for genomic analysis in autopolyploids (e.g., autotetraploids and autohexaploids). It provides an optimized Expectation-Maximization (EM) framework to simultaneously estimate genome-wide Linkage Disequilibrium (LD) and Double Reduction (DR) rates using dosage-based genotypic data.

## Key Features

- **Joint Estimation Model**: Accurately estimates fundamental 2nd-order LD (r and D' equivalents) while dynamically accounting for locus-specific double reduction during polyploid meiosis.
- **High-Order LD Statistics**: Captures complex multi-allelic co-segregations by evaluating higher-order LD components.
- **Multi-Ploidy Support**: Comprehensively supports both autotetraploid (4x) and autohexaploid (6x) genomic datasets.
- **Ultra-Fast Parallel Scanning**: Equipped with `AutoLD_optimized`, a sliding-window engine built on native multi-core architectures, capable of rapidly scanning chromosomal-scale data.
- **Built-in Simulations**: Includes robust simulation modules for power analysis and statistical validation.

## Installation

You can install the development version of AutoLD directly from GitHub. We highly recommend building the vignettes to access the comprehensive real-data tutorials.

```R
# Install remotes package if not already installed
if (!require("remotes", quietly = TRUE)) {
    install.packages("remotes")
}

# Install AutoLD from GitHub (SDUT-Sysbio)
remotes::install_github("SDUT-Sysbio/AutoLD", build_vignettes = TRUE).

```
## Quick Start

1. Computer Simulation
AutoLD provides built-in functions to simulate both autotetraploid and autohexaploid genotypes under specific LD and Double Reduction conditions. This is highly useful for power analysis and method validation.

```R
library(AutoLD)

# --- Example A: Autotetraploid (4x) Simulation ---
# Simulate a 2-locus dataset (500 individuals)
para <- c(pA = 0.6, pB = 0.4,alphaA = 1/6, alphaB = 1/7,Deab=0.05, DAb=0.02, DaB=0.01, DAB=0.002)
sim_data_4x <- sim_2locus_4x_cLD(para=para, n=500, type="full")

# Estimate LD and DR from the simulated 4x data
sim_result_4x <- auto4_est(sim_data_4x)
print(sim_result_4x)

# --- Example B: Autohexaploid (6x) Simulation ---
# Simulate a 2-locus dataset (500 individuals)
para <- c(pA = 0.6,pB = 0.4,alphaA = 1/6,alphaB = 1/5,Deab = 0.03,DAb = 0.002,DaB = 0.005,DAAb = 0.003,DaBB = 0.002,DAB = 0.002,DAAB = 0.003,DABB = 0.003,DAABB = 0.001)
sim_data_6x <- sim_2locus_6x_cLD(para=para, n=500, type="full")

# Estimate LD and DR from the simulated 6x data
sim_result_6x <- auto6_est(sim_data_6x)
print(sim_result_6x)

```
2. Real Data Analysis
Here is a minimal example of how to perform a genome-wide LD scan using the built-in autotetraploid Arabidopsis arenosa dataset.

```R
file_path <- system.file("extdata", "Arenosa_example.RDS", package = "AutoLD")
Arenosa <- readRDS(file_path)
geno_df <- Arenosa$Arenosa_gen
info_df <- Arenosa$Arenosa_gen_info

# Run the high-speed sliding window scan
results_chr1 <- AutoLD_optimized(
  chr_name  = "LR999451.1", 
  geno_df   = geno_df, 
  info_df   = info_df, 
  method    = "full", 
  window_kb = 200, 
  n_cores   = 4,
  ploidy    = 4
)

# View the core metrics
head(results_chr1[, c("POS1", "POS2", "Dist_bp", "aA", "Deab_r", "Deab_DS")])

```
### Documentation

Detailed documentation is embedded within the package. For a complete walkthrough of real data processing, statistical significance testing, and downstream visualizations, please refer to the package vignette:

```R
# Open the comprehensive HTML tutorial
vignette("AutoLD_RealData_Tutorial")
```
To view the manual for specific functions:

```R
?AutoLD_optimized
?sim_2locus_4x_cLD
?sim_2locus_6x_cLD

```
### Result Interpretation
The output data frame provides highly interpretable metrics for downstream population genetics analysis:

Dist_bp: Physical distance (bp), essential for LD decay modeling.

aA / aB: Estimated Double Reduction rates, revealing abnormal meiotic homologous pairing intensity.

Deab_r: Standardized LD on the correlation scale (analogous to r).

Deab_DS: Standardized LD on the D-prime scale (analogous to D').

P_Value: Rigorous Likelihood Ratio Test (LRT) statistics assessing the significance of both double reduction and linkage.

```
### Citation
If you use AutoLD in your research, please cite our upcoming publication:

Jiang, L. et al. (2026). AutoLD: A novel statistical framework for linkage disequilibrium and double reduction estimation in autopolyploids. (Submitted / In Preparation).

Contact & Authors
Jiang-Lab Systems Biology Group (SDUT-Sysbio)

Shandong University of Technology

Maintainer: Libo Jiang (libojiang@bjfu.edu.cn)

For bug reports, feature requests, or questions regarding usage, please open an issue on the GitHub Issues page.
