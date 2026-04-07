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
remotes::install_github("SDUT-Sysbio/AutoLD", build_vignettes = TRUE)
