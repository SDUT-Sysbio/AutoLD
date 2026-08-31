# AutoLD

**AutoLD** is an R/C++ framework for resolving dosage-dependent two-locus linkage disequilibrium (LD) in autopolyploid genomes.

The current `main` implementation is organized around the substantially reconstructed framework used in the current manuscript. Higher-order LD in AutoLD refers to **orthogonal dosage contrasts within a two-locus gametic dosage distribution**, not to interactions among three or more loci, epistasis, or regulatory networks.

The implementation corresponding to the earlier bioRxiv/preprint version should be retained separately under the Git tag/release **`preprint-v1`**.

## Core functionality

- ploidy-aware segregation transition kernel and fixed-point equilibrium solver;
- EM reconstruction of two-locus gametic dosage distributions;
- orthogonal dosage basis and `D_rs` tensor decomposition;
- effective marginal segregation-distortion estimation;
- fixed-margin feasible normalization using linear programming;
- marginal-preserving permutation calibration;
- 4x/6x simulation workflows;
- chromosome-scale empirical genome scans and candidate-only permutation testing.

## Installation

```r
install.packages(c("Rcpp", "data.table", "lpSolve", "remotes"))
remotes::install_github("SDUT-Sysbio/AutoLD")
library(AutoLD)
```

## Single-pair analysis

```r
fit <- AutoLD(
  dosage_A = c(0, 1, 2, 3, 4),
  dosage_B = c(0, 1, 2, 3, 4),
  ploidy = 4,
  type = "full",
  calc_norm = TRUE,
  permutations = 0
)
print(fit)
```

For practical workflows, see:

- `examples/01_simulation_example.R`
- `examples/02_real_data_Aarenosa.R`

## Empirical input convention

For `AutoLD_GenomeScan()`, `geno_df` is a marker-by-individual dosage table and `info_df` contains matching marker rows with at least `CHROM` and `POS` columns. Tetraploid full-dosage data should use integer states 0-4, with missing genotypes represented as `NA`.

## Reproducibility note

This first package encapsulation preserves the numerical/statistical logic of the supplied core analysis scripts and primarily reorganizes them into package form. Additional figure-generation, supplementary validation, and dataset-specific preprocessing scripts can be added incrementally without changing the public core API.

## Help and documentation

After installation, the main help pages can be opened directly in R:

```r
help(package = "AutoLD")
?AutoLD
?AutoLD_GenomeScan
?simulate_from_eta
?evaluate_AutoLD_pipeline
?AutoLD_Permutation_Test
?thin_snps_by_dist
```

The package also includes lower-level documentation for the segregation kernel,
orthogonal LD tensor, and core encoded-data estimator. The scripts in `examples/`
show the main simulation and empirical-analysis workflows.
