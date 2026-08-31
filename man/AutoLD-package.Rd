\name{AutoLD-package}
\alias{AutoLD-package}
\alias{AutoLD}
\docType{package}
\title{AutoLD: Dosage-Dependent Linkage Disequilibrium in Autopolyploid Genomes}
\description{
AutoLD implements a population-genetic framework for reconstructing two-locus gametic dosage distributions in autopolyploids and decomposing them into orthogonal dosage-dependent linkage disequilibrium components. The package includes single-pair estimation, simulation, genome scanning, fixed-margin normalization, and marginal-preserving permutation calibration.
}
\details{
The main user-facing entry points are \code{AutoLD()} for a single marker pair, \code{simulate_from_eta()} and \code{evaluate_AutoLD_pipeline()} for simulation, and \code{AutoLD_GenomeScan()} for empirical chromosome-scale analyses.

In AutoLD, higher-order LD refers to higher-order dosage contrasts within a two-locus distribution. It does not denote association among more than two loci or epistasis.
}
\seealso{
\code{\link{AutoLD}}, \code{\link{AutoLD_GenomeScan}}, \code{\link{simulate_from_eta}}, \code{\link{evaluate_AutoLD_pipeline}}
}
\examples{
D <- matrix(0, 2, 2)
D[1, 1] <- 0.05
D[1, 2] <- 0.03

sim <- simulate_from_eta(
  v_ploidy = 4, pA = 0.6, pB = 0.5,
  eta1_A = 0.15, eta1_B = 0.05,
  D_matrix = D, n = 500, seed = 1
)

fit <- AutoLD(sim$zygote_A, sim$zygote_B, ploidy = 4,
              calc_norm = FALSE, permutations = 0)
print(fit)
}
