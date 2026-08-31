\name{AutoLD_Permutation_Test}
\alias{AutoLD_Permutation_Test}
\title{Marginal-Preserving Permutation Calibration for AutoLD}
\description{
Calibrates first-order, higher-order, and full-tensor AutoLD statistics by permuting one locus across individuals while preserving the single-locus dosage distributions.
}
\usage{
AutoLD_Permutation_Test(geno, c_ploidy,
  type = c("full", "partial"), B = 1000L,
  strata = NULL, seed = NULL, sd_tol = 1e-8,
  min_convergence_rate = 0.95, max_attempt_factor = 2,
  obs_res = NULL)
}
\arguments{
  \item{geno}{Data frame containing integer column \code{class}, encoding the two-locus observed state.}
  \item{c_ploidy}{Gametic ploidy, equal to half the somatic ploidy.}
  \item{type}{\code{"full"} or \code{"partial"}.}
  \item{B}{Target number of successfully converged permutations. At least 999 is recommended for empirical inference.}
  \item{strata}{Optional vector defining strata within which locus B is permuted.}
  \item{seed}{Optional random seed.}
  \item{sd_tol}{Minimum null standard deviation required for a component to be considered estimable.}
  \item{min_convergence_rate}{Minimum acceptable fraction of successful permutation fits.}
  \item{max_attempt_factor}{Maximum attempted permutations as a multiple of \code{B}.}
  \item{obs_res}{Optional previously fitted AutoLD estimate for the observed data, avoiding redundant estimation.}
}
\value{
A list containing observed tensor components, calibrated Z statistics, \code{P_D11}, omnibus higher-order and full-tensor P values, component-wise raw and MaxT-adjusted P values, null moments, estimability indicators, and convergence diagnostics.
}
\details{
The permutation destroys cross-locus association while retaining the empirical dosage distribution at each locus. For each successful permutation, the AutoLD estimation pipeline is rerun.
}
