\name{estimate_LD_decoupled}
\alias{estimate_LD_decoupled}
\title{Core AutoLD Estimator for Encoded Two-Locus Data}
\description{
Low-level estimator used by the high-level AutoLD interfaces. It reconstructs the joint gametic distribution by EM, computes marginal distributions and orthogonal LD components, and optionally applies fixed-margin normalization.
}
\usage{
estimate_LD_decoupled(geno, c_ploidy, type = "full", calc_norm = TRUE)
}
\arguments{
  \item{geno}{Data frame with integer column \code{class} encoding the two-locus observation.}
  \item{c_ploidy}{Gametic ploidy.}
  \item{type}{\code{"full"} or \code{"partial"}.}
  \item{calc_norm}{Whether to calculate sign-specific feasible normalization.}
}
\value{
A list containing effective sample size, marginal segregation estimates, allele frequencies, conventional dosage LD summaries, EM convergence status, raw \eqn{D_{rs}} components, and optional normalized components and feasible bounds.
}
\seealso{\code{\link{AutoLD}}}
