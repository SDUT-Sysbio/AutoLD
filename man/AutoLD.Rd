\name{AutoLD}
\alias{AutoLD}
\alias{print.AutoLD}
\title{Fit AutoLD to a Single Pair of Autopolyploid Markers}
\description{
Fits the AutoLD model to two dosage-resolved loci, reconstructs the two-locus gametic distribution by EM, estimates orthogonal LD components, and optionally performs fixed-margin normalization and marginal-preserving permutation calibration.
}
\usage{
AutoLD(dosage_A, dosage_B, ploidy = 4L,
       type = c("full", "partial"), calc_norm = TRUE,
       permutations = 0L, seed = NULL)

\method{print}{AutoLD}(x, ...)
}
\arguments{
  \item{dosage_A, dosage_B}{Integer allele-dosage vectors for loci A and B. For full data, values must lie between 0 and the somatic ploidy. Missing values are allowed and are removed pairwise.}
  \item{ploidy}{Even somatic ploidy. The current implementation is primarily intended for autotetraploid and autohexaploid analyses.}
  \item{type}{Observation type: \code{"full"} for resolved dosage states or \code{"partial"} for collapsed 0/intermediate/ploidy classes.}
  \item{calc_norm}{Logical; whether to compute sign-specific feasible normalization of each LD component under fixed reconstructed margins.}
  \item{permutations}{Number of marginal-preserving permutations. Set to 0 to skip permutation testing. Values of at least 999 are recommended for inference.}
  \item{seed}{Optional random seed used for permutation testing.}
  \item{x}{An object of class \code{"AutoLD"}.}
  \item{...}{Additional arguments passed to the print method; currently unused.}
}
\value{
An object of class \code{"AutoLD"}, containing the call, ploidy, observation type, effective sample size, an \code{estimate} list with ESD, first-order LD, tensor components and optional normalized components, and a \code{permutation} list when permutation testing is requested.
}
\details{
The first-order component \eqn{D_{11}} represents the primary dosage association. Components with \eqn{r>1} or \eqn{s>1} represent higher-order dosage contrasts within the same two-locus gametic distribution.
}
\examples{
D4 <- matrix(0, 2, 2)
D4[1, 1] <- 0.05
D4[1, 2] <- 0.03
sim4 <- simulate_from_eta(4, 0.6, 0.5, 0.15, 0.05, D4,
                          n = 1000, type = "full", seed = 123)
fit4 <- AutoLD(sim4$zygote_A, sim4$zygote_B,
               ploidy = 4, calc_norm = TRUE, permutations = 0)
fit4
}
