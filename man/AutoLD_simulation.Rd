\name{AutoLD_simulation}
\alias{simulate_unified_zygotes}
\alias{simulate_from_eta}
\alias{evaluate_AutoLD_pipeline}
\title{Simulation Utilities for AutoLD}
\description{
Generate two-locus autopolyploid data from specified gametic margins and LD tensors, or from allele frequencies and segregation parameters, and evaluate AutoLD over repeated simulations.
}
\usage{
simulate_unified_zygotes(QA, QB, D_matrix, c_ploidy,
  n = 500, type = "full", hwd_f = 0, seed = NULL)

simulate_from_eta(v_ploidy, pA, pB, eta1_A, eta1_B, D_matrix,
  n = 500, type = "full", hwd_f = 0, seed = NULL)

evaluate_AutoLD_pipeline(reps = 100, n = 500, c_ploidy = 2,
  type = "full", pA = 0.5, pB = 0.5, eta_A = 0,
  eta_B = 0, D_true_matrix = NULL, hwd_f = 0,
  calc_norm = FALSE, calc_pvalue = FALSE, B_perm = 200,
  alpha = 0.05, n_cores = 1, seed = 123)
}
\arguments{
  \item{QA, QB}{Gametic marginal dosage distributions for loci A and B.}
  \item{D_matrix, D_true_matrix}{Matrix of true orthogonal LD components. Dimensions are \code{c_ploidy} by \code{c_ploidy}.}
  \item{c_ploidy}{Gametic ploidy.}
  \item{v_ploidy}{Somatic ploidy.}
  \item{pA, pB}{Reference-allele frequencies at loci A and B.}
  \item{eta1_A, eta1_B, eta_A, eta_B}{Single-IBDR segregation-vector component used to generate the theoretical marginal distributions in the current 4x/6x simulation interface.}
  \item{n}{Number of simulated individuals per replicate.}
  \item{type}{\code{"full"} or \code{"partial"}.}
  \item{hwd_f}{Probability that the second gamete is copied from the first gamete, used to induce a controlled Hardy-Weinberg disequilibrium perturbation in simulations.}
  \item{seed}{Random seed.}
  \item{reps}{Number of simulation replicates.}
  \item{calc_norm}{Whether to calculate normalized tensor components.}
  \item{calc_pvalue}{Whether to perform permutation testing in each replicate.}
  \item{B_perm}{Number of permutations per replicate when \code{calc_pvalue = TRUE}.}
  \item{alpha}{Significance level used to summarize empirical power or false-positive rates.}
  \item{n_cores}{Number of CPU cores used across simulation replicates.}
}
\value{
\code{simulate_unified_zygotes()} and \code{simulate_from_eta()} return a list containing observed classes, zygotic dosages, true gametic distribution, true margins, and true tensor components. \code{evaluate_AutoLD_pipeline()} invisibly returns converged replicate results and a component-level summary table.
}
\examples{
D4 <- matrix(0, 2, 2)
D4[1, 1] <- 0.05
D4[1, 2] <- 0.03
sim4 <- simulate_from_eta(4, 0.6, 0.5, 0.15, 0.05, D4,
                          n = 1000, seed = 123)
head(sim4$zygote_A)

\dontrun{
evaluate_AutoLD_pipeline(
  reps = 100, n = 2000, c_ploidy = 2,
  pA = 0.6, pB = 0.5,
  eta_A = 0.15, eta_B = 0.05,
  D_true_matrix = D4, n_cores = 4
)
}
}
