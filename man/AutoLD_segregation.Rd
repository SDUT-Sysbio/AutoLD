\name{AutoLD_segregation}
\alias{transition_kernel_biallelic}
\alias{solve_equilibrium_Q}
\alias{get_equilibrium_Q_one_eta}
\alias{get_equilibrium_Q_4x_closed}
\alias{validate_4x_solver}
\alias{project_ESD_general}
\alias{project_ESD_4x}
\alias{project_ESD_6x}
\title{Segregation Kernels, Equilibrium Margins, and ESD Projection}
\description{
Functions implementing the ploidy-aware biallelic meiotic transition kernel, fixed-point equilibrium gametic distributions, autotetraploid closed-form validation, and projection of empirical gametic margins onto the segregation-equilibrium family.
}
\usage{
transition_kernel_biallelic(v, eta)
solve_equilibrium_Q(v, p, trans_matrix, tol = 1e-12, max_iter = 2000)
get_equilibrium_Q_one_eta(v, p, eta1)
get_equilibrium_Q_4x_closed(p, a)
validate_4x_solver()
project_ESD_general(Q_emp)
project_ESD_4x(Q_emp)
project_ESD_6x(Q_emp)
}
\arguments{
  \item{v}{Even somatic ploidy.}
  \item{eta}{Segregation vector on the probability simplex.}
  \item{p}{Reference-allele frequency.}
  \item{trans_matrix}{Transition matrix returned by \code{transition_kernel_biallelic()}.}
  \item{tol}{Fixed-point convergence tolerance.}
  \item{max_iter}{Maximum fixed-point iterations.}
  \item{eta1}{Probability assigned to one IBDR pair in the current one-parameter helper.}
  \item{a}{Classical autotetraploid double-reduction parameter used in the closed-form equilibrium expression.}
  \item{Q_emp}{Empirical/reconstructed gametic dosage distribution.}
}
\value{
Depending on the function, a transition matrix, equilibrium gametic distribution, validation output, or a list containing the projected effective segregation-distortion estimate, allele frequency, and projection loss.
}
\details{
For empirical interpretation, ESD is an effective marginal distortion parameter and should not automatically be interpreted as a direct cytological double-reduction rate.
}
