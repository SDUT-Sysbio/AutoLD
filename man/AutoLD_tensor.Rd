\name{AutoLD_tensor}
\alias{build_orthogonal_basis}
\alias{extract_Drs_tensor}
\alias{convert_b_to_a_4x}
\alias{calc_Deab_comp}
\title{Orthogonal Dosage Basis and LD Tensor Utilities}
\description{
Construct probability-weighted orthogonal dosage bases, project a two-locus gametic probability table onto the AutoLD tensor, convert the autotetraploid linear-mixture parameter to the classical scale, and calculate a conventional dosage-based gametic LD summary.
}
\usage{
build_orthogonal_basis(Q, c_ploidy)
extract_Drs_tensor(P_AB, c_ploidy)
convert_b_to_a_4x(b)
calc_Deab_comp(P_AB)
}
\arguments{
  \item{Q}{Gametic marginal dosage probability vector.}
  \item{c_ploidy}{Gametic ploidy.}
  \item{P_AB}{Joint two-locus gametic dosage probability matrix.}
  \item{b}{Autotetraploid linear-mixture parameter.}
}
\value{
A basis matrix, LD tensor matrix, converted scalar parameter, or conventional gametic LD summary, depending on the function.
}
\details{
The orthogonal basis is constructed with respect to the supplied marginal distribution. Components \eqn{D_{rs}} are therefore contrasts conditional on the reconstructed marginal dosage structure.
}
