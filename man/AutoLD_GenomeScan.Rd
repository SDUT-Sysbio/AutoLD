\name{AutoLD_GenomeScan}
\alias{AutoLD_GenomeScan}
\alias{AutoLD_GenomeScan_can}
\title{Genome-Scale AutoLD Scan Within a Chromosome}
\description{
Evaluates marker pairs within a physical-distance window, or a supplied candidate-pair set, using the AutoLD estimator. Optional permutation testing can be triggered by conventional or higher-order screening thresholds.
}
\usage{
AutoLD_GenomeScan(chr_name, geno_df, info_df,
  window_kb = 500, ploidy = 4, type = "full",
  calc_norm = TRUE, calc_pvalue = FALSE,
  pvalue_trigger_r2 = NULL, pvalue_trigger_HO = NULL,
  B_perm = 200, n_cores = 1, min_N_eff = 30)

AutoLD_GenomeScan_can(chr_name, geno_df, info_df,
  window_kb = 500, ploidy = 4, type = "full",
  calc_norm = TRUE, calc_pvalue = FALSE,
  pvalue_trigger_r2 = NULL, pvalue_trigger_HO = NULL,
  B_perm = 200, n_cores = 1, min_N_eff = 30,
  candidate_pairs = NULL, force_pvalue = FALSE)
}
\arguments{
  \item{chr_name}{Chromosome identifier to analyze; must match values in \code{info_df$CHROM}.}
  \item{geno_df}{Marker-by-individual dosage data frame or matrix. Rows must correspond to rows of \code{info_df}.}
  \item{info_df}{Marker information with at least \code{CHROM} and \code{POS} columns.}
  \item{window_kb}{Maximum physical distance, in kb, between marker pairs in a full scan.}
  \item{ploidy}{Even somatic ploidy.}
  \item{type}{\code{"full"} or \code{"partial"}.}
  \item{calc_norm}{Whether to calculate normalized tensor components.}
  \item{calc_pvalue}{Whether to perform permutation calibration for eligible pairs.}
  \item{pvalue_trigger_r2}{Optional threshold on \code{r2_comp} that triggers permutation testing.}
  \item{pvalue_trigger_HO}{Optional threshold on the maximum absolute higher-order normalized component that triggers permutation testing.}
  \item{B_perm}{Number of permutations per tested pair.}
  \item{n_cores}{Number of CPU cores. Unix-like systems use forked parallelism; Windows uses a cluster.}
  \item{min_N_eff}{Minimum number of complete two-locus observations required for a pair.}
  \item{candidate_pairs}{Optional data frame of candidate pairs. It should contain \code{POS_A} and \code{POS_B}, and optionally \code{CHROM} or \code{CHR}.}
  \item{force_pvalue}{Logical; in candidate mode, force permutation testing regardless of screening thresholds.}
}
\value{
A data frame with marker positions, distance, effective sample size, ESD estimates, conventional LD summaries, raw and optionally normalized \eqn{D_{rs}} components, and permutation-calibrated P values when requested.
}
\examples{
\dontrun{
scan <- AutoLD_GenomeScan(
  chr_name = unique(info$CHROM)[1],
  geno_df = geno,
  info_df = info,
  window_kb = 500,
  ploidy = 4,
  calc_norm = TRUE,
  calc_pvalue = FALSE,
  n_cores = 4,
  min_N_eff = 50
)
}
}
