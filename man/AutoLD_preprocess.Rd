\name{AutoLD_preprocess}
\alias{thin_snps_by_dist}
\alias{add_autold_screening}
\title{Preprocessing and Screening Utilities for AutoLD Genome Scans}
\description{
Utilities for physical-distance SNP thinning and for assigning conventional/high-order screening classes to AutoLD genome-scan results.
}
\usage{
thin_snps_by_dist(geno_df, info_df, min_dist_kb = 0.2)

add_autold_screening(x, r2_threshold = 0.05, ho_threshold = 0.80)
}
\arguments{
  \item{geno_df}{Marker-by-individual genotype or dosage data frame.}
  \item{info_df}{Marker metadata with \code{CHROM} and \code{POS}.}
  \item{min_dist_kb}{Minimum physical distance, in kb, between consecutively retained markers within a chromosome.}
  \item{x}{Genome-scan result containing \code{r2_comp} and normalized tensor columns such as \code{D_norm_11}, \code{D_norm_12}, and so forth.}
  \item{r2_threshold}{Screening threshold for conventional first-order association.}
  \item{ho_threshold}{Screening threshold for the maximum absolute normalized higher-order component.}
}
\value{
\code{thin_snps_by_dist()} returns a list with thinned genotype data, thinned marker metadata, and retained original row indices. \code{add_autold_screening()} returns the input data with \code{H_HO} and \code{ScreenClass} columns.
}
