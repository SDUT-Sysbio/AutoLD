#' Estimate Linkage Disequilibrium for Diploids
#'
#' This function estimates linkage disequilibrium (LD) parameter D and allele frequencies
#' for diploids using an Expectation-Maximization (EM) algorithm.
#'
#' @param dose A numeric vector or factor indicating the 9 possible joint genotype classes for two loci.
#' @return A list containing:
#' \itemize{
#'   \item \code{p}: Estimated allele frequency for locus A.
#'   \item \code{q}: Estimated allele frequency for locus B.
#'   \item \code{D}: Estimated Linkage Disequilibrium parameter.
#'   \item \code{iter}: Number of EM iterations until convergence.
#' }
#' @export
#' @examples
#' \dontrun{
#' true_para <- c(0.6, 0.4, 0.05) # p, q, D
#' sim_data <- sim_2locus_2x_cLD(true_para, n=1000)
#' est_res <- auto2_est(sim_data)
#' print(est_res)
#' }
auto2_est <- function(dose){
  nts <- rep(0, 9)
  nt_table <- table(dose)
  nts[as.numeric(names(nt_table))] <- as.numeric(nt_table)
  N_mat <- matrix(nts, nrow=3, ncol=3)
  nn <- sum(N_mat)

  p <- 0.5; q <- 0.5; D <- 0.0
  loop_k <- 1; diff <- 1

  while(diff > 1e-7 && loop_k < 200) {
    p_old <- p; q_old <- q; D_old <- D
    p11 <- p*q + D; p10 <- p*(1-q) - D
    p01 <- (1-p)*q - D; p00 <- (1-p)*(1-q) + D

    denom <- (p11 * p00 + p10 * p01)
    theta <- ifelse(denom == 0, 0, (p11 * p00) / denom)

    count_AB <- 2 * N_mat[3, 3] + N_mat[3, 2] + N_mat[2, 3] + theta * N_mat[2, 2]
    p11 <- count_AB / (2 * nn)

    count_Ab <- 2 * N_mat[3, 1] + N_mat[3, 2] + N_mat[2, 1] + (1 - theta) * N_mat[2, 2]
    p10 <- count_Ab / (2 * nn)

    count_aB <- 2 * N_mat[1, 3] + N_mat[1, 2] + N_mat[2, 3] + (1 - theta) * N_mat[2, 2]
    p01 <- count_aB / (2 * nn)

    p00 <- 1 - p11 - p10 - p01
    p <- max(min(p11 + p10, 0.9999), 0.0001)
    q <- max(min(p11 + p01, 0.9999), 0.0001)
    D <- p11 * p00 - p10 * p01

    diff <- max(abs(p - p_old), abs(q - q_old), abs(D - D_old))
    loop_k <- loop_k + 1
  }
  return(list(p=p, q=q, D=D, iter=loop_k))
}

#' Estimate Linkage Disequilibrium for Autotetraploids
#'
#' This function estimates LD parameters and double reduction rates for autotetraploids.
#' It supports both full genotype resolution (25 classes) and partial resolution (9 classes).
#'
#' @param geno A data frame or list containing a \code{class} column representing the joint genotype dosage.
#' @param type A character string specifying the data resolution. Use \code{"full"} for 25-class data or \code{"partial"} for 9-class data. Default is \code{"full"}.
#' @return A named numeric vector containing estimated parameters:
#' allele frequencies (\code{pA}, \code{pB}), double reduction rates (\code{aA}, \code{aB}),
#' and LD components (\code{Deab}, \code{DAb}, \code{DaB}, \code{DAB}).
#' @export
#' @examples
#' \dontrun{
#' para <- c(pA=0.6, pB=0.4, aA=0.1, aB=0.14, Deab=0.05, DAb=0.02, DaB=0.01, DAB=0.002)
#' gen4 <- sim_2locus_4x_cLD(para=para, n=500, type="full")
#' res <- auto4_est(geno=gen4, type="full")
#' print(res)
#' }
auto4_est <- function(geno, type="full"){
  dose <- geno$class
  if(type=="full"){
    dp <- auto4_est_LD(dose=dose)
    a1 <- solve_a_from_QAA_4x(dp[1],dp[7])
    a2 <- solve_a_from_QAA_4x(dp[2],dp[8])
  } else {
    dp <- auto4_est_LDp(dose=dose)
    a1 <- solve_a_from_QAA_4x(dp[1],dp[7])
    a2 <- solve_a_from_QAA_4x(dp[2],dp[8])
  }
  init.par <- c(dp[1:2], a1, a2, dp[3:6])
  names(init.par) <- c("pA","pB","aA","aB","Deab","DAb","DaB","DAB")
  return(init.par)
}

#' Estimate Linkage Disequilibrium for Autohexaploids
#'
#' This function estimates LD parameters and double reduction rates for autohexaploids.
#' It handles the high complexity of hexaploid genetics, supporting both full (49 classes)
#' and partial (9 classes) genotype data.
#'
#' @param geno A data frame or list containing a \code{class} column representing the joint genotype dosage.
#' @param type A character string specifying the data resolution. Use \code{"full"} for 49-class data or \code{"partial"} for 9-class data. Default is \code{"full"}.
#' @return A named numeric vector containing estimated parameters:
#' allele frequencies (\code{pA}, \code{pB}), double reduction rates (\code{aA}, \code{aB}),
#' and higher-order LD components (\code{Deab}, \code{DAb}, \code{DaB}, \code{DAAb}, \code{DaBB}, \code{DAB}, \code{DAAB}, \code{DABB}, \code{DAABB}).
#' @export
auto6_est <- function(geno, type="full"){
  dose <- geno$class
  if(type=="full"){
    dp <- auto6_est_LD(dose=dose)
  } else {
    dp <- auto6_est_LDp(dose=dose)
  }
  a1 <- estimate_a_from_gamete6x_freq(p=dp$estp[1], dp$gamA)$a
  a2 <- estimate_a_from_gamete6x_freq(p=dp$estp[2], dp$gamB)$a

  init.par <- c(dp$estp[1:2], a1, a2, dp$estp[c(-c(1:4,6:7))])
  names(init.par) <- c("pA","pB","aA","aB","Deab","DAb","DaB","DAAb","DaBB","DAB","DAAB","DABB","DAABB")
  return(init.par)
}
