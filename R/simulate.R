#' Simulate Two-Locus Genotype Data for Diploids
#'
#' Generates simulated joint genotype classes for two loci in a diploid population
#' based on specified allele frequencies and linkage disequilibrium (LD).
#'
#' @param para A numeric vector of length 3 containing the true parameters in the following order:
#' \code{p} (allele frequency of locus A), \code{q} (allele frequency of locus B), and \code{D} (LD parameter).
#' @param n An integer specifying the sample size. Default is 500.
#' @return A numeric vector representing the simulated genotype classes (ranging from 1 to 9).
#' @export
#' @examples
#' \dontrun{
#' # Define parameters: p = 0.6, q = 0.4, D = 0.05
#' true_para <- c(0.6, 0.4, 0.05)
#' sim_data <- sim_2locus_2x_cLD(para = true_para, n = 1000)
#' table(sim_data)
#' }
sim_2locus_2x_cLD <- function(para, n=500){
  p <- para[1]; q <- para[2]; D <- para[3]
  p11 <- p*q + D; p10 <- p*(1-q) - D
  p01 <- (1-p)*q - D; p00 <- (1-p)*(1-q) + D
  if(any(c(p11, p10, p01, p00) < 0)) stop("Invalid D given allele frequencies")

  M <- matrix(0, 3, 3)
  M[1, 1] <- p00^2; M[3, 3] <- p11^2; M[1, 3] <- p01^2; M[3, 1] <- p10^2
  M[2, 1] <- 2 * p10 * p00; M[2, 3] <- 2 * p11 * p01
  M[1, 2] <- 2 * p01 * p00; M[3, 2] <- 2 * p11 * p10
  M[2, 2] <- 2 * (p11 * p00 + p10 * p01)

  prob_vec <- as.numeric(M)
  cls <- sample(1:9, size=n, replace=TRUE, prob=prob_vec)
  return(cls)
}

#' Simulate Two-Locus Genotype Data for Autotetraploids
#'
#' Generates simulated joint genotype data for two loci in autotetraploids. This function
#' explicitly incorporates allele frequencies, double reduction rates, and various
#' linkage disequilibrium (LD) components.
#'
#' @param para A numeric vector containing true parameters in the following exact order:
#' \code{pA, pB, aA, aB, Deab, DAb, DaB, DAB}. Here \code{aA} and \code{aB} denote the
#' double reduction rates for locus A and B, respectively.
#' @param n An integer specifying the sample size to simulate. Default is 500.
#' @param type A character string specifying the output resolution format: \code{"full"} (25 distinct classes)
#' or \code{"partial"} (collapsed into 9 classes). Default is \code{"full"}.
#' @param hwdA Numeric. Hardy-Weinberg Disequilibrium (HWD) parameter for locus A. Default is 0.
#' @param hwdB Numeric. Hardy-Weinberg Disequilibrium (HWD) parameter for locus B. Default is 0.
#' @return A data frame containing a single column \code{class} that records the simulated
#' joint genotype categories for each individual.
#' @export
#' @examples
#' \dontrun{
#' # Define parameters for an autotetraploid population
#' my_para <- c(pA=0.6, pB=0.4, aA=0.1, aB=0.1, Deab=0.05, DAb=0.01, DaB=0.01, DAB=0.005)
#'
#' # Generate 1000 samples with full 25-class resolution
#' sim_data_full <- sim_2locus_4x_cLD(para = my_para, n = 1000, type = "full")
#' head(sim_data_full)
#' }
sim_2locus_4x_cLD <- function(para, n=500, type="full", hwdA=0, hwdB=0){
  zyp <- DF4_cLD(para=para, hwdA=hwdA, hwdB=hwdB)
  cls <- sample(1:25, size=n, replace=TRUE, prob=zyp)
  if(type=="full"){
    res <- data.frame(class=cls)
  }else{
    gene_c <- cls25_to_9class4(cls=cls)
    res <- data.frame(class=gene_c[,5])
  }
  return(res)
}

cls25_to_9class4 <- function(cls){
  cls <- as.integer(cls)
  if(any(is.na(cls)) || any(cls < 1 | cls > 25)){
    stop("cls must be integers in [1, 25].")
  }
  doseA <- (cls - 1) %/% 5
  doseB <- (cls - 1) %% 5
  collapse3 <- function(d){
    out <- d
    out[d %in% c(1,2,3)] <- 1
    out[d == 4] <- 2
    out
  }
  A3 <- collapse3(doseA)
  B3 <- collapse3(doseB)
  class9 <- 3*A3 + B3 + 1
  cbind(doseA = doseA, doseB = doseB, A3 = A3, B3 = B3, class9 = class9)
}

DF4_cLD <- function(para, hwdA, hwdB){
  gfA <- gamete_margin_DR_4x(p = para[1], a = para[3], D = hwdA)
  gfB <- gamete_margin_DR_4x(p = para[2], a = para[4], D = hwdB)
  gem <- DF4(gfA=gfA, gfB=gfB, pall=para[-c(1:4)])
  gp0 <- auto_DF4_from_gp(gp=gem)
  gp0
}

auto_DF4_from_gp <- function(gp){
  PAABB <- gp[1]; PAABb <- gp[2]; PAAbb <- gp[3]
  PAaBB <- gp[4]; PAaBb <- gp[5]; PAabb <- gp[6]
  PaaBB <- gp[7]; PaaBb <- gp[8]; Paabb <- gp[9]

  zy1 <- (PAABB)^2; zy2 <- 2*PAABB*PAABb; zy3 <- 2*PAABB*PAAbb+(PAABb)^2
  zy4 <- 2*PAABb*PAAbb; zy5 <- (PAAbb)^2; zy6 <- 2*PAABB*PAaBB
  zy7 <- 2*PAABB*PAaBb+2*PAaBB*PAABb
  zy8 <- 2*PAABb*PAaBb+2*PAaBB*PAAbb+2*PAabb*PAABB
  zy9 <- 2*PAAbb*PAaBb+2*PAabb*PAABb; zy10 <- 2*PAabb*PAAbb
  zy11 <- (PAaBB)^2+2*PAABB*PaaBB
  zy12 <- 2*PAaBB*PAaBb+2*PAABB*PaaBb+2*PAABb*PaaBB
  zy13 <- (PAaBb)^2+2*PAABb*PaaBb+2*PAABB*Paabb+2*PAAbb*PaaBB+2*PAabb*PAaBB
  zy14 <- 2*PAabb*PAaBb+2*PAABb*Paabb+2*PAAbb*PaaBb
  zy15 <- (PAabb)^2+2*PAAbb*Paabb; zy16 <- 2*PAaBB*PaaBB
  zy17 <- 2*PaaBB*PAaBb+2*PaaBb*PAaBB
  zy18 <- 2*PaaBb*PAaBb+2*Paabb*PAaBB+2*PaaBB*PAabb
  zy19 <- 2*Paabb*PAaBb+2*PaaBb*PAabb; zy20 <- 2*Paabb*PAabb
  zy21 <- (PaaBB)^2; zy22 <- 2*PaaBB*PaaBb
  zy23 <- (PaaBb)^2+2*pmax(PaaBB,0)*pmax(Paabb,0)
  zy24 <- 2*PaaBb*Paabb; zy25 <- (Paabb)^2

  tmp <- c(zy1,zy2,zy3,zy4,zy5,zy6,zy7,zy8,zy9,zy10,zy11,zy12,zy13,zy14,zy15,zy16,zy17,zy18,zy19,zy20,zy21,zy22,zy23,zy24,zy25)
  names(tmp) <- NULL
  tmp
}

DF4 <- function(gfA, gfB, pall){
  pA <- gfA[1] + gfA[2]/2; pa <- 1 - pA
  pB <- gfB[1] + gfB[2]/2; pb <- 1 - pB
  DA    <- gfA[1] - pA^2; DB    <- gfB[1] - pB^2
  Deab  <- pall[1]; DAb   <- pall[2]; DaB   <- pall[3]; DAB   <- pall[4]

  PAABB <- pA^2*pB^2+pA^2*DB+pB^2*DA+2*pA*pB*Deab+2*pB*DAb+2*pA*DaB+DA*DB+Deab^2+DAB
  PAABb <- 2*(pA^2*pB*pb - pA^2*DB + pB*pb*DA + (pA*pb - pA*pB)*Deab + (pb-pB)*DAb - 2*pA*DaB) - 2*(DA*DB+Deab^2+DAB)
  PAAbb <- pA^2*pb^2 + pA^2*DB + pb^2*DA - 2*pA*pb*Deab - 2*pb*DAb + 2*pA*DaB + DA*DB + Deab^2 + DAB
  PAaBB <- 2*(pA*pa*pB^2 + pA*pa*DB - pB^2*DA + (pa*pB-pA*pB)*Deab - 2*pB*DAb + (pa-pA)*DaB) - 2*(DA*DB+Deab^2+DAB)
  PAaBb <- 2*(pA*pa*pB*pb - pA*pa*DB - pB*pb*DA + (pA*pB+pa*pb)*Deab - (pA*pb+pa*pB)*Deab + (pB-pb)*DAb + (pA-pa)*DaB) +
    2*(pA*pa*pB*pb - pA*pa*DB - pB*pb*DA + (pB-pb)*DAb + (pA-pa)*DaB) + 4*(DA*DB+Deab^2+DAB)
  PAabb <- 2*(pA*pa*pb^2 + pA*pa*DB - pb^2*DA + (pA*pb-pa*pb)*Deab + 2*pb*DAb + (pa-pA)*DaB) - 2*(DA*DB+Deab^2+DAB)
  PaaBB <- pa^2*pB^2 + pa^2*DB + pB^2*DA - 2*pa*pB*Deab + 2*pB*DAb - 2*pa*DaB + DA*DB + Deab^2 + DAB
  PaaBb <- 2*(pa^2*pB*pb - pa^2*DB + pB*pb*DA + (pa*pB-pa*pb)*Deab + (pb-pB)*DAb + 2*pa*DaB) - 2*(DA*DB+Deab^2+DAB)
  Paabb <- pa^2*pb^2 + pb^2*DA + pa^2*DB + 2*pa*pb*Deab - 2*pb*DAb - 2*pa*DaB + DA*DB + Deab^2 + DAB

  tmp <- c(PAABB,PAABb,PAAbb,PAaBB,PAaBb,PAabb,PaaBB,PaaBb,Paabb)
  names(tmp) <- c("AA_BB","AA_Bb","AA_bb","Aa_BB","Aa_Bb","Aa_bb","aa_BB","aa_Bb","aa_bb")
  tmp
}

#' Simulate Unified Haplotype Data for Autotetraploids
#'
#' Generates joint genotype classes for autotetraploids using a unified haplotype-based approach.
#' This is particularly useful for theoretical validations and comparisons.
#'
#' @param pall A numeric vector of length 3 containing the true parameters:
#' \code{pA} (allele frequency of A), \code{pB} (allele frequency of B), and \code{D} (overall LD).
#' @param n An integer specifying the sample size. Default is 1000.
#' @param mode A character string specifying the output format: \code{"full"} (25 classes)
#' or \code{"partial"} (9 classes). Default is \code{"partial"}.
#' @return A numeric vector representing the simulated genotype classes for each individual.
#' @export
#' @examples
#' \dontrun{
#' hap_para <- c(0.6, 0.4, 0.12) # pA, pB, D
#' hap_data <- sim_4x_unified_hap(pall = hap_para, n = 1000, mode = "partial")
#' table(hap_data)
#' }
sim_4x_unified_hap <- function(pall, n=1000, mode="partial") {
  gametes <- std_4x_gametes_hap(pall)
  if(mode == "full") {
    prob_vec <- std_4x_zygotes_25_hap(gametes)
    return(sample(1:25, n, replace=TRUE, prob=prob_vec))
  } else if(mode == "partial") {
    prob_vec <- std_4x_zygotes_9_hap(gametes)
    return(sample(1:9, n, replace=TRUE, prob=prob_vec))
  } else {
    stop("Invalid mode. Please use 'full' or 'partial'.")
  }
}

std_4x_gametes_hap <- function(pall){
  pA <- pall[1]; pB <- pall[2]; D <- pall[3]
  pa <- 1-pA; pb <- 1-pB
  h_AB <- pA*pB + D; h_Ab <- pA*pb - D
  h_aB <- pa*pB - D; h_ab <- pa*pb + D
  if(any(c(h_AB, h_Ab, h_aB, h_ab) < 0)) stop("Invalid D for allele frequencies")

  g <- numeric(9)
  g[1] <- h_AB^2; g[2] <- 2 * h_AB * h_Ab; g[3] <- h_Ab^2
  g[4] <- 2 * h_AB * h_aB; g[5] <- 2 * h_AB * h_ab + 2 * h_Ab * h_aB
  g[6] <- 2 * h_Ab * h_ab; g[7] <- h_aB^2
  g[8] <- 2 * h_aB * h_ab; g[9] <- h_ab^2
  return(g)
}

std_4x_zygotes_25_hap <- function(g){
  z <- numeric(25)
  z[1] <- g[1]^2; z[2] <- 2*g[1]*g[2]; z[3] <- 2*g[1]*g[3] + g[2]^2
  z[4] <- 2*g[2]*g[3]; z[5] <- g[3]^2; z[6] <- 2*g[1]*g[4]
  z[7] <- 2*g[1]*g[5] + 2*g[4]*g[2]; z[8] <- 2*g[2]*g[5] + 2*g[4]*g[3] + 2*g[6]*g[1]
  z[9] <- 2*g[3]*g[5] + 2*g[6]*g[2]; z[10] <- 2*g[6]*g[3]
  z[11] <- g[4]^2 + 2*g[1]*g[7]; z[12] <- 2*g[4]*g[5] + 2*g[1]*g[8] + 2*g[2]*g[7]
  z[13] <- g[5]^2 + 2*g[2]*g[8] + 2*g[1]*g[9] + 2*g[3]*g[7] + 2*g[6]*g[4]
  z[14] <- 2*g[6]*g[5] + 2*g[2]*g[9] + 2*g[3]*g[8]; z[15] <- g[6]^2 + 2*g[3]*g[9]
  z[16] <- 2*g[4]*g[7]; z[17] <- 2*g[7]*g[5] + 2*g[8]*g[4]
  z[18] <- 2*g[8]*g[5] + 2*g[9]*g[4] + 2*g[7]*g[6]; z[19] <- 2*g[9]*g[5] + 2*g[8]*g[6]
  z[20] <- 2*g[9]*g[6]; z[21] <- g[7]^2; z[22] <- 2*g[7]*g[8]
  z[23] <- g[8]^2 + 2*g[7]*g[9]; z[24] <- 2*g[8]*g[9]; z[25] <- g[9]^2
  return(z)
}

std_4x_zygotes_9_hap <- function(g){
  z <- numeric(9)
  z[1] <- g[1]^2
  z[2] <- 2*g[1]*g[2] + 2*g[1]*g[3] + g[2]^2 + 2*g[2]*g[3]; z[3] <- g[3]^2
  z[4] <- 2*g[1]*g[4] + g[4]^2 + 2*g[1]*g[7] + 2*g[4]*g[7]
  z[5] <- 2*g[1]*g[5] + 2*g[4]*g[2] + 2*g[2]*g[5] + 2*g[4]*g[3] + 2*g[6]*g[1] + 2*g[3]*g[5] + 2*g[6]*g[2] + 2*g[4]*g[5] + 2*g[1]*g[8] + 2*g[2]*g[7] + g[5]^2 + 2*g[2]*g[8] + 2*g[1]*g[9] + 2*g[3]*g[7] + 2*g[6]*g[4] + 2*g[6]*g[5] + 2*g[2]*g[9] + 2*g[3]*g[8] + 2*g[7]*g[5] + 2*g[8]*g[4] + 2*g[8]*g[5] + 2*g[9]*g[4] + 2*g[7]*g[6] + 2*g[9]*g[5] + 2*g[8]*g[6]
  z[6] <- 2*g[6]*g[3] + g[6]^2 + 2*g[3]*g[9] + 2*g[9]*g[6]; z[7] <- g[7]^2
  z[8] <- 2*g[7]*g[8] + g[8]^2 + 2*g[7]*g[9] + 2*g[8]*g[9]; z[9] <- g[9]^2
  return(z)
}

#' Simulate Unified Haplotype Data for Autohexaploids
#'
#' Generates simulated joint genotype data for autohexaploids using a unified haplotype-based
#' probability model. It handles the high dimensionality of hexaploid genetics by generating
#' data in either 49-class or 9-class resolutions.
#'
#' @param pall A numeric vector of length 3 containing the true parameters:
#' \code{pA} (allele frequency of A), \code{pB} (allele frequency of B), and \code{D} (overall LD).
#' @param n An integer specifying the sample size. Default is 1000.
#' @param mode A character string specifying the output format: \code{"full"} (49 distinct classes)
#' or \code{"partial"} (collapsed into 9 classes). Default is \code{"partial"}.
#' @return A numeric vector representing the simulated joint genotype classes.
#' @export
#' @examples
#' \dontrun{
#' hex_para <- c(0.6, 0.4, 0.05) # pA, pB, D
#' hex_data <- sim_6x_unified_hap(pall = hex_para, n = 1000, mode = "full")
#' table(hex_data)
#' }
sim_6x_unified_hap <- function(pall, n=1000, mode="partial") {
  gametes <- std_6x_gametes_hap(pall)
  prob_49 <- std_6x_zygotes_49_hap(gametes)
  if(mode == "full") {
    return(sample(1:49, n, replace=TRUE, prob=prob_49))
  } else if(mode == "partial") {
    prob_9 <- std_6x_collapse_to_9_hap(prob_49)
    return(sample(1:9, n, replace=TRUE, prob=prob_9))
  } else {
    stop("Invalid mode. Use 'full' or 'partial'.")
  }
}

std_6x_gametes_hap <- function(pall){
  pA <- pall[1]; pa <- 1-pA; pB <- pall[2]; pb <- 1-pB; D <- pall[3]
  hAB <- pA*pB+D; hAb <- pA*pb-D; haB <- pa*pB-D; hab <- pa*pb+D
  g <- numeric(16)
  g[1] <- hAB^3; g[2] <- 3*hAB^2*hAb; g[3] <- 3*hAB*hAb^2; g[4] <- hAb^3
  g[5] <- 3*hAB^2*haB; g[6] <- 3*hAB^2*hab + 6*hAB*hAb*haB
  g[7] <- 3*hAb^2*haB + 6*hAB*hAb*hab; g[8] <- 3*hAb^2*hab
  g[9] <- 3*hAB*haB^2; g[10] <- 3*hAb*haB^2 + 6*hAB*haB*hab
  g[11] <- 3*hAB*hab^2 + 6*hAb*haB*hab; g[12] <- 3*hAb*hab^2
  g[13] <- haB^3; g[14] <- 3*haB^2*hab; g[15] <- 3*haB*hab^2; g[16] <- hab^3
  return(g)
}

std_6x_zygotes_49_hap <- function(gp){
  PAAABBB=gp[1]; PAAABBb=gp[2]; PAAABbb=gp[3]; PAAAbbb=gp[4]
  PAAaBBB=gp[5]; PAAaBBb=gp[6]; PAAaBbb=gp[7]; PAAabbb=gp[8]
  PAaaBBB=gp[9]; PAaaBBb=gp[10]; PAaaBbb=gp[11]; PAaabbb=gp[12]
  PaaaBBB=gp[13]; PaaaBBb=gp[14]; PaaaBbb=gp[15]; Paaabbb=gp[16]
  zy <- numeric(49)
  zy[1] <- (PAAABBB)^2; zy[2] <- 2*PAAABBB*PAAABBb; zy[3] <- 2*PAAABBB*PAAABbb+(PAAABBb)^2
  zy[4] <- 2*PAAABBB*PAAAbbb+2*PAAABBb*PAAABbb; zy[5] <- 2*PAAABBb*PAAAbbb+(PAAABbb)^2
  zy[6] <- 2*PAAABbb*PAAAbbb; zy[7] <- (PAAAbbb)^2; zy[8] <- 2*PAAABBB*PAAaBBB
  zy[9] <- 2*PAAABBB*PAAaBBb+2*PAAaBBB*PAAABBb; zy[10] <- 2*PAAABBB*PAAaBbb+2*PAAaBBB*PAAABbb+2*PAAABBb*PAAaBBb
  zy[11] <- 2*PAAABBB*PAAabbb+2*PAAaBBB*PAAAbbb+2*PAAABBb*PAAaBbb+2*PAAaBBb*PAAABbb
  zy[12] <- 2*PAAABBb*PAAabbb+2*PAAaBBb*PAAAbbb+2*PAAABbb*PAAaBbb
  zy[13] <- 2*PAAABbb*PAAabbb+2*PAAaBbb*PAAAbbb; zy[14] <- 2*PAAAbbb*PAAabbb
  zy[15] <- 2*PAAABBB*PAaaBBB+(PAAaBBB)^2; zy[16] <- 2*PAAABBB*PAaaBBb+2*PAAaBBB*PAAaBBb+2*PAaaBBB*PAAABBb
  zy[17] <- (PAAaBBb)^2+2*PAAABBB*PAaaBbb+2*PAAaBBB*PAAaBbb+2*PAaaBBB*PAAABbb+2*PAAABBb*PAaaBBb
  zy[18] <- 2*PAAABBB*PAaabbb+2*PAAaBBB*PAAabbb+2*PAaaBBB*PAAAbbb+2*PAAABBb*PAaaBbb+2*PAAaBBb*PAAaBbb+2*PAaaBBb*PAAABbb
  zy[19] <- (PAAaBbb)^2+2*PAAABBb*PAaabbb+2*PAAaBBb*PAAabbb+2*PAaaBBb*PAAAbbb+2*PAAABbb*PAaaBbb
  zy[20] <- 2*PAAABbb*PAaabbb+2*PAaaBbb*PAAAbbb+2*PAAaBbb*PAAabbb; zy[21] <- (PAAabbb)^2+2*PAAAbbb*PAaabbb
  zy[22] <- 2*PAAABBB*PaaaBBB+2*PAAaBBB*PAaaBBB; zy[23] <- 2*PAAABBB*PaaaBBb+2*PAAaBBB*PAaaBBb+2*PAaaBBB*PAAaBBb+2*PaaaBBB*PAAABBb
  zy[24] <- 2*PAAABBB*PaaaBbb+2*PAAaBBB*PAaaBbb+2*PAaaBBB*PAAaBbb+2*PaaaBBB*PAAABbb+2*PAAABBb*PaaaBBb+2*PAAaBBb*PAaaBBb
  zy[25] <- 2*PAAABBB*Paaabbb+2*PAAaBBB*PAaabbb+2*PAaaBBB*PAAabbb+2*PaaaBBB*PAAAbbb+2*PAAABBb*PaaaBbb+2*PAAaBBb*PAaaBbb+2*PAaaBBb*PAAaBbb+2*PaaaBBb*PAAABbb
  zy[26] <- 2*PAAABBb*Paaabbb+2*PAAaBBb*PAaabbb+2*PAaaBBb*PAAabbb+2*PaaaBBb*PAAAbbb+2*PAAABbb*PaaaBbb+2*PAAaBbb*PAaaBbb
  zy[27] <- 2*PAAABbb*Paaabbb+2*PAAaBbb*PAaabbb+2*PAaaBbb*PAAabbb+2*PaaaBbb*PAAAbbb; zy[28] <- 2*PAAAbbb*Paaabbb+2*PAAabbb*PAaabbb
  zy[29] <- (PAaaBBB)^2+2*PAAaBBB*PaaaBBB; zy[30] <- 2*PAAaBBB*PaaaBBb+2*PAaaBBB*PAaaBBb+2*PaaaBBB*PAAaBBb
  zy[31] <- (PAaaBBb)^2+2*PAAaBBB*PaaaBbb+2*PAaaBBB*PAaaBbb+2*PaaaBBB*PAAaBbb+2*PAAaBBb*PaaaBBb
  zy[32] <- 2*PAAaBBB*Paaabbb+2*PAaaBBB*PAaabbb+2*PaaaBBB*PAAabbb+2*PAAaBBb*PaaaBbb+2*PAaaBBb*PAaaBbb+2*PaaaBBb*PAAaBbb
  zy[33] <- (PAaaBbb)^2+2*PAAaBBb*Paaabbb+2*PAaaBBb*PAaabbb+2*PaaaBBb*PAAabbb+2*PAAaBbb*PaaaBbb
  zy[34] <- 2*PAAaBbb*Paaabbb+2*PAaaBbb*PAaabbb+2*PaaaBbb*PAAabbb; zy[35] <- (PAaabbb)^2+2*PAAabbb*Paaabbb
  zy[36] <- 2*PAaaBBB*PaaaBBB; zy[37] <- 2*PAaaBBB*PaaaBBb+2*PaaaBBB*PAaaBBb
  zy[38] <- 2*PAaaBBB*PaaaBbb+2*PaaaBBB*PAaaBbb+2*PAaaBBb*PaaaBBb; zy[39] <- 2*PAaaBBB*Paaabbb+2*PaaaBBB*PAaabbb+2*PAaaBBb*PaaaBbb+2*PaaaBBb*PAaaBbb
  zy[40] <- 2*PAaaBBb*Paaabbb+2*PaaaBBb*PAaabbb+2*PAaaBbb*PaaaBbb; zy[41] <- 2*PAaaBbb*Paaabbb+2*PaaaBbb*PAaabbb
  zy[42] <- 2*PAaabbb*Paaabbb; zy[43] <- (PaaaBBB)^2; zy[44] <- 2*PaaaBBB*PaaaBBb
  zy[45] <- (PaaaBBb)^2+2*PaaaBBB*PaaaBbb; zy[46] <- 2*PaaaBBB*Paaabbb+2*PaaaBBb*PaaaBbb
  zy[47] <- (PaaaBbb)^2+2*PaaaBBb*Paaabbb; zy[48] <- 2*PaaaBbb*Paaabbb; zy[49] <- (Paaabbb)^2
  return(zy)
}

std_6x_collapse_to_9_hap <- function(prob_49){
  prob_mat <- matrix(prob_49, nrow=7, ncol=7)
  r1 <- prob_mat[1,]; r2 <- colSums(prob_mat[2:6,]); r3 <- prob_mat[7,]
  mat_3x7 <- rbind(r1, r2, r3)
  c1 <- mat_3x7[,1]; c2 <- rowSums(mat_3x7[,2:6]); c3 <- mat_3x7[,7]
  mat_3x3 <- cbind(c1, c2, c3)
  return(as.numeric(mat_3x3))
}


#' Simulate Two-Locus Genotype Data for Autohexaploids
#'
#' Generates simulated joint genotype data for two loci in autohexaploids. This function
#' explicitly incorporates allele frequencies, double reduction rates, and the complex,
#' higher-order linkage disequilibrium (LD) components unique to hexaploids.
#'
#' @param para A numeric vector containing true parameters in the exact order:
#' \code{pA, pB, aA, aB}, followed by the comprehensive LD components
#' (\code{Deab, DAb, DaB, DAAb, DaBB, DAB, DAAB, DABB, DAABB}). Here \code{aA} and \code{aB}
#' denote the double reduction rates for locus A and B, respectively.
#' @param n An integer specifying the sample size to simulate. Default is 500.
#' @param type A character string specifying the output resolution format: \code{"full"}
#' (49 distinct classes) or \code{"partial"} (collapsed into 9 classes). Default is \code{"full"}.
#' @param hwdA Numeric. Hardy-Weinberg Disequilibrium (HWD) parameter for locus A. Default is 0.
#' @param hwdB Numeric. Hardy-Weinberg Disequilibrium (HWD) parameter for locus B. Default is 0.
#' @return A data frame containing a single column \code{class} that records the simulated
#' joint genotype categories for each individual.
#' @export
#' @examples
#' \dontrun{
#' # Define parameters for an autohexaploid population
#' # Order: pA, pB, aA, aB, and the 9 LD components
#' my_para6x <- c(pA=0.6, pB=0.4, aA=0.05, aB=0.08, Deab=0.05, DAb=0.01, DaB=0.01,
#'                DAAb=0.002, DaBB=0.002, DAB=0.005, DAAB=0.001, DABB=0.001, DAABB=0.0005)
#'
#' # Generate 1000 samples with full 49-class resolution
#' sim_data_6x <- sim_2locus_6x_cLD(para = my_para6x, n = 1000, type = "full")
#' head(sim_data_6x)
#' }
sim_2locus_6x_cLD <- function(para, n=500,type="full",hwdA=0,hwdB=0){

  prob_vec <- DF6_cLD(para=para,hwdA=hwdA,hwdB=hwdB)
  if (any(prob_vec < 0)) {
    cat(sprintf("\n[Warning] Negative probability detected! Min Value: %.6f\n", min(prob_vec)))
    prob_vec[prob_vec < 0] <- 0
  }

  sum_p <- sum(prob_vec)
  if (sum_p == 0) {
    stop("Error: All probabilities are 0 (or negative turned to 0). Parameters conflict heavily!")
  } else {
    prob_vec <- prob_vec/sum_p
  }
  cls <- sample(1:49, size=n, replace=TRUE, prob=prob_vec)
  if(type=="full"){
    res <- data.frame(class=cls)
  }else{
    gene_c <- cls49_to_9class6(cls=cls)
    res <- data.frame(class=gene_c[,5])
  }
  return(res)
}


cls49_to_9class6 <- function(cls){
  cls <- as.integer(cls)
  if(any(is.na(cls)) || any(cls < 1 | cls > 49)){
    stop("cls must be integers in [1, 49] for hexaploid.")
  }
  doseA <- (cls - 1) %/% 7
  doseB <- (cls - 1) %% 7
  collapse3_6x <- function(d){
    out <- d
    out[d %in% 1:5] <- 1
    out[d == 6] <- 2
    out
  }

  A3 <- collapse3_6x(doseA)
  B3 <- collapse3_6x(doseB)

  # 3) 3 x 3 -> 1..9
  class9 <- 3*A3 + B3 + 1

  cbind(
    doseA = doseA,
    doseB = doseB,
    A3 = A3,
    B3 = B3,
    class9 = class9
  )
}


DF6_cLD <- function(para,hwdA,hwdB){

  gfA <- gamete_margin_DR_6xD(p = para[1], a = para[3], D = hwdA)  # A位点配子边际
  gfB <- gamete_margin_DR_6xD(p = para[2], a = para[4], D = hwdB)  # B位点配子边际

  names(gfA) <- names(gfB) <- c("aaa","Aaa","AAa","AAA")
  Da <- DA2_DA3_from_gamete6x(q=gfA)
  Db <- DA2_DA3_from_gamete6x(q=gfB)

  gem <- DF6(pall=c(para[1],para[2],Da[2],Db[2],para[5],Da[3],Db[3],para[-c(1:5)]))
  gp0 <- auto_DF6_from_gp(gp=gem)
  gp0
}

DF6 <- function(pall){


  pA <- pall[1]
  pB <- pall[2]

  HM <- H(pA=pA,pB=pB)
  nH <- dim(HM)[1]

  DM <- c(1,pall[-c(1:2)])

  PMM <- HM*matrix(rep(DM,nH),nrow=nH,byrow=T)

  tmp <- rowSums(PMM)

  return(tmp)
}



DA2_DA3_from_gamete6x <- function(q){
  # q: named vector with names c("aaa","Aaa","AAa","AAA")
  q0 <- unname(q["aaa"])
  q1 <- unname(q["Aaa"])
  q2 <- unname(q["AAa"])
  q3 <- unname(q["AAA"])

  # first moments
  p  <- (q1 + 2*q2 + 3*q3) / 3
  p2 <- (q2 + 3*q3) / 3
  p3 <- q3

  DA2 <- p2 - p^2
  DA3 <- p3 - 3*p*p2 + 2*p^3

  return(c(p = p, DA2 = DA2, DA3 = DA3))
}


gamete_margin_DR_6x <- function(p, a){

  p <- as.numeric(p)
  a <- as.numeric(a)

  stopifnot(p >= 0, p <= 1, a >= 0, a <= 1)

  # 6x zygote genotype distribution: k = number of A copies in the individual (0..6)
  k <- 0:6
  Pk <- dbinom(k, size = 6, prob = p)  # choose(6,k) p^k (1-p)^(6-k)

  # helper: Hypergeometric part (no DR): sample 3 chromosomes from 6 without replacement
  hg_prob <- function(kA){
    # i = A copies in gamete (0..3)
    denom <- choose(6, 3)
    i <- 0:3
    num <- choose(kA, i) * choose(6 - kA, 3 - i)
    # invalid combinations become 0 automatically because choose(n, m)=0 for m>n in R? (actually choose gives 0 for m>n often)
    num / denom
  }

  # helper: DR part (minimal mechanism):
  # pick one homolog twice (2 copies), plus one homolog from remaining 5 (1 copy)
  dr_prob <- function(kA){
    # X = allele on duplicated homolog (A=1 with prob kA/6)
    # Y = allele on single homolog among remaining 5
    # i = 2X + Y, i in {0,1,2,3}
    i <- 0:3
    out <- rep(0, 4)

    # handle edge cases kA=0..6 cleanly
    # P(X=1) = kA/6 ; P(X=0) = (6-kA)/6
    pX1 <- kA / 6
    pX0 <- (6 - kA) / 6

    # If X=0, remaining A count is kA, sample 1 from remaining 5 => P(Y=1|X=0)=kA/5
    # If X=1, remaining A count is kA-1, sample 1 from remaining 5 => P(Y=1|X=1)=(kA-1)/5
    pY1_X0 <- ifelse(5 > 0, kA / 5, 0)
    pY0_X0 <- 1 - pY1_X0

    pY1_X1 <- ifelse(5 > 0, (kA - 1) / 5, 0)
    pY0_X1 <- 1 - pY1_X1

    # i=0: X=0, Y=0
    out[1] <- pX0 * pY0_X0
    # i=1: X=0, Y=1
    out[2] <- pX0 * pY1_X0
    # i=2: X=1, Y=0
    out[3] <- pX1 * pY0_X1
    # i=3: X=1, Y=1
    out[4] <- pX1 * pY1_X1

    # numerical guard
    pmax(out, 0)
  }

  # gamete A-count distribution marginalizing over k
  # qi[i+1] corresponds to i=0..3
  qi <- rep(0, 4)
  for (idx in seq_along(k)){
    kA <- k[idx]
    Pi <- (1 - a) * hg_prob(kA) + a * dr_prob(kA)  # length 4, i=0..3
    qi <- qi + Pk[idx] * Pi
  }

  # map i=3,2,1,0 to AAA, AAa, Aaa, aaa
  q <- c(AAA = qi[4], AAa = qi[3], Aaa = qi[2], aaa = qi[1])
  q <- q / sum(q)
  return(q)
}



gamete_margin_DR_6xD <- function(p, a,D){

  p <- as.numeric(p)
  a <- as.numeric(a)
  D <- as.numeric(D)

  stopifnot(p >= 0, p <= 1, a >= 0, a <= 1)

  # 6x zygote genotype distribution: k = number of A copies in the individual (0..6)
  k <- 0:6
  Pk <- dbinom(k, size = 6, prob = p)  # choose(6,k) p^k (1-p)^(6-k)

  # helper: Hypergeometric part (no DR): sample 3 chromosomes from 6 without replacement
  hg_prob <- function(kA){
    # i = A copies in gamete (0..3)
    denom <- choose(6, 3)
    i <- 0:3
    num <- choose(kA, i) * choose(6 - kA, 3 - i)
    # invalid combinations become 0 automatically because choose(n, m)=0 for m>n in R? (actually choose gives 0 for m>n often)
    num / denom
  }

  # helper: DR part (minimal mechanism):
  # pick one homolog twice (2 copies), plus one homolog from remaining 5 (1 copy)
  dr_prob <- function(kA){
    # X = allele on duplicated homolog (A=1 with prob kA/6)
    # Y = allele on single homolog among remaining 5
    # i = 2X + Y, i in {0,1,2,3}
    i <- 0:3
    out <- rep(0, 4)

    # handle edge cases kA=0..6 cleanly
    # P(X=1) = kA/6 ; P(X=0) = (6-kA)/6
    pX1 <- kA / 6
    pX0 <- (6 - kA) / 6

    # If X=0, remaining A count is kA, sample 1 from remaining 5 => P(Y=1|X=0)=kA/5
    # If X=1, remaining A count is kA-1, sample 1 from remaining 5 => P(Y=1|X=1)=(kA-1)/5
    pY1_X0 <- ifelse(5 > 0, kA / 5, 0)
    pY0_X0 <- 1 - pY1_X0

    pY1_X1 <- ifelse(5 > 0, (kA - 1) / 5, 0)
    pY0_X1 <- 1 - pY1_X1

    # i=0: X=0, Y=0
    out[1] <- pX0 * pY0_X0
    # i=1: X=0, Y=1
    out[2] <- pX0 * pY1_X0
    # i=2: X=1, Y=0
    out[3] <- pX1 * pY0_X1
    # i=3: X=1, Y=1
    out[4] <- pX1 * pY1_X1

    # numerical guard
    pmax(out, 0)
  }

  # gamete A-count distribution marginalizing over k
  # qi[i+1] corresponds to i=0..3
  qi <- rep(0, 4)
  for (idx in seq_along(k)){
    kA <- k[idx]
    Pi <- (1 - a) * hg_prob(kA) + a * dr_prob(kA)  # length 4, i=0..3
    qi <- qi + Pk[idx] * Pi
  }

  # map i=3,2,1,0 to AAA, AAa, Aaa, aaa
  q_base <- c(AAA = qi[4], AAa = qi[3], Aaa = qi[2], aaa = qi[1])

  # normalize safeguard (should already sum to 1)
  q_final <- c(
    aaa = q_base["aaa"] + D,
    Aaa = q_base["Aaa"] - D,
    AAa = q_base["AAa"] - D,
    AAA = q_base["AAA"] + D
  )

  # 边界保护
  q_final[q_final < 0] <- 0
  q_final <- q_final / sum(q_final)
  #q <- q / sum(q)
  return(q_final)
}



H <- function(pA,pB){

  p <- c(pA,1-pA); q <- c(pB,1-pB)

  hi <- c(1,3,3,1,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,1,3,3,1)

  g6i1 <- rbind(c(1,1,1),c(1,1,1),c(1,1,1),c(1,1,1),
                c(1,1,2),c(1,1,2),c(1,1,2),c(1,1,2),c(1,1,2),c(1,1,2),c(1,1,2),c(1,1,2),
                c(1,2,2),c(1,2,2),c(1,2,2),c(1,2,2),c(1,2,2),c(1,2,2),c(1,2,2),c(1,2,2),
                c(2,2,2),c(2,2,2),c(2,2,2),c(2,2,2))

  g6i2 <- rbind(c(1,1,1),c(1,1,2),c(1,2,2),c(2,2,2),
                c(1,1,1),c(1,1,2),c(1,2,1),c(2,1,1),c(1,2,2),c(2,1,2),c(2,2,1),c(2,2,2),
                c(1,1,1),c(1,1,2),c(1,2,1),c(2,1,1),c(1,2,2),c(2,1,2),c(2,2,1),c(2,2,2),
                c(1,1,1),c(1,1,2),c(1,2,2),c(2,2,2))
  si_A2_B2 <- c(1,2,1,3,2,3)
  si_ab2 <- c(1,2,3,3,2,3)
  p1 <- matrix(p[g6i1],ncol=3);q1 <- matrix(q[g6i2],ncol=3)
  pq <- c(p1[,1]*p1[,2]*p1[,3]*q1[,1]*q1[,2]*q1[,3])*hi

  DA2i <- apply(cbind(g6i1,q1),1,comDA2_B2,ii=si_A2_B2,p=p)*hi
  DB2i <- apply(cbind(g6i2,p1),1,comDA2_B2,ii=si_A2_B2,p=q)*hi
  Deabi <- apply(cbind(g6i1,g6i2),1,comDab1,p=p,q=q)*hi
  DA3i <- -(-1)^apply(g6i1,1,sum)*apply(matrix(q[g6i2],ncol=3),1,prod)*hi
  DB3i <- -(-1)^apply(g6i2,1,sum)*apply(matrix(p[g6i1],ncol=3),1,prod)*hi
  DAbi <- apply(cbind(g6i1,g6i2),1,comDAb,p=p,q=q)*hi
  DaBi <- apply(cbind(g6i2,g6i1),1,comDAb,p=q,q=p)*hi
  DAAbi <- apply(cbind(g6i1,g6i2),1,comDAAb,p=q)*hi
  DaBBi <- apply(cbind(g6i2,g6i1),1,comDAAb,p=p)*hi
  DABi <- apply(cbind(g6i1,g6i2),1,comDAB,p=p,q=q)*hi
  DAABi <- apply(cbind(g6i1,g6i2),1,comDAAB,p=q)*hi
  DABBi <- apply(cbind(g6i2,g6i1),1,comDAAB,p=p)*hi
  DAABBi <- (-1)^apply(cbind(g6i1,g6i2),1,sum)*hi

  PM <- cbind(pq,DA2i,DB2i,Deabi,DA3i,DB3i,DAbi,DaBi,DAAbi,DaBBi,DABi,DAABi,DABBi,DAABBi)

  PM1 <- rbind(PM[1:5,],colSums(PM[6:8,]),colSums(PM[9:11,]),PM[12,],PM[13,],colSums(PM[14:16,]),
               colSums(PM[17:19,]),PM[20,],PM[21,],PM[22,],PM[23,],PM[24,])

  return(PM1)
}


auto_DF6_from_gp <- function(gp){

  PAAABBB <- gp[1]
  PAAABBb <- gp[2]
  PAAABbb <- gp[3]
  PAAAbbb <- gp[4]
  PAAaBBB <- gp[5]
  PAAaBBb <- gp[6]
  PAAaBbb <- gp[7]
  PAAabbb <- gp[8]
  PAaaBBB <- gp[9]
  PAaaBBb <- gp[10]
  PAaaBbb <- gp[11]
  PAaabbb <- gp[12]
  PaaaBBB <- gp[13]
  PaaaBBb <- gp[14]
  PaaaBbb <- gp[15]
  Paaabbb <- gp[16]

  zy1 <- (PAAABBB)^2
  zy2 <- 2*PAAABBB*PAAABBb
  zy3 <- 2*PAAABBB*PAAABbb+(PAAABBb)^2
  zy4 <- 2*PAAABBB*PAAAbbb+2*PAAABBb*PAAABbb
  zy5 <- 2*PAAABBb*PAAAbbb+(PAAABbb)^2
  zy6 <- 2*PAAABbb*PAAAbbb
  zy7 <- (PAAAbbb)^2

  zy8 <- 2*PAAABBB*PAAaBBB
  zy9 <- 2*PAAABBB*PAAaBBb+2*PAAaBBB*PAAABBb
  zy10 <- 2*PAAABBB*PAAaBbb+2*PAAaBBB*PAAABbb+2*PAAABBb*PAAaBBb
  zy11 <- 2*PAAABBB*PAAabbb+2*PAAaBBB*PAAAbbb+2*PAAABBb*PAAaBbb+2*PAAaBBb*PAAABbb
  zy12 <- 2*PAAABBb*PAAabbb+2*PAAaBBb*PAAAbbb+2*PAAABbb*PAAaBbb
  zy13 <- 2*PAAABbb*PAAabbb+2*PAAaBbb*PAAAbbb
  zy14 <- 2*PAAAbbb*PAAabbb

  zy15 <- 2*PAAABBB*PAaaBBB+(PAAaBBB)^2
  zy16 <- 2*PAAABBB*PAaaBBb+2*PAAaBBB*PAAaBBb+2*PAaaBBB*PAAABBb
  zy17 <- (PAAaBBb)^2+2*PAAABBB*PAaaBbb+2*PAAaBBB*PAAaBbb+2*PAaaBBB*PAAABbb+2*PAAABBb*PAaaBBb
  zy18 <- 2*PAAABBB*PAaabbb+2*PAAaBBB*PAAabbb+2*PAaaBBB*PAAAbbb+2*PAAABBb*PAaaBbb+2*PAAaBBb*PAAaBbb+2*PAaaBBb*PAAABbb
  zy19 <- (PAAaBbb)^2+2*PAAABBb*PAaabbb+2*PAAaBBb*PAAabbb+2*PAaaBBb*PAAAbbb+2*PAAABbb*PAaaBbb
  zy20 <- 2*PAAABbb*PAaabbb+2*PAaaBbb*PAAAbbb+2*PAAaBbb*PAAabbb
  zy21 <- (PAAabbb)^2+2*PAAAbbb*PAaabbb

  zy22 <- 2*PAAABBB*PaaaBBB+2*PAAaBBB*PAaaBBB
  zy23 <- 2*PAAABBB*PaaaBBb+2*PAAaBBB*PAaaBBb+2*PAaaBBB*PAAaBBb+2*PaaaBBB*PAAABBb
  zy24 <- 2*PAAABBB*PaaaBbb+2*PAAaBBB*PAaaBbb+2*PAaaBBB*PAAaBbb+2*PaaaBBB*PAAABbb+2*PAAABBb*PaaaBBb+2*PAAaBBb*PAaaBBb
  zy25 <- 2*PAAABBB*Paaabbb+2*PAAaBBB*PAaabbb+2*PAaaBBB*PAAabbb+2*PaaaBBB*PAAAbbb+2*PAAABBb*PaaaBbb+2*PAAaBBb*PAaaBbb+2*PAaaBBb*PAAaBbb+2*PaaaBBb*PAAABbb
  zy26 <- 2*PAAABBb*Paaabbb+2*PAAaBBb*PAaabbb+2*PAaaBBb*PAAabbb+2*PaaaBBb*PAAAbbb+2*PAAABbb*PaaaBbb+2*PAAaBbb*PAaaBbb
  zy27 <- 2*PAAABbb*Paaabbb+2*PAAaBbb*PAaabbb+2*PAaaBbb*PAAabbb+2*PaaaBbb*PAAAbbb
  zy28 <- 2*PAAAbbb*Paaabbb+2*PAAabbb*PAaabbb

  zy29 <- (PAaaBBB)^2+2*PAAaBBB*PaaaBBB
  zy30 <- 2*PAAaBBB*PaaaBBb+2*PAaaBBB*PAaaBBb+2*PaaaBBB*PAAaBBb
  zy31 <- (PAaaBBb)^2+2*PAAaBBB*PaaaBbb+2*PAaaBBB*PAaaBbb+2*PaaaBBB*PAAaBbb+2*PAAaBBb*PaaaBBb
  zy32 <- 2*PAAaBBB*Paaabbb+2*PAaaBBB*PAaabbb+2*PaaaBBB*PAAabbb+2*PAAaBBb*PaaaBbb+2*PAaaBBb*PAaaBbb+2*PaaaBBb*PAAaBbb
  zy33 <- (PAaaBbb)^2+2*PAAaBBb*Paaabbb+2*PAaaBBb*PAaabbb+2*PaaaBBb*PAAabbb+2*PAAaBbb*PaaaBbb
  zy34 <- 2*PAAaBbb*Paaabbb+2*PAaaBbb*PAaabbb+2*PaaaBbb*PAAabbb
  zy35 <- (PAaabbb)^2+2*PAAabbb*Paaabbb

  zy36 <- 2*PAaaBBB*PaaaBBB
  zy37 <- 2*PAaaBBB*PaaaBBb+2*PaaaBBB*PAaaBBb
  zy38 <- 2*PAaaBBB*PaaaBbb+2*PaaaBBB*PAaaBbb+2*PAaaBBb*PaaaBBb
  zy39 <- 2*PAaaBBB*Paaabbb+2*PaaaBBB*PAaabbb+2*PAaaBBb*PaaaBbb+2*PaaaBBb*PAaaBbb
  zy40 <- 2*PAaaBBb*Paaabbb+2*PaaaBBb*PAaabbb+2*PAaaBbb*PaaaBbb
  zy41 <- 2*PAaaBbb*Paaabbb+2*PaaaBbb*PAaabbb
  zy42 <- 2*PAaabbb*Paaabbb

  zy43 <- (PaaaBBB)^2
  zy44 <- 2*PaaaBBB*PaaaBBb
  zy45 <- (PaaaBBb)^2+2*PaaaBBB*PaaaBbb
  zy46 <- 2*PaaaBBB*Paaabbb+2*PaaaBBb*PaaaBbb
  zy47 <- (PaaaBbb)^2+2*PaaaBBb*Paaabbb
  zy48 <- 2*PaaaBbb*Paaabbb
  zy49 <- (Paaabbb)^2

  tmp <- c(zy1,zy2,zy3,zy4,zy5,zy6,zy7,zy8,zy9,zy10,
           zy11,zy12,zy13,zy14,zy15,zy16,zy17,zy18,zy19,
           zy20,zy21,zy22,zy23,zy24,zy25,zy26,zy27,zy28,
           zy29,zy30,zy31,zy32,zy33,zy34,zy35,zy36,zy37,
           zy38,zy39,zy40,zy41,zy42,zy43,zy44,zy45,
           zy46,zy47,zy48,zy49)
  #sum(tmp)
  return(tmp)
}
