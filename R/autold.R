#' Fit AutoLD to a single pair of markers
#'
#' @param dosage_A,dosage_B Integer allele dosages for the two loci.
#' @param ploidy Even somatic ploidy (currently intended for 4x and 6x analyses).
#' @param type "full" for resolved dosage or "partial" for collapsed 0/mid/v classes.
#' @param calc_norm Whether to compute fixed-margin feasible normalization.
#' @param permutations Number of marginal-preserving permutations. Set 0 to skip.
#' @param seed Optional random seed.
AutoLD <- function(dosage_A, dosage_B, ploidy = 4L,
                   type = c("full", "partial"), calc_norm = TRUE,
                   permutations = 0L, seed = NULL) {
  type <- match.arg(type)
  if (ploidy %% 2L != 0L || ploidy < 4L) stop("ploidy must be an even integer >= 4.")
  if (length(dosage_A) != length(dosage_B)) stop("dosage_A and dosage_B must have equal length.")

  valid <- !is.na(dosage_A) & !is.na(dosage_B)
  A <- as.integer(dosage_A[valid])
  B <- as.integer(dosage_B[valid])
  if (!length(A)) stop("No complete two-locus observations.")
  if (any(A < 0L | A > ploidy) || any(B < 0L | B > ploidy)) {
    stop("Dosage values must lie between 0 and ploidy.")
  }

  c_ploidy <- ploidy / 2L
  if (type == "partial") {
    A <- ifelse(A == 0L, 0L, ifelse(A == ploidy, 2L, 1L))
    B <- ifelse(B == 0L, 0L, ifelse(B == ploidy, 2L, 1L))
    geno <- data.frame(class = as.integer(A * 3L + B))
  } else {
    geno <- data.frame(class = as.integer(A * (ploidy + 1L) + B))
  }

  estimate <- estimate_LD_decoupled(geno, c_ploidy, type = type, calc_norm = calc_norm)
  perm <- NULL
  if (permutations > 0L) {
    perm <- AutoLD_Permutation_Test(
      geno = geno, c_ploidy = c_ploidy, type = type,
      B = as.integer(permutations), seed = seed, obs_res = estimate
    )
  }

  out <- list(
    call = match.call(),
    ploidy = ploidy,
    type = type,
    n = length(A),
    estimate = estimate,
    permutation = perm
  )
  class(out) <- "AutoLD"
  out
}

print.AutoLD <- function(x, ...) {
  cat("AutoLD two-locus analysis\n")
  cat(sprintf("  Ploidy: %dx | type: %s | N: %d\n", x$ploidy, x$type, x$n))
  cat(sprintf("  ESD A: %.6g | ESD B: %.6g\n",
              x$estimate$eta_A_est, x$estimate$eta_B_est))
  cat(sprintf("  r_comp: %.6g | r2_comp: %.6g\n",
              x$estimate$r_comp, x$estimate$r2_comp))
  if (!is.null(x$permutation)) {
    cat(sprintf("  P_D11: %.6g | P_HO: %.6g | P_All: %.6g\n",
                x$permutation$P_D11,
                x$permutation$P_HigherOrder,
                x$permutation$P_AllTensor))
  }
  invisible(x)
}
