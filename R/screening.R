#' Add AutoLD screening summaries to a genome-scan result
add_autold_screening <- function(x, r2_threshold = 0.05, ho_threshold = 0.80) {
  x <- as.data.frame(x)
  if (!"r2_comp" %in% names(x)) stop("x must contain r2_comp.")
  ho_cols <- grep("^D_norm_[0-9]+[0-9]+$", names(x), value = TRUE)
  ho_cols <- setdiff(ho_cols, "D_norm_11")
  if (!length(ho_cols)) stop("No normalized higher-order D columns were found.")

  mat <- abs(as.matrix(x[, ho_cols, drop = FALSE]))
  mat[!is.finite(mat)] <- NA_real_
  H_HO <- apply(mat, 1L, function(z) if (all(is.na(z))) NA_real_ else max(z, na.rm = TRUE))
  x$H_HO <- H_HO

  high_d11 <- is.finite(x$r2_comp) & x$r2_comp >= r2_threshold
  high_ho <- is.finite(x$H_HO) & x$H_HO >= ho_threshold
  x$ScreenClass <- ifelse(
    high_d11 & high_ho, "Both-high",
    ifelse(high_d11, "D11-only", ifelse(high_ho, "HO-only", "Neither"))
  )
  x
}
