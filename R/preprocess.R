#' Thin SNPs by minimum physical distance
#'
#' Retains markers sequentially within each chromosome so that adjacent retained
#' markers are separated by at least min_dist_kb.
thin_snps_by_dist <- function(geno_df, info_df, min_dist_kb = 0.2) {
  if (!all(c("CHROM", "POS") %in% names(info_df))) {
    stop("info_df must contain CHROM and POS columns.")
  }
  if (nrow(geno_df) != nrow(info_df)) {
    stop("geno_df and info_df must have the same number of marker rows.")
  }
  if (!is.numeric(min_dist_kb) || length(min_dist_kb) != 1L || min_dist_kb < 0) {
    stop("min_dist_kb must be a non-negative numeric scalar.")
  }

  min_bp <- min_dist_kb * 1000
  keep_global <- integer(0)

  for (chr in unique(info_df$CHROM)) {
    idx <- which(info_df$CHROM == chr)
    idx <- idx[order(info_df$POS[idx])]
    if (!length(idx)) next

    keep_chr <- idx[1]
    last_pos <- info_df$POS[idx[1]]
    if (length(idx) > 1L) {
      for (k in idx[-1]) {
        if ((info_df$POS[k] - last_pos) >= min_bp) {
          keep_chr <- c(keep_chr, k)
          last_pos <- info_df$POS[k]
        }
      }
    }
    keep_global <- c(keep_global, keep_chr)
  }

  keep_global <- keep_global[order(match(info_df$CHROM[keep_global], unique(info_df$CHROM)),
                                   info_df$POS[keep_global])]
  list(
    geno = geno_df[keep_global, , drop = FALSE],
    info = info_df[keep_global, , drop = FALSE],
    kept_index = keep_global
  )
}
