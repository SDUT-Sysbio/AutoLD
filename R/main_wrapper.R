#' Optimized Parallel Pipeline for Chromosome-wide LD Estimation
#'
#' @param chr_name Character. The name of the chromosome to analyze (e.g., "LR999451.1").
#' @param geno_df A data frame or matrix of genotype dosages. Rows must represent markers (SNPs)
#' and columns must represent samples. Values should be integer dosages corresponding to the
#' ploidy level (e.g., 0, 1, 2, 3, 4 for tetraploids) or \code{NA} for missing data.
#' The number of rows must exactly match the number of rows in \code{info_df}.
#' @param info_df A data frame containing variant information. It must include at least two
#' mandatory columns: \code{CHROM} (chromosome ID) and \code{POS} (physical position in base pairs).
#' The row order must strictly align with \code{geno_df}.
#' @param method Character. The resolution model: \code{"full"} or \code{"partial"}. Default is \code{"full"}.
#' @param window_kb Numeric. The sliding window size in kilobases (kb). Default is 200.
#' @param n_cores Integer. Number of CPU cores for parallel processing. Default is 1.
#' @param ploidy Integer. The ploidy level (e.g., 4 or 6). Default is 4.
#' @return A data frame containing pairwise LD estimations, p-values, and standardized statistics.
#' @importFrom parallel makeCluster clusterExport parLapply stopCluster mclapply
#' @importFrom stats optim optimize pchisq dbinom
#' @export
#' @export
AutoLD_optimized <- function(chr_name, geno_df, info_df, method = "full", window_kb = 200, n_cores = 1, ploidy = 4) {

  message(paste0("[", chr_name, "] Starting analysis (Ploidy=", ploidy, ")..."))

  idx_chr <- which(info_df$CHROM == chr_name)
  if (length(idx_chr) < 2) {
    warning(paste("Chromosome", chr_name, "skipped: < 2 SNPs."))
    return(NULL)
  }

  sub_info <- info_df[idx_chr, ]
  sub_geno <- t(as.matrix(geno_df[idx_chr, ]))

  ord <- order(sub_info$POS)
  sub_info <- sub_info[ord, ]
  sub_geno <- sub_geno[, ord]
  positions <- sub_info$POS

  n_snps <- length(positions)
  message(paste0("[", chr_name, "] Loaded ", n_snps, " SNPs."))

  window_bp <- window_kb * 1000
  end_idxs <- findInterval(positions + window_bp, positions)
  idxs <- 1:n_snps
  pair_counts <- end_idxs - idxs

  valid_mask <- pair_counts > 0
  if (!any(valid_mask)) return(NULL)

  valid_i <- idxs[valid_mask]
  valid_counts <- pair_counts[valid_mask]

  i_vec <- rep(valid_i, valid_counts)
  j_vec <- i_vec + sequence(valid_counts)
  tasks_mat <- cbind(i_vec, j_vec)
  n_pairs <- nrow(tasks_mat)

  message(paste0("[", chr_name, "] Total pairs: ", n_pairs, " (Applying chunked parallelization)"))

  n_chunks <- n_cores * 4
  chunk_factor <- cut(seq_len(n_pairs), breaks = n_chunks, labels = FALSE)
  task_chunks <- split(seq_len(n_pairs), chunk_factor)


  process_chunk <- function(row_indices) {
    chunk_res <- vector("list", length(row_indices))

    for (k in seq_along(row_indices)) {
      real_row <- row_indices[k]
      i <- tasks_mat[real_row, 1]
      j <- tasks_mat[real_row, 2]

      g1_raw <- sub_geno[, i]
      g2_raw <- sub_geno[, j]
      valid <- !is.na(g1_raw) & !is.na(g2_raw)

      if (sum(valid) < 100) next

      g1 <- g1_raw[valid]
      g2 <- g2_raw[valid]

      if(method == "full"){
        sample_class_ids <- g1 + 1 + (g2 * (ploidy+1))
      } else {
        sample_class_ids <- g1 + 1 + (g2 * 3)
      }

      res_list <- tryCatch({
        if (ploidy == 4) {
          auto4_est(geno=data.frame(class=sample_class_ids), type=method)
        } else {
          auto6_est(geno=data.frame(class=sample_class_ids), type=method)
        }
      }, error = function(e) return(NULL))

      if (is.null(res_list)) next

      res_vec <- unlist(res_list)

      pv <- LD_test_optimized_final(res=res_vec, geno=sample_class_ids, ploidy=ploidy, type=method)
      pv_vec <- as.numeric(pv)
      if (!is.null(names(pv))) {
        names(pv_vec) <- paste0("P_", names(pv))
      } else {
        names(pv_vec) <- paste0("P_Test", 1:length(pv_vec))
      }

      norm_df <- normalize_LD(res_vec, ploidy = ploidy)
      if (is.null(norm_df) || nrow(norm_df) == 0) next

      stats_vec <- c()
      for(m in 1:nrow(norm_df)) {
        pname <- norm_df$Parameter[m]
        stats_vec[paste0(pname, "_r")]  <- norm_df$D_over_SD[m]
        stats_vec[paste0(pname, "_DS")] <- norm_df$D_Scaled_Clamped[m]
      }

      chunk_res[[k]] <- c(
        list(CHROM = chr_name, POS1 = positions[i], POS2 = positions[j], Dist_bp = positions[j] - positions[i], N_Samples = sum(valid)),
        as.list(pv_vec),
        as.list(res_vec),
        as.list(stats_vec)
      )
    }

    chunk_res <- chunk_res[!sapply(chunk_res, is.null)]
    if (length(chunk_res) == 0) return(NULL)

    return(data.table::rbindlist(chunk_res, fill = TRUE))
  }

  if (Sys.info()[['sysname']] == "Windows") {
    cl <- makeCluster(n_cores)
    clusterExport(cl, varlist = c("sub_geno", "positions", "chr_name", "method", "ploidy", "tasks_mat",
                                  "auto4_est", "auto6_est", "auto4_est_LDp",
                                  "normalize_LD", "LD_test_optimized_final"), envir = environment())
    res_list_final <- parLapply(cl, task_chunks, process_chunk)
    stopCluster(cl)
  } else {
    # Linux/Mac
    res_list_final <- mclapply(task_chunks, process_chunk, mc.cores = n_cores)
  }

  res_list_final <- res_list_final[!sapply(res_list_final, is.null)]
  if(length(res_list_final) == 0) return(NULL)

  final_df <- data.table::rbindlist(res_list_final, fill = TRUE)

  return(as.data.frame(final_df))
}
