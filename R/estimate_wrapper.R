# ==============================================================================
# 完整终极版：estimate_wrapper.R (纯净正交解耦 + MLE + 收敛监控)
# ==============================================================================

# 辅助函数 1: 计算类别
make_counts <- function(dose, n_classes) {
  cls <- as.integer(as.character(dose))
  if (any(cls < 0 | cls >= n_classes, na.rm = TRUE)) {
    stop("Dose class strictly out of bounds [0, n_classes - 1]")
  }
  tab <- tabulate(cls + 1, nbins = n_classes)
  return(tab)
}


# 核心模块 1: 参数化分离畸变估计 (极大似然 MLE)
estimate_etas_from_data <- function(geno, c_ploidy, type = "full") {
  v <- 2 * c_ploidy
  S_B <- ifelse(type == "full", v + 1, 3)
  locus_A_dosage <- geno$class %/% S_B
  locus_B_dosage <- geno$class %% S_B
  levels_vec <- if (type == "full") 0:v else 0:2
  
  fit_single <- function(dosage_vec) {
    dosage_vec <- dosage_vec[!is.na(dosage_vec)]
    if(length(dosage_vec) == 0) return(c(p = 0.5, eta = 0))
    
    obs_counts <- table(factor(dosage_vec, levels = levels_vec))
    n_k <- as.numeric(obs_counts)
    if(sum(n_k) == 0) return(c(p = 0.5, eta = 0))
    
    if (type == "full") {
      # 完整模式：1D 轮廓似然搜索
      p_hat <- sum((0:v) * n_k) / (v * sum(n_k))
      if(p_hat < 1e-4 || p_hat > 1 - 1e-4) return(c(p = p_hat, eta = 0))
      
      loss_1d <- function(eta_cand) {
        Q_theo <- get_equilibrium_Q_one_eta(v, p = p_hat, eta1 = eta_cand)
        Z_theo <- pmax(convolve_gametes(Q_theo), 1e-12)
        return(-sum(n_k * log(Z_theo)))
      }
      opt <- optimize(loss_1d, interval = c(0, 0.33), tol = 1e-6)
      return(c(p = p_hat, eta = opt$minimum))
      
    } else {
      # Coarsened 模式：2D 空间联合 MLE
      loss_2d <- function(params) {
        p_cand <- params[1]; eta_cand <- params[2]
        Q_theo <- get_equilibrium_Q_one_eta(v, p = p_cand, eta1 = eta_cand)
        Z_theo <- convolve_gametes(Q_theo)
        Z_partial <- pmax(c(Z_theo[1], sum(Z_theo[2:v]), Z_theo[v+1]), 1e-12)
        return(-sum(n_k * log(Z_partial)))
      }
      
      p_init <- (n_k[2] * 0.5 + n_k[3] * 1.0) / sum(n_k)
      p_init <- pmax(1e-3, pmin(1 - 1e-3, p_init))
      
      opt <- optim(c(p_init, 0.0), loss_2d, method = "L-BFGS-B",
                   lower = c(1e-4, 0), upper = c(1 - 1e-4, 0.33))
      return(c(p = opt$par[1], eta = opt$par[2]))
    }
  }
  
  res_A <- fit_single(locus_A_dosage)
  res_B <- fit_single(locus_B_dosage)
  
  return(list(eta_A = res_A["eta"], eta_B = res_B["eta"],
              p_A = res_A["p"], p_B = res_B["p"]))
}

# 核心模块 2: 纯净 EM 引擎包装器
run_pure_EM <- function(geno, c_ploidy, type = "full", tol = 1e-8) {
  v_ploidy <- c_ploidy * 2
  n_classes <- if (type == "full") (v_ploidy + 1)^2 else 9
  nts <- make_counts(geno$class, n_classes)
  
  em_res <- if(type == "full") {
    EM_full_cpp(nts, c_ploidy, tol = tol)
  } else {
    EM_partial_cpp(nts, c_ploidy, tol = tol)
  }
  
  em_vector <- em_res$pi
  P_AB <- matrix(em_vector, nrow = c_ploidy + 1, ncol = c_ploidy + 1, byrow = TRUE)
  P_AB[P_AB < 0] <- 0
  P_AB <- P_AB / sum(P_AB)
  
  return(list(
    P_AB = P_AB,
    converged = em_res$converged,
    monotone = em_res$likelihood_monotone,
    iterations = em_res$iterations,
    logLik = em_res$logLik
  ))
}


estimate_LD_decoupled <- function(geno, c_ploidy, type = "full", calc_norm = TRUE) {
  em_out <- run_pure_EM(geno, c_ploidy, type)
  P_AB <- em_out$P_AB
  is_converged <- em_out$converged
  N_eff <- sum(make_counts(geno$class, if (type == "full") (2*c_ploidy + 1)^2 else 9))
  
  QA <- rowSums(P_AB); QB <- colSums(P_AB)
  Phi_A <- build_orthogonal_basis(QA, c_ploidy)
  Phi_B <- build_orthogonal_basis(QB, c_ploidy)
  
  D_tensor <- extract_Drs_tensor(P_AB, c_ploidy)
  
  D_norm_tensor <- matrix(NA_real_, c_ploidy, c_ploidy)
  D_plus_tensor <- matrix(NA_real_, c_ploidy, c_ploidy)
  D_minus_tensor <- matrix(NA_real_, c_ploidy, c_ploidy)
  
  if (is_converged && calc_norm) {
    for(r in 1:c_ploidy) {
      for(s in 1:c_ploidy) {
        norm_res <- normalize_tensor_LD(D_tensor[r,s], r, s, Phi_A, Phi_B, QA, QB)
        D_norm_tensor[r,s] <- norm_res$D_norm
        D_plus_tensor[r,s] <- norm_res$D_plus
        D_minus_tensor[r,s] <- norm_res$D_minus
      }
    }
  }
  
  etas <- estimate_etas_from_data(geno, c_ploidy, type)
  dose <- 0:c_ploidy
  mean_A <- sum(dose * QA); mean_B <- sum(dose * QB)
  cov_gam <- sum(outer(dose, dose, "*") * P_AB) - (mean_A * mean_B)
  var_A <- sum((dose - mean_A)^2 * QA); var_B <- sum((dose - mean_B)^2 * QB)
  r_gam <- cov_gam / sqrt(var_A * var_B + 1e-12)
  
  res_list <- list(
    N_eff = N_eff,
    eta_A_est = etas$eta_A, eta_B_est = etas$eta_B,
    pA = etas$p_A, pB = etas$p_B,
    D_comp = cov_gam / c_ploidy, r_comp = r_gam, r2_comp = r_gam^2,
    em_converged = is_converged
  )
  
  for(r in 1:c_ploidy) {
    for(s in 1:c_ploidy) {
      res_list[[paste0("D_", r, s)]] <- D_tensor[r,s]
      res_list[[paste0("D_norm_", r, s)]] <- D_norm_tensor[r,s]
      res_list[[paste0("D_plus_", r, s)]] <- D_plus_tensor[r,s]
      res_list[[paste0("D_minus_", r, s)]] <- D_minus_tensor[r,s]
    }
  }
  return(res_list)
}