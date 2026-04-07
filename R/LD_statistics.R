
#' @noRd
calc_loglik_core <- function(params, geno_counts, prob_func,hwdA,hwdB) {
  if (any(params < -3.0) || any(params > 3.0)) return(-Inf)
  probs <- tryCatch(prob_func(params,hwdA=hwdA,hwdB=hwdB), error = function(e) NULL)
  if (is.null(probs) || any(is.nan(probs)) || any(is.na(probs))) return(-Inf)
  probs[probs < 1e-12] <- 1e-12
  return(sum(geno_counts * log(probs)))
}

# Initial optimized implementation of the LD test
LD_test_optimized <- function(res, geno, ploidy=4, type="full",optim_method="BFGS",hwdA=0,hwdB=0) {

  if (ploidy == 4) {
    if (type == "full")    { zyn_g <- DF4_cLD;  ii <- 25 }
    else                   { zyn_g <- DF4_cLDp; ii <- 9  }
  } else if (ploidy == 6) {
    if (type == "full")    { zyn_g <- DF6_cLD;  ii <- 49 }
    else                   { zyn_g <- DF6_cLDp; ii <- 9  }
  } else {
    stop("Invalid ploidy/type input")
  }

  nt <- table(geno)
  nts <- rep(0, ii)
  idx <- as.numeric(names(nt))
  valid_idx <- idx <= ii & idx >= 1
  if(any(valid_idx)) nts[idx[valid_idx]] <- as.numeric(nt)[valid_idx]


  logL_full <- calc_loglik_core(res, nts, zyn_g,hwdA = 0,hwdB=0)

  para_n <- length(res)
  pv_list <- c()

  target_indices <- 3:para_n

  for (i in target_indices) {

    fn_null <- function(params_sub) {
      current_params <- numeric(para_n)

      current_params[-i] <- params_sub
      current_params[i] <- 0

      return(-calc_loglik_core(current_params, nts, zyn_g,hwdA=hwdA,hwdB=hwdB))
    }

    start_params <- res[-i]

    fit_null <- optim(par = start_params, fn = fn_null, method = optim_method, control = list(maxit = 1000))

    logL_null <- -fit_null$value

    LR_stat <- 2 * (logL_full - logL_null)
    if (LR_stat < 0) LR_stat <- 0

    p_val <- pchisq(LR_stat, df = 1, lower.tail = FALSE)
    pv_list <- c(pv_list, p_val)
  }

  names(pv_list) <- names(res)[target_indices]
  return(pv_list)
}

## Legacy implementation of optimized LD test using specified optimization methods
LD_test_polished <- function(res, geno, ploidy=4, type="full", optim_method="BFGS",hwdA=0,hwdB=0) {


  if (ploidy == 4) {
    if (type == "full")    { zyn_g <- DF4_cLD;  ii <- 25 }
    else                   { zyn_g <- DF4_cLDp; ii <- 9  }
  } else if (ploidy == 6) {
    if (type == "full")    { zyn_g <- DF6_cLD;  ii <- 49 }
    else                   { zyn_g <- DF6_cLDp; ii <- 9  }
  } else {
    stop("Invalid ploidy/type input")
  }

  nt <- table(geno)
  nts <- rep(0, ii)
  idx <- as.numeric(names(nt))
  valid_idx <- idx <= ii & idx >= 1
  if(any(valid_idx)) nts[idx[valid_idx]] <- as.numeric(nt)[valid_idx]

  fn_full <- function(params) {
    return(-calc_loglik_core(params, nts, zyn_g,hwdA=hwdA,hwdB=hwdB))
  }


  fit_full <- optim(par = res, fn = fn_full, method = optim_method, control = list(maxit = 1000))

  res_polished <- fit_full$par    # 这是更精确的 MLE
  logL_full <- -fit_full$value    # 这是 Optim 能找到的最高山峰

  para_n <- length(res_polished)
  pv_list <- c()
  target_indices <- 3:para_n

  for (i in target_indices) {

    fn_null <- function(params_sub) {
      current_params <- numeric(para_n)

      current_params[-i] <- params_sub
      current_params[i] <- 0

      return(-calc_loglik_core(current_params, nts, zyn_g,hwdA = 0,hwdB=0))
    }

    start_params <- res_polished[-i]

    fit_null <- optim(par = start_params, fn = fn_null, method = optim_method, control = list(maxit = 1000))

    logL_null <- -fit_null$value

    diff_L <- logL_full - logL_null

    if (diff_L < 0) {
      LR_stat <- 0
    } else {
      LR_stat <- 2 * diff_L
    }

    p_val <- pchisq(LR_stat, df = 1, lower.tail = FALSE)
    pv_list <- c(pv_list, p_val)
  }

  names(pv_list) <- names(res)[target_indices]
  return(pv_list)
}


LD_test_ignore_DR <- function(res, geno, ploidy=4, type="full", optim_method="BFGS",hwdA=0,hwdB=0) {

  if (ploidy == 4) {
    if (type == "full")    { zyn_g <- DF4_cLD;  ii <- 25 }
    else                   { zyn_g <- DF4_cLDp; ii <- 9  }
    alpha_idx <- c(3, 4)
  } else if (ploidy == 6) {
    if (type == "full")    { zyn_g <- DF6_cLD;  ii <- 49 }
    else                   { zyn_g <- DF6_cLDp; ii <- 9  }
    alpha_idx <- c(3, 4)
  } else {
    stop("Invalid ploidy/type input")
  }

  nt <- table(geno)
  nts <- rep(0, ii)
  idx <- as.numeric(names(nt))
  valid_idx <- idx <= ii & idx >= 1
  if(any(valid_idx)) nts[idx[valid_idx]] <- as.numeric(nt)[valid_idx]

  start_params_no_alpha <- res[-alpha_idx]

  fn_full_no_alpha <- function(params_sub) {

    full_params <- numeric(length(res))

    cnt <- 1
    for(k in 1:length(full_params)){
      if(k %in% alpha_idx){
        full_params[k] <- 0 # 强制双减为 0
      } else {
        full_params[k] <- params_sub[cnt]
        cnt <- cnt + 1
      }
    }

    return(-calc_loglik_core(full_params, nts, zyn_g,hwdA = 0,hwdB=0))
  }

  fit_full <- optim(par = start_params_no_alpha, fn = fn_full_no_alpha,
                    method = optim_method, control = list(maxit = 1000))

  logL_full <- -fit_full$value
  res_full_optimized <- fit_full$par


  pv_list <- c()
  bias_list <- c()

  target_indices_original <- 5:length(res)

  for (i_orig in target_indices_original) {

    i_new <- i_orig - 2

    fn_null_no_alpha <- function(params_sub_null) {

      current_params_no_alpha <- numeric(length(start_params_no_alpha))

      current_params_no_alpha[-i_new] <- params_sub_null
      current_params_no_alpha[i_new] <- 0

      full_params <- numeric(length(res))
      cnt <- 1
      for(k in 1:length(full_params)){
        if(k %in% alpha_idx){
          full_params[k] <- 0
        } else {
          full_params[k] <- current_params_no_alpha[cnt]
          cnt <- cnt + 1
        }
      }
      return(-calc_loglik_core(full_params, nts, zyn_g,hwdA=hwdA,hwdB=hwdB))
    }

    start_params_null <- res_full_optimized[-i_new]

    fit_null <- optim(par = start_params_null, fn = fn_null_no_alpha,
                      method = optim_method, control = list(maxit = 2000))

    logL_null <- -fit_null$value

    LR_stat <- 2 * (logL_full - logL_null)

    if (LR_stat < -1e-6) {
      warning(paste("Negative LR detected:", LR_stat, "Check optimization convergence."))
      LR_stat <- 0
    } else if (LR_stat < 0) {
      LR_stat <- 0
    }

    p_val <- pchisq(LR_stat, df = 1, lower.tail = FALSE)
    pv_list <- c(pv_list, p_val)

  }

  names(pv_list) <- names(res)[target_indices_original]

  return(list(p_val = pv_list, est_params = res_full_optimized))
}

#' Highly Optimized Likelihood Ratio Test for LD (Accounting for Double Reduction)
#'
#' The final, most highly optimized implementation of the Likelihood Ratio Test (LRT)
#' to evaluate the statistical significance of LD in autopolyploids. This function
#' natively accounts for double reduction rates estimated from the data.
#'
#' @param res A numeric vector of parameter estimates under the alternative hypothesis (H1).
#' @param geno A numeric vector or data frame column of joint genotype classes.
#' @param ploidy An integer specifying the ploidy level (e.g., 4 or 6). Default is 4.
#' @param type A character string specifying the resolution model: \code{"full"} or \code{"partial"}. Default is \code{"full"}.
#' @param hwdA Numeric. Hardy-Weinberg Disequilibrium parameter for locus A. Default is 0.
#' @param hwdB Numeric. Hardy-Weinberg Disequilibrium parameter for locus B. Default is 0.
#' @return A named numeric vector of p-values corresponding to the tested LD parameters.
#' @export

LD_test_optimized_final <- function(res, geno, ploidy=4, type="full",hwdA=0,hwdB=0) {

  if (ploidy == 4) {
    if (type == "full")    { zyn_g <- DF4_cLD;  ii <- 25 }
    else                   { zyn_g <- DF4_cLDp; ii <- 9  }
  } else if (ploidy == 6) {
    if (type == "full")    { zyn_g <- DF6_cLD;  ii <- 49 }
    else                   { zyn_g <- DF6_cLDp; ii <- 9  }
  } else {
    stop("Invalid ploidy/type input")
  }

  nt <- table(geno)
  nts <- rep(0, ii)
  idx <- as.numeric(names(nt))
  valid_idx <- idx <= ii & idx >= 1
  if(any(valid_idx)) nts[idx[valid_idx]] <- as.numeric(nt)[valid_idx]


  logL_full <- calc_loglik_core(res, nts, zyn_g,hwdA = hwdB,hwdB=hwdB)

  para_n <- length(res)
  pv_list <- c()

  target_indices <- 3:para_n

  for (i in target_indices) {

    fn_null <- function(params_sub) {
      current_params <- numeric(para_n)

      current_params[-i] <- params_sub
      current_params[i] <- 0

      return(-calc_loglik_core(current_params, nts, zyn_g,hwdA=hwdA,hwdB=hwdB))
    }

    start_params <- res[-i]

    fit_null <- try(optim(par = start_params, fn = fn_null, method = "BFGS", control = list(maxit = 1000)),silent = T)
    if(inherits(fit_null, "try-error")){
      fit_null <- try(optim(par = start_params, fn = fn_null, method = "Nelder-Mead", control = list(maxit = 1000)),silent = T)
    }
    if(inherits(fit_null, "try-error")){
      pv_list <- c(pv_list, NA)
    }else{
      logL_null <- -fit_null$value

      LR_stat <- 2 * (logL_full - logL_null)
      if (LR_stat < 0) LR_stat <- 0

      p_val <- pchisq(LR_stat, df = 1, lower.tail = FALSE)
      pv_list <- c(pv_list, p_val)
    }


  }

  names(pv_list) <- names(res)[target_indices]
  return(pv_list)
}


#' Normalize Linkage Disequilibrium Parameters
#'
#' This function standardizes the raw Linkage Disequilibrium (LD) parameter D into
#' comparable metrics such as D' (D prime) and r-squared, adjusting for the specific
#' ploidy level and estimated allele frequencies.
#'
#' @param res_vec A numeric vector of estimated parameters from the estimation functions
#' (e.g., auto4_est or auto6_est).
#' @param ploidy An integer specifying the ploidy level (e.g., 4 or 6). Default is 6.
#' @return A data frame containing the standardized LD statistics, including D_over_SD
#' and D_Scaled_Clamped for various LD components.
#' @export

normalize_LD <- function(res_vec, ploidy=6) {

  get_val <- function(name) {
    if(name %in% names(res_vec)) return(as.numeric(res_vec[name]))
    return(NA_real_)
  }
  pA <- get_val("pA"); if(is.na(pA)) return(NULL)
  pB <- get_val("pB"); if(is.na(pB)) return(NULL)
  DA_val <- get_val("aA"); if(is.na(DA_val)) DA_val <- 0
  DB_val <- get_val("aB"); if(is.na(DB_val)) DB_val <- 0

  full_map <- data.frame(
    Name = c("Deab", "DAb", "DaB", "DAB", "DAAb", "DaBB", "DAAB", "DABB", "DAABB"),
    Order_A = c(1, 2, 1, 2,  3, 1, 3, 2, 3),
    Order_B = c(1, 1, 2, 2,  1, 3, 2, 3, 3)
  )
  if (ploidy == 4) {
    target_names <- c("Deab", "DAb", "DaB", "DAB")
    param_map <- full_map[full_map$Name %in% target_names, ]
  } else {
    param_map <- full_map
  }


  get_SD_robust <- function(p, D_hwd) {
    var_base <- p * (1 - p) + D_hwd
    if (is.na(var_base) || var_base <= 1e-9) {
      var_base <- p * (1 - p)
    }
    if (var_base <= 1e-9) return(NA_real_)
    return(sqrt(var_base))
  }
  sd_A <- get_SD_robust(pA, DA_val)
  sd_B <- get_SD_robust(pB, DB_val)

  calc_D_over_SD <- function(D_val) {
    if (is.na(sd_A) || is.na(sd_B)) return(NA_real_)
    return(D_val / (sd_A * sd_B))
  }

  calc_D_Scaled_Components <- function(D_val, ordA, ordB) {
    if (is.na(D_val)) return(list(raw=NA_real_, clamped=NA_real_))
    if (abs(D_val) < 1e-12) return(list(raw=0, clamped=0))

    get_freq_proxy <- function(p, D_hwd, ord) {
      base_p <- p; base_q <- 1 - p
      if(ord == 2) { f_pos <- base_p^2 + D_hwd; f_neg <- base_q^2 + D_hwd }
      else if (ord == 3) { f_pos <- base_p^3; f_neg <- base_q^3 }
      else { f_pos <- base_p; f_neg <- base_q }

      f_pos <- max(0, f_pos)
      f_neg <- max(0, f_neg)
      s <- f_pos + f_neg

      if (s < 1e-9) {
        if(ord == 2) { f_pos <- base_p^2; f_neg <- base_q^2 }
        else if(ord == 3) { f_pos <- base_p^3; f_neg <- base_q^3 }
        else { f_pos <- base_p; f_neg <- base_q }
        s <- f_pos + f_neg
      }

      return(c(pos=f_pos/s, neg=f_neg/s))
    }

    fA <- get_freq_proxy(pA, DA_val, ordA)
    fB <- get_freq_proxy(pB, DB_val, ordB)

    limit <- NA_real_
    if (D_val > 0) limit <- min(fA['pos'] * fB['neg'], fA['neg'] * fB['pos'])
    else limit <- min(fA['pos'] * fB['pos'], fA['neg'] * fB['neg'])

    if (is.na(limit) || limit <= 1e-9) return(list(raw=NA_real_, clamped=NA_real_))

    raw_val <- D_val / limit
    clamped_val <- max(min(raw_val, 1), -1)
    return(list(raw=raw_val, clamped=clamped_val))
  }

  available_params <- intersect(names(res_vec), param_map$Name)
  results <- data.frame(Parameter = available_params, Value_Raw = NA, D_over_SD = NA, D_Scaled_Clamped = NA)

  for(i in 1:nrow(results)) {
    p_name <- results$Parameter[i]
    val <- as.numeric(res_vec[p_name])
    info <- param_map[param_map$Name == p_name, ]

    results$Value_Raw[i] <- val
    results$D_over_SD[i] <- calc_D_over_SD(val)
    d_out <- calc_D_Scaled_Components(val, info$Order_A, info$Order_B)
    results$D_Scaled_Clamped[i] <- d_out$clamped
  }
  return(results)
}









