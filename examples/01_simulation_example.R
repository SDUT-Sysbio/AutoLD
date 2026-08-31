## AutoLD simulation examples
library(AutoLD)

n_cores <- max(1L, parallel::detectCores(logical = FALSE) - 1L)

## ------------------------------------------------------------
## 1. Autotetraploid: D11 + higher-order D12
## ------------------------------------------------------------
D4 <- matrix(0, nrow = 2, ncol = 2)
D4[1, 1] <- 0.05
D4[1, 2] <- 0.03

## One simulated dataset and one-pair fit
sim4 <- simulate_from_eta(
  v_ploidy = 4,
  pA = 0.6, pB = 0.5,
  eta1_A = 0.15, eta1_B = 0.05,
  D_matrix = D4,
  n = 1000,
  type = "full",
  seed = 123
)

fit4 <- AutoLD(
  sim4$zygote_A,
  sim4$zygote_B,
  ploidy = 4,
  type = "full",
  calc_norm = TRUE,
  permutations = 0
)
print(fit4)

## Repeated evaluation
res4_full <- evaluate_AutoLD_pipeline(
  reps = 100,
  n = 1000,
  c_ploidy = 2,
  type = "full",
  pA = 0.6, pB = 0.5,
  eta_A = 0.15, eta_B = 0.05,
  D_true_matrix = D4,
  calc_norm = TRUE,
  calc_pvalue = FALSE,
  n_cores = n_cores,
  seed = 123
)

## Coarsened/partial dosage example
res4_partial <- evaluate_AutoLD_pipeline(
  reps = 100,
  n = 1000,
  c_ploidy = 2,
  type = "partial",
  pA = 0.6, pB = 0.5,
  eta_A = 0.15, eta_B = 0.05,
  D_true_matrix = D4,
  calc_norm = TRUE,
  calc_pvalue = FALSE,
  n_cores = n_cores,
  seed = 123
)

## To evaluate permutation-calibrated power, use for example:
## calc_pvalue = TRUE, B_perm = 1000.

## ------------------------------------------------------------
## 2. Autohexaploid: D11 + D22 + D31
## ------------------------------------------------------------
D6 <- matrix(0, nrow = 3, ncol = 3)
D6[1, 1] <- 0.05
D6[2, 2] <- 0.03
D6[3, 1] <- 0.02

res6_full <- evaluate_AutoLD_pipeline(
  reps = 100,
  n = 1000,
  c_ploidy = 3,
  type = "full",
  pA = 0.6, pB = 0.5,
  eta_A = 0.15, eta_B = 0.05,
  D_true_matrix = D6,
  calc_norm = TRUE,
  calc_pvalue = FALSE,
  n_cores = n_cores,
  seed = 123
)
