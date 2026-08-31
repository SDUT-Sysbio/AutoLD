
library(parallel)

source("tensor_core.R")
source("simulate.R")
source("segregation_kernel.R")
source("estimate_wrapper.R")
source("LD_statistics.R")

library(Rcpp)
sourceCpp("em_core.cpp") 

# 构造一个真实的 4x 张量矩阵 (注入 D11=0.05 和 高阶 D12=0.03)
my_D_4x <- matrix(0, nrow=2, ncol=2)
my_D_4x[1, 1] <- 0.05
my_D_4x[1, 2] <- 0.03

res_fast_4x <- evaluate_AutoLD_pipeline(reps = 100, n = 1000, c_ploidy = 2, type = "full",pA = 0.6, pB = 0.5, eta_A = 0.15, eta_B = 0.05, 
  D_true_matrix = my_D_4x,calc_norm = T, calc_pvalue = T, n_cores = 20)

res_fast_4xp <- evaluate_AutoLD_pipeline(reps = 100, n = 1000, c_ploidy = 2, type = "partial",pA = 0.6, pB = 0.5, eta_A = 0.15, eta_B = 0.05, 
                                        D_true_matrix = my_D_4x,calc_norm = T, calc_pvalue = F, n_cores = 20)

#构造一个真实的 6x 张量矩阵 (3x3)
my_D_6x <- matrix(0, nrow=3, ncol=3)

my_D_6x[1, 1] <- 0.05
my_D_6x[2, 2] <- 0.03
my_D_6x[3, 1] <- 0.02

res_fast_6x <- evaluate_AutoLD_pipeline(reps = 100, n = 1000, c_ploidy = 3, type = "full", pA = 0.6, pB = 0.5, eta_A = 0.15, eta_B = 0.05, 
  D_true_matrix = my_D_6x,  calc_norm = TRUE,  calc_pvalue = FALSE, n_cores = 1)







