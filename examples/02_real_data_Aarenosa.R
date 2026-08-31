## AutoLD empirical genome-scan example: autotetraploid Arabidopsis arenosa
library(AutoLD)

## Replace these two paths with the local filtered dosage and marker-information files.
geno_file <- "arenosa_dosage_0to4.txt"
info_file <- "arenosa_markers_filter.txt"

Arenosa_gen <- read.table(geno_file, header = TRUE, check.names = FALSE)
Arenosa_info <- read.table(info_file, header = TRUE)

stopifnot(nrow(Arenosa_gen) == nrow(Arenosa_info))
stopifnot(all(c("CHROM", "POS") %in% names(Arenosa_info)))

## Manuscript analysis: minimum inter-marker distance = 0.2 kb.
Arenosa_sparse <- thin_snps_by_dist(
  geno_df = Arenosa_gen,
  info_df = Arenosa_info,
  min_dist_kb = 0.2
)

output_dir <- "AutoLD_Results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

all_chroms <- unique(Arenosa_sparse$info$CHROM)
n_cores <- max(1L, parallel::detectCores(logical = FALSE) - 1L)

for (target_chr in all_chroms) {
  message("\n=== ", target_chr, " ===")

  ## Stage 1: fast scan, no permutations.
  fast_file <- file.path(output_dir, paste0("AutoLD_fast_", target_chr, ".rds"))
  if (file.exists(fast_file)) {
    fast <- readRDS(fast_file)
  } else {
    fast <- AutoLD_GenomeScan(
      chr_name = target_chr,
      geno_df = Arenosa_sparse$geno,
      info_df = Arenosa_sparse$info,
      window_kb = 500,
      ploidy = 4,
      type = "full",
      calc_norm = TRUE,
      calc_pvalue = FALSE,
      n_cores = n_cores,
      min_N_eff = 50
    )
    fast <- add_autold_screening(fast, r2_threshold = 0.05, ho_threshold = 0.80)
    saveRDS(fast, fast_file)
  }

  ## Candidate universe used for candidate-level permutation calibration.
  candidates <- fast[fast$ScreenClass != "Neither", c("CHR", "POS_A", "POS_B"), drop = FALSE]
  if (!nrow(candidates)) next

  ## Stage 2: run B=1000 marginal-preserving permutations only for screened pairs.
  strict_file <- file.path(output_dir, paste0("AutoLD_permuted_", target_chr, ".rds"))
  if (!file.exists(strict_file)) {
    strict <- AutoLD_GenomeScan_can(
      chr_name = target_chr,
      geno_df = Arenosa_sparse$geno,
      info_df = Arenosa_sparse$info,
      window_kb = 500,
      ploidy = 4,
      type = "full",
      calc_norm = TRUE,
      calc_pvalue = TRUE,
      B_perm = 1000,
      n_cores = n_cores,
      min_N_eff = 50,
      candidate_pairs = candidates,
      force_pvalue = TRUE
    )
    strict <- add_autold_screening(strict, r2_threshold = 0.05, ho_threshold = 0.80)
    saveRDS(strict, strict_file)
  }

  rm(fast, candidates)
  gc(verbose = FALSE)
}

## Combine candidate-level permutation results.
perm_files <- list.files(output_dir, pattern = "^AutoLD_permuted_.*\\.rds$", full.names = TRUE)
all_perm <- data.table::rbindlist(lapply(perm_files, readRDS), fill = TRUE)

## Clean HO-only set used for conservative empirical interpretation:
## screening: r2_comp < 0.05 and H_HO >= 0.80
## calibration: P_HO <= 0.05 and P_D11 >= 0.05
HO_clean <- all_perm[
  ScreenClass == "HO-only" &
  Permutation_status == "success" &
  P_HO_Omnibus <= 0.05 &
  P_D11 >= 0.05
]

saveRDS(HO_clean, file.path(output_dir, "Arenosa_HO_only_clean.rds"))
message("Retained HO-only clean pairs: ", nrow(HO_clean))
