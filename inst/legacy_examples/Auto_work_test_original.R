





source("Paper_util.R")
source("Auto_genome_scan.R")



#Arenosa
Arenosa_gen <- read.table("~/下载/arenosa_dosage_0to4.txt",header=T)
Arenosa_gen_info <- read.table("~/下载/arenosa_markers_filter.txt",header=T)

Arensa_sparse <- thin_snps_by_dist(geno_df=Arenosa_gen, info_df=Arenosa_gen_info, min_dist_kb = 0.2)

all_chroms <- unique(Arensa_sparse$info$CHROM)
target_chr <- all_chroms[1]

#chr1
cat("\n>>> Stage 1: Fast scan for LD landscape...\n")
results_chr1_fast <- AutoLD_GenomeScan(chr_name = target_chr,geno_df = Arensa_sparse$geno,info_df = Arensa_sparse$info,
                                       window_kb = 500,ploidy = 4, type = "full",calc_norm = TRUE, calc_pvalue = FALSE, 
                                       n_cores = 50, min_N_eff = 50)


# 阶段 2：热点靶向显著性检验 (开启 P 值)

cat("\n>>> Stage 2: Strict inference with Permutations...\n")
results_chr1_strict <- AutoLD_GenomeScan(chr_name = target_chr,geno_df = Arensa_sparse$geno,info_df = Arensa_sparse$info,
                                         window_kb = 500,ploidy = 4,type = "full",calc_norm = TRUE,calc_pvalue = TRUE,           
                                         pvalue_trigger_r2 = 0.05,pvalue_trigger_HO = 0.80, B_perm = 1000,n_cores = 50,min_N_eff = 50)



#P_HO_for_FDR <- ifelse(
#  results_chr1_strict$Permutation_tested & results_chr1_strict$Permutation_status == "success",
#  results_chr1_strict$P_HO_Omnibus,
#  1.0  
#)

#results_chr1_strict$P_HO_FDR <- p.adjust(P_HO_for_FDR, method = "BH")


output_dir <- "./AutoLD_Results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
  cat(sprintf("Created output directory: %s\n", output_dir))
}

# 2. 获取所有唯一的染色体名称
all_chroms <- unique(Arensa_sparse$info$CHROM)
cat(sprintf("Total chromosomes to process: %d\n", length(all_chroms)))

# 3. 开启全基因组循环
for (target_chr in all_chroms) {
  cat(rep("=", 60), "\n")
  cat(sprintf(">>> Processing Chromosome: %s\n", target_chr))
  
  file_stage2 <- file.path(output_dir, sprintf("AutoLD_Stage2_Strict_%s.rds", target_chr))
  if (file.exists(file_stage2)) {
    cat(sprintf("[Skip] Stage 2 for %s already completed. Moving on.\n", target_chr))
  } else {
    cat("\n>>> Stage 2: Strict inference with Permutations...\n")
    results_chr_strict <- AutoLD_GenomeScan(
      chr_name = target_chr,
      geno_df = Arensa_sparse$geno,
      info_df = Arensa_sparse$info,
      window_kb = 500,
      ploidy = 4,
      type = "full",
      calc_norm = TRUE,
      calc_pvalue = TRUE,
      pvalue_trigger_r2 = 0.05,
      pvalue_trigger_HO = 0.80,  # 极其精准的高阶触发阈值
      B_perm = 1000,             # 1000次置换，最低P值可达 0.001
      n_cores = 60,
      min_N_eff = 50
    )
    
    # 直接在此步为成功置换的位点计算保守 FDR，省去后期拼接的麻烦
    if (nrow(results_chr_strict) > 0 && "Permutation_tested" %in% names(results_chr_strict)) {
      P_HO_for_FDR <- ifelse(
        results_chr_strict$Permutation_tested & results_chr_strict$Permutation_status == "success",
        results_chr_strict$P_HO_Omnibus,
        1.0
      )
      results_chr_strict$P_HO_FDR <- p.adjust(P_HO_for_FDR, method = "BH")
    }
    
    saveRDS(results_chr_strict, file_stage2)
    cat(sprintf("Saved Stage 2 results to: %s\n", file_stage2))
  }
  
  # 清理当前循环产生的巨大 Dataframe，并强制 R 进行深度内存垃圾回收
  rm(list = intersect(ls(), c("results_chr_fast", "results_chr_strict")))
  gc(verbose = FALSE)
  
  cat(sprintf(">>> Finished Chromosome: %s\n", target_chr))
}

cat(rep("=", 60), "\n")
cat(">>> ALL CHROMOSOMES PROCESSED SUCCESSFULLY! <<<\n")


# 读取并拼合所有的 Stage 2 结果
stage2_files <- list.files("./AutoLD_Results", pattern = "Stage2_Strict_.*\\.rds$", full.names = TRUE)

all_stage2_results <- data.table::rbindlist(lapply(stage2_files, readRDS), fill = TRUE)

top_candidates <- all_stage2_results[P_HO_FDR < 0.05]
