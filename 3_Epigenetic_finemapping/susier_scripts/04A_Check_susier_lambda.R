# ============================================================
# Check SuSiE lambda diagnostics across traits
# Original script:
# /home/plopez/scratch/B_cell_roadmap_2/susie/scripts/checks/01_Check_susier_lambda.R
# ============================================================

library(data.table)
library(tidyverse)

# Hardcoded paths and variables
gwas_selected <- c("pedBALL_ICR", "pedBALLer_ICR", "pedBALLhyper_ICR","pedBALL_Alyssa","BALL_Alyssa", "LYMPHOMA_UKBB","FL_ILH","HL_ILH","DLBCL_ILH","MZL_ILH","CLL_ILH","MM_ICR")
path.results = "~/scratch/B_cell_roadmap_2/susie/results/"
outputdir <- "~/scratch/B_cell_roadmap_2/susie/analysis/checks/"

# ------------------------------------------------------------
# 1. Loop over traits and collect SuSiE lambda diagnostics
# ------------------------------------------------------------
trait <- "DLBCL_ILH"
for(i in 1:length(gwas_selected)){
  
  trait <- gwas_selected[i]
  print(trait)
  
  path_lambdas <- paste0(path.results, trait, "/diagnostic_tables") 
  lambda.files.list <- list.files(path_lambdas, pattern = "*.tsv", full.names = T)
  
  # Read all lambda values
  for(x in 1:length(lambda.files.list)){
    lambda.file <- lambda.files.list[x]
    lambda.df1 <- data.table::fread(lambda.file)
    lambda.df1 <- lambda.df1[lambda.df1$chr != "chr_block",]
    if(x == 1){
      lambda.df <- lambda.df1
    } else{ lambda.df <- rbind(lambda.df, lambda.df1)}

  }
  
  # Diagnostics of lambda values
  # lambda.df$warning <- ifelse(lambda.df$lambda_clean > 0.5, "BAD", ".")
  lambda.df$traitid <- trait
  # lambda.df$total_LDblocks <- nrow(lambda.df)
  # lambda.df$over0.2_LDblocks <- nrow(lambda.df[lambda.df$lambda_clean >= 0.2,])
  # lambda.df$under0.2_LDblocks <- nrow(lambda.df[lambda.df$lambda_clean < 0.2,])
  # lambda.df$over0.1_LDblocks <- nrow(lambda.df[lambda.df$lambda_clean >= 0.1,])
  # lambda.df$under0.1_LDblocks <- nrow(lambda.df[lambda.df$lambda_clean < 0.1,])
  # lambda.df$n_warnings <- nrow(lambda.df[lambda.df$lambda_clean > 0.5,])
  
  # Merge all traits
  if(i == 1){ lambda <- lambda.df} else { lambda <- rbind(lambda, lambda.df)}
}

lambda$lambda_raw = as.numeric(lambda$lambda_raw)
lambda$lambda_clean = as.numeric(lambda$lambda_clean)
lambda$lambda_raw = round(lambda$lambda_raw, digits = 4)
lambda$lambda_clean = round(lambda$lambda_clean, digits = 4)
# lambda.summary <- unique(lambda[,c("traitid", "total_LDblocks","n_warnings", "over0.2_LDblocks", "under0.2_LDblocks",
#                                    "over0.1_LDblocks", "under0.1_LDblocks")])



data.table::fwrite(lambda, file=paste0(outputdir, "/lambda_allblocks_v2.tsv"), sep = "\t")
# data.table::fwrite(lambda.summary, file=paste0(outputdir, "/lambda_summary_pertrait.tsv"), sep = "\t")


# Join PDFs ----
for(trait in gwas_selected){
  pdfs <- list.files(paste0("~/scratch/B_cell_roadmap_2/susie/results/", trait, "/plots/finemapping"), full.names = T)
  pdfs <- grep("LDblock1032_13_60692242_63397672_LYMPHOMA_UKBB_finemapping_susieR_cs95.pdf|LDblock236_3_26996850_30160118_LYMPHOMA_UKBB_finemapping_susieR_cs95.pdf", pdfs,
               value = T, invert = T)
  pdfs <- grep("LDblock1119_15_54844072_58093065_CLL_ILH|LDblock236_3_26996850_30160118_CLL_ILH_|LDblock484_5_152257892_156029821_CLL_ILH|LDblock757_9_19882539_22992379_CLL_ILH|LDblock877_11_1221461_3860500_CLL_ILH", pdfs,
               value = T, invert = T)
  qpdf::pdf_combine(
  input = pdfs,
  output = paste0("~/scratch/B_cell_roadmap_2/susie/results/", trait, "/plots/Allblocks_finemapping_CV1_CV5.pdf")
)
}

# Join PDFs ----
iter <- "100"
cov <- "95"
for(trait in gwas_selected[10]){
  print(trait)
  pdfs <- list.files(paste0("~/scratch/B_cell_roadmap_2/susie/results/", trait, "/plots/finemapping_iter_", iter, "_cov", cov), full.names = T)
  # pdfs <- grep("LDblock1032_13_60692242_63397672_LYMPHOMA_UKBB_finemapping_susieR_cs95.pdf|LDblock236_3_26996850_30160118_LYMPHOMA_UKBB_finemapping_susieR_cs95.pdf", pdfs,
  #              value = T, invert = T)
  # pdfs <- grep("LDblock1119_15_54844072_58093065_CLL_ILH|LDblock236_3_26996850_30160118_CLL_ILH_|LDblock484_5_152257892_156029821_CLL_ILH|LDblock757_9_19882539_22992379_CLL_ILH|LDblock877_11_1221461_3860500_CLL_ILH", pdfs,
  #              value = T, invert = T)
  qpdf::pdf_combine(
    input = pdfs,
    output = paste0("~/scratch/B_cell_roadmap_2/susie/results/", trait, "/plots/Allblocks_finemapping_CV1_CV5_iter", iter, "_cov", cov,".pdf")
  )
}


