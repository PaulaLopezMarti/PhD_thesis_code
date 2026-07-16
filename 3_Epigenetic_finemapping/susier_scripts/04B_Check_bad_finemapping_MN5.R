#!/usr/bin/env Rscript

# ============================================================
# Script: 04B_Check_bad_finemapping_MN5.R
# Purpose:
#   QC and summarise SuSiE fine-mapping results across LD blocks
#   for multiple GWAS traits and epigenetics models.
#
# Inputs:
#   trait (GWAS ID) and test configuration (e.g. CV5_cov90)
#
# Outputs:
#   Per-trait LD block QC tables with exclusion flags and diagnostics
#
# Original script:
# /home/plopez/scratch/B_cell_roadmap_2/susie/scripts/checks/02_Check_bad_finemapping_MN5.R
# ============================================================

args = commandArgs(trailingOnly=TRUE)

##################################################

## Usage:
## module load R/4.3.2-cairo (use cairo to print pdf)
## /gpfs/apps/MN5/GPP/R/4.3.2/INTEL/bin/Rscript /gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/susie/scripts/checks/02_Check_bad_finemapping_MN5.R HL_ILH

######################
## Prepare commands ----
# combis.selected = c("comb2CC", "comb2DD")
# command.file = "~/scratch/B_cell_roadmap_2/susie/jobs/finemap_eval/02_check_bad_finemapping_commands.sh"
# cat("", file = command.file)
# bcellmalignancies <- c( "MM_ICR", "CLL_ILH", "pedBALL_ICR",
#                         "DLBCL_ILH","FL_ILH", "HL_ILH", "MZL_ILH", "LYMPHOMA_UKBB",
#                         "pedBALLer_ICR", "pedBALLhyper_ICR",
#                         "BALL_Alyssa","pedBALL_Alyssa")
# for(trait in bcellmalignancies){
#   a <- list.files(paste0("~/scratch/B_cell_roadmap_2/susie/results/", trait, "/finemapping/"), pattern="CV5")
#   for(test in a){
#     if(strsplit(test, "_")[[1]][3] %in% combis.selected){
#       command <- paste0("/gpfs/apps/MN5/GPP/R/4.3.2/INTEL/bin/Rscript /gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/susie/scripts/checks/02_Check_bad_finemapping_MN5.R ", trait, " ",test, "\n")
#       cat(command,file = command.file, append = T)
#       # cat(command)
#     }
#   }
#   command <- paste0("/gpfs/apps/MN5/GPP/R/4.3.2/INTEL/bin/Rscript /gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/susie/scripts/checks/02_Check_bad_finemapping_MN5.R ", trait, " CV5_cov90\n")
#   cat(command,file = command.file, append = T)
# }
######################

library(tidyverse)
library(data.table)

# ============================================================
# 1. Parse input arguments & Global configuration
# ============================================================
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>=1) {
  trait <- args[1]
  test <- args[2]
}

lambda.threshold = 0.4
pval.prom = 5e-5
# gwas_selected <- c("pedBALL_ICR", "pedBALLer_ICR", "pedBALLhyper_ICR","pedBALL_Alyssa","BALL_Alyssa", "LYMPHOMA_UKBB","FL_ILH","HL_ILH","DLBCL_ILH","MZL_ILH","CLL_ILH","MM_ICR")

# Modify test name to save
test.save = gsub("_iter_100", "", test)
message("Processing test: ", test.save)
message("Trait: ", trait)

# Hardcoded paths
lambda.file = "/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/susie/analysis/checks/lambda_allblocks_v2.tsv"
results.dir = "/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/susie/results/"
output.dir = "/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/susie/analysis/checks/select_LDblocks/"
output.dir.alltraits = "/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/susie/analysis/checks/LDblocks_alltraits/"

# ============================================================
# 2. Read lambda values
# ============================================================
blocks.info.all <- data.table::fread(lambda.file)


# ============================================================
# 3. Select trait 
# ============================================================

# Filter by trait
blocks.info <- blocks.info.all[blocks.info.all$traitid == trait,]
blocks.info <- blocks.info[,c("LDblock","chr_block","start_block", "end_block", "lambda_raw", "lambda_clean")]

# Flag poorly behaved LD blocks based on lambda inflation
blocks.info$lambda_over_0.4 <- ifelse(blocks.info$lambda_clean > lambda.threshold, "BAD", "ok")

# Initialize QC tracking columns for CV1 and CV5 runs
blocks.info[, `:=`(
  convIBSS_CV1 = ".",
  convIBSS_CV5 = ".",
  min_pval = ".",
  nCS_CV1 = ".",
  nCS_CV5 = ".",
  pass5e5_CV1 = ".",
  pass5e5_CV5 = ".",
  sum_ppi_CV1 = ".",
  sum_ppi_CV5 = "."
)]

# output_dir <- paste0("~/scratch/B_cell_roadmap_2/susie/analysis/summarise_finemapping/", trait)

# Fine-mapping results
path_results <- paste0(results.dir, trait, "/finemapping/", test)
files_results <- list.files(path_results, full.names = T)

# Warnings path
path_warnings <- paste0(results.dir, trait, "/warnings")
files_warnings <- list.files(path_warnings, full.names = T)

# ============================================================
# 4. Iterate over LD blocks (fine-mapping results per block)
# ============================================================
for(i in 1:length(files_results)){
  
  # print(paste0(i, " out of ", length(files_results)))
  
  # ----------------------------------------------------------
  # 4.1 Parse LD block identifiers from filename and read results
  # ----------------------------------------------------------
  ldidnum <- stringr::str_split(basename(files_results[i]), '_')[[1]][1]
  ldid <- stringr::str_split(basename(files_results[i]), '_')[[1]][2]
  l_value <- gsub(".tsv", "",stringr::str_split(basename(files_results[i]), '_')[[1]][5]) 
  
  # Read fine-mapping results
  df <- data.table::fread(files_results[i])
  
  # ----------------------------------------------------------
  # 4.2 Check convergence warnings (CV1 and CV5)
  # ----------------------------------------------------------
  warning_files <- grep(paste0("_", ldidnum, "_"), files_warnings, value = T)
  cv1 <- readLines(grep("CV1_", warning_files, value=T))
  cv5 <- readLines(grep("CV5_", warning_files, value=T))
  
  ## Flag if LD block had a convergence warning
  if(length(cv1) == 0){ 
    blocks.info[blocks.info$LDblock == ldidnum,]$convIBSS_CV1 <- "ok"
  } else{
      blocks.info[blocks.info$LDblock == ldidnum,]$convIBSS_CV1 <- "BAD" }
  if(length(cv5) == 0){ 
    blocks.info[blocks.info$LDblock == ldidnum,]$convIBSS_CV5 <- "ok"
  } else{
      blocks.info[blocks.info$LDblock == ldidnum,]$convIBSS_CV5 <- "BAD" }
  
  # ----------------------------------------------------------
  # 4.3 Load CV1 counterpart of current LD block
  # ----------------------------------------------------------
  file_selected.cv1 <- gsub("CV5", "CV1", files_results[i])
  file_selected.cv1 <- gsub("L5", "L1", file_selected.cv1)
  df.cv1 <- data.table::fread(file_selected.cv1)
  
  # Store PIP summaries
  blocks.info[blocks.info$LDblock == ldidnum,]$sum_ppi_CV1 <- sum(df.cv1$pip)
  blocks.info[blocks.info$LDblock == ldidnum,]$sum_ppi_CV5 <- sum(df$pip)
  
  # Store minimum P-value
  blocks.info[blocks.info$LDblock == ldidnum,]$min_pval <- min(df$P)
  
  
  # ============================================================
  # 5. CS-level QC checks (CV1)
  # ============================================================
  df.cv1.inCS <- df.cv1[df.cv1$cs != -1,]
  
  blocks.info[blocks.info$LDblock == ldidnum,]$nCS_CV1 <- length(unique(df.cv1.inCS$cs))
  
  # if there is any identified CS, then check if they contain at least one promising SNP
  if(dim(df.cv1.inCS)[1] != 0){
    
    prom4.all <- c()
    
    for(cs.n in unique(df.cv1.inCS$cs)){
      df.f <- df.cv1.inCS[df.cv1.inCS$cs == cs.n,]
      prom4 <- any(df.f$P < pval.prom) 
      prom4.all <- c(prom4.all, prom4)
      
    }
    
    if(any(prom4.all == F)){ 
      blocks.info[blocks.info$LDblock == ldidnum,]$pass5e5_CV1 <- "no"
    } else { 
        blocks.info[blocks.info$LDblock == ldidnum,]$pass5e5_CV1 <- "yes" }
    } 
  
  # ============================================================
  # 6. CS-level QC checks (CV5)
  # ============================================================
  df.inCS <- df[df$cs != -1,]
  
  blocks.info[blocks.info$LDblock == ldidnum,]$nCS_CV5 <- length(unique(df.inCS$cs))
  
  # if there is any identified CS, then check if they contain at least one promising SNP
  if(dim(df.inCS)[1] != 0){
    
    prom4.all <- c()
    
    for(cs.n in unique(df.inCS$cs)){
      df.f <- df.inCS[df.inCS$cs == cs.n,]
      prom4 <- any(df.f$P < pval.prom) 
      prom4.all <- c(prom4.all, prom4)
      
    }
    if(any(prom4.all == F)){ 
      blocks.info[blocks.info$LDblock == ldidnum,]$pass5e5_CV5 <- "no"
    } else {
        blocks.info[blocks.info$LDblock == ldidnum,]$pass5e5_CV5 <- "yes" }
  } 
  

  # ============================================================
  # 7. Define LD block exclusion rules
  # ============================================================
  blocks.info <- blocks.info %>%
    mutate(
      exclude = if_else(
        lambda_over_0.4 == "BAD" | # Exclude high lambda inflation blocks
          (convIBSS_CV1 == "BAD" & convIBSS_CV5 == "BAD") | # Exclude if both CV1 and CV5 fail convergence
          (pass5e5_CV1 == "no" & pass5e5_CV5 == "no"),  # Exclude if both CVs fail CS significance
        "BAD",
        "ok"
      ),
      
      # Decide fallback use of CV1: if CV5 did not converged or didnt contain any relevant SNP
      use_cv1 = if_else(convIBSS_CV5 == "BAD" | pass5e5_CV5 == "no", "yes", "no")
    )
  
}
  
# ============================================================
# 8. Save per-trait QC table
# ============================================================
output.file <- paste0(output.dir, trait, "_", test.save, "_LDblocks_selection.tsv")
data.table::fwrite(blocks.info, file=output.file, sep="\t")

#   # ============================================================
#   # 9. Summarise LD block QC results for each trait
#   # Aggregates LD block QC across traits
#   # and produces a global summary table stratified by p-value bins
#   # ============================================================
#   # 1. Create p-value bins for stratified QC summaries
#   df2 <- blocks.info[, pval_group := fifelse(min_pval < 5e-8, "<5e-8",
#                                     fifelse(min_pval < 10^-6, "<10^-6",
#                                             fifelse(min_pval < pval.prom, "<pval.prom", ">=pval.prom")))]
#  
#   # 2. Summary statistics per p-value group
#   df3 <- df2 %>% group_by(pval_group) %>% 
#     summarise(
#       # number of LD blocks
#       nBlocks = length(unique(LDblock)),
#       
#       # lambda inflation filter
#       lambda_over_0.4 = sum(lambda_over_0.4 == "BAD", na.rm = TRUE),
#       
#       # convergence failures
#       no_convIBSS_CV1 = sum(convIBSS_CV1 == "BAD", na.rm = TRUE),
#       no_convIBSS_CV5 = sum(convIBSS_CV5 == "BAD", na.rm = TRUE),
#       noconv_anyCV = sum(convIBSS_CV1 == "BAD" & convIBSS_CV5 == "BAD", na.rm = TRUE),
#       
#       # CS significance failures
#       no_pass5e5_CV1 = sum(pass5e5_CV1 == "no", na.rm = TRUE),
#       no_pass5e5_CV5 = sum(pass5e5_CV5 == "no", na.rm = TRUE),
#       
#       # exclusion summary
#       total_excluded =  sum(exclude == "BAD", na.rm = TRUE),
#       
#       # final retained blocks
#       final_blocks =  sum(exclude == "ok", na.rm = TRUE),
#       
#       # CV selection logic
#       use_cv1 =  sum(exclude == "ok" & use_cv1 == "yes", na.rm = TRUE)
#     )
#   
#   # 3. Summarise CS-level statistics (only retained blocks)
#   df4 <- df2 %>%
#     filter(exclude == "ok") %>% # only LD blocking passing QC
#     group_by(pval_group) %>%
#     summarise(
#       # blocks with at least one credible set (CV-dependent)
#       with_CS = sum(ifelse(use_cv1 == "yes", nCS_CV1 != 0, nCS_CV5 != 0)),
#       
#       # total number of credible sets across blocks
#       total_CS = sum(ifelse(use_cv1 == "yes", nCS_CV1, nCS_CV5))
#     )
#   
#   # 4. Merge QC summaries
#   df3 <- merge(df3, df4)
#   
#   # order by pval group
#   df3$pval_group <- factor(df3$pval_group , levels = c(  "<5e-8", "<1e-6","<pval.prom"))
#   df3 <- df3[order(df3$pval_group), ]
#   
#   # 5. Add total row (sum across all p-value groups)
#   df.total <- df3 %>% bind_rows(df3 %>% summarise( pval_group = "Total",
#                                                    across(where(is.numeric), \(x) sum(x, na.rm = TRUE)) ) )
#   df.total$trait_name <- trait
#   df.total <- df.total[, c("trait_name", setdiff(names(df.total), "trait_name"))]
#   
#   # 6. Store results across traits
#   if(trait == gwas_selected[1]){ df.total.all <- df.total} else{ df.total.all <- rbind(df.total.all, df.total)}
#   
# }
#   
#   
# # Save results all traits
# data.table::fwrite(df.total.all, paste0(output.dir.alltraits, "LDblocks_", test.save, "_finemapping_QC_summary_alltraits.tsv"), sep="\t")
  





