#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# ==============================================================================
# Summarise fine-mapping results for epitests
# ==============================================================================
#
# Original script:
# ~/scratch/B_cell_roadmap_2/susie/scripts/checks/Finemapping_evaluation/02_Finemapping_Evaluation_compare_epitest_plots_MN5.R
#
# For a given trait and fine-mapping configuration, the script:
#   1. Reads fine-mapping results.
#   2. Loads QC information for LD blocks.
#   3. Optionally restricts analyses to "promising" GWAS loci
#      (minimum p-value < 1e-6).
#   4. Summarises credible set (CS) characteristics.
#   5. Counts SNPs exceeding different PIP thresholds.
#   6. Produces summary tables stratified by cell type and epitest.
#
# ==============================================================================

######################
# Prepare commands ----
# command.file = "~/scratch/B_cell_roadmap_2/susie/jobs/finemap_eval/02_finemapping_evaluation_commands.sh"
# cat("", file = command.file)
# bcellmalignancies <- c( "MM_ICR", "CLL_ILH", "pedBALL_ICR",
#                         "DLBCL_ILH","FL_ILH", "HL_ILH", "MZL_ILH", "LYMPHOMA_UKBB",
#                         "pedBALLer_ICR", "pedBALLhyper_ICR",
#                         "BALL_Alyssa","pedBALL_Alyssa")
# for(trait in bcellmalignancies){
#   a <- list.files(paste0("~/scratch/B_cell_roadmap_2/susie/analysis/summarise_finemapping/", trait, "/"), pattern="comb3")
#   a <- grep("annot", a, value=T, invert = T)
#   a <- grep("iter", a, value=T, invert = T)
#   for(test in a){
#     test <- gsub("_finemapping.tsv", "", test)
#     test <- gsub(paste0(trait, "_"), "", test)
#     command <- paste0("/gpfs/apps/MN5/GPP/R/4.3.2/INTEL/bin/Rscript /gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/susie/scripts/checks/Finemapping_evaluation/02_Finemapping_Evaluation_compare_epitest_plots_MN5.R ", trait, " ",test, "\n")
#     cat(command,file = command.file, append = T)
#     # cat(command)
#   }
# }

######################

suppressMessages(library(tidyverse))

###############################################################################
# 1. Parse arguments and Global settings
###############################################################################

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>=1) {
  trait <- args[1]
  test <- args[2]
}

message("Parsing ", trait, " | ", test)

# Variables
prom.threshold <- 10^-6
results.file <- paste0("/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/susie/analysis/summarise_finemapping/", trait, "/" , trait,"_", test, "_finemapping.tsv")
output_dir <- paste0("/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/susie/analysis/tables_summary_finemapping/", trait)
if(!dir.exists(output_dir)) {dir.create(output_dir, recursive = T)}
  
comb = strsplit(test, "_")[[1]][1]
cv = strsplit(test, "_")[[1]][2]
comb2 = strsplit(test, "_")[[1]][3]

# LD blocks file
if(comb != "standard"){
  blocks.info.file <- paste0("/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/susie/analysis/checks/select_LDblocks/", trait,"_", cv, "_", comb, "_", comb2, "_cov90_LDblocks_selection.tsv")
} else {
  blocks.info.file <- paste0("/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/susie/analysis/checks/select_LDblocks/", trait,"_", cv, "_cov90_LDblocks_selection.tsv")
}

##################################################
## 2. Load fine-mapping results and blocks info
##################################################
print(paste0("reading ", results.file))
df.all <- data.table::fread(results.file)

message("LD blocks info: ", blocks.info.file)
blocks.info <- data.table::fread(blocks.info.file)
blocks.info$ldid <- paste0(gsub("chr", "", blocks.info$chr_block), "_", blocks.info$start_block, "_", blocks.info$end_block)
blocks.info.promising <- blocks.info[blocks.info$min_pval < prom.threshold,]

##################################################
## 3. Summarise results
##################################################
for(snps.selected in c("allSNPs", "promSNPs")){
  print(snps.selected)
  if(snps.selected == "allSNPs"){ df <- df.all} else if(snps.selected == "promSNPs") { df <- df.all[df.all$ld %in% blocks.info.promising$ldid,]}
  
  ################################################
  ## 3A. Credible-set level summaries
  ################################################
  
  df.f2 <- df[df$cs_og != -1,]
  
  ## For each credible set calculate:
  ##   - number of SNPs
  ##   - number of SNPs above several PIP thresholds
  df.c <- df.f2 %>% 
    group_by(ld, celltype, epitest, cs_og) %>% 
    summarise(n_snps = length(unique(pid)),
              n_snps_0.9 = sum(pip > 0.9, na.rm = TRUE),
              n_snps_0.75 = sum(pip > 0.75, na.rm = TRUE),
              n_snps_0.5 = sum(pip > 0.5, na.rm = TRUE),
              n_snps_0.1 = sum(pip > 0.1, na.rm = TRUE))
  
  ## Aggregate CS statistics at the COO × epitest level.
  df.c2 <- df.c %>% 
    group_by(celltype, epitest) %>% 
    summarise(mean_snps.inCS = mean(n_snps),
              n_cs = n(),
              sum_snps_0.9.inCS = sum(n_snps_0.9),
              sum_snps_0.75.inCS = sum(n_snps_0.75),
              sum_snps_0.5.inCS = sum(n_snps_0.5),
              sum_snps_0.1.inCS = sum(n_snps_0.1))
  
  ## Number of LD blocks containing at least one CS
  df.c2$n_ldblocks_withCS <- length(unique(df.f2$ld))
  
  ################################################
  ## 3B. Genome-wide SNP summaries
  ################################################
  
  ## Repeat the same calculations using all SNPs,
  ## regardless of credible set assignment
  df.c.all <- df %>% 
    group_by(ld, celltype, epitest) %>% 
    summarise(n_snps = length(unique(pid)),
              n_snps_0.9 = sum(pip > 0.9, na.rm = TRUE),
              n_snps_0.75 = sum(pip > 0.75, na.rm = TRUE),
              n_snps_0.5 = sum(pip > 0.5, na.rm = TRUE),
              n_snps_0.1 = sum(pip > 0.1, na.rm = TRUE))

  df.c2.all <- df.c.all %>% 
    group_by(celltype, epitest) %>% 
    summarise(n_ldblocks = n(),
              sum_snps_0.9 = sum(n_snps_0.9),
              sum_snps_0.75 = sum(n_snps_0.75),
              sum_snps_0.5 = sum(n_snps_0.5),
              sum_snps_0.1 = sum(n_snps_0.1))
  
  ################################################
  ## 3C. Combine summaries
  ################################################
  
  df.sum <- merge(df.c2, df.c2.all, by=c("celltype", "epitest"))
  
  ## Add metadata
  df.sum$trait_id <- trait
  df.sum$CV <- cv
  df.sum$selectedSNPS <- snps.selected
  
  ################################################
  ## 3D. Save results
  ################################################
  output.file <- paste0(output_dir,"/", trait,"_", test,"_", snps.selected, "_summarytable.tsv")
  message("Saving results in...", output.file)
  data.table::fwrite(df.sum, file=output.file, sep="\t")
}
