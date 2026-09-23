#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# ==============================================================================
# Convert SuSiE fine-mapping results into mCOGS input format
# ==============================================================================
#
# Original script:
# ~/scratch/Programs/mCOGS/my_scripts/clean/04_Convert_SusieR_ouput_for_mCOGS.R
#
# This script processes SuSiE fine-mapping results and prepares them for
# downstream mCOGS analyses.
#
# Main tasks:
#   1. Load SuSiE fine-mapping results for a GWAS trait.
#   2. Apply LD-block quality-control decisions from previous QC steps.
#   3. Select whether to use CV1 or CV5 fine-mapping results per LD block.
#   4. Extract posterior inclusion probabilities (PIPs) and credible sets.
#   5. Merge SuSiE results with ABF fine-mapping results for LD blocks
#      not suitable for SuSiE.
#   6. Generate mCOGS-compatible input files.
#
# Input:
#   - SuSiE fine-mapping results
#   - ABF fine-mapping results
#   - LD-block QC summaries
#
# Output:
#   - Harmonized fine-mapping results for evaluation
#   - mCOGS input files containing PIPs and credible sets
#
# Example:
# module load R/4.3.2-cairo
#
# /gpfs/apps/MN5/GPP/R/4.3.2-cairo/INTEL/bin/Rscript \
#   04_Convert_SusieR_ouput_for_mCOGS.R HL_ILH
#
# ==============================================================================

######################
## Prepare commands ----
# gwas_selected <- c("pedBALL_ICR", "pedBALLer_ICR", "pedBALLhyper_ICR","pedBALL_Alyssa","BALL_Alyssa", "LYMPHOMA_UKBB","FL_ILH","HL_ILH","DLBCL_ILH","MZL_ILH","CLL_ILH","MM_ICR")
# 
# for(trait in gwas_selected){
#   command <- paste0("/gpfs/apps/MN5/GPP/R/4.3.2/INTEL/bin/Rscript /gpfs/scratch/bsc08/bsc008246/Programs/mCOGS/my_scripts/clean/04_Convert_SusieR_ouput_for_mCOGS.R ", trait, "\n")
#   cat(command)
# }
######################

###############################################################################
# Parse arguments and Global settings
###############################################################################

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>=1) {
  trait <- args[1]
  selected.comb <- args[2]
}

library(tidyverse)


###############################################################################
# Function: format_results()
#
# Processes all LD blocks from one SuSiE fine-mapping run.
#
# Main tasks:
#   1. Apply LD-block QC decisions.
#   2. Choose CV1 or CV5 results where appropriate.
#   3. Reformat SuSiE outputs.
#   4. Save full fine-mapping results.
#   5. Merge SuSiE PIPs with ABF PIPs.
#   6. Export mCOGS-compatible file.
###############################################################################

format_results <- function(
    finemapping_dir, 
    trait_name, 
    files_results, 
    cv.s, 
    cell_type.s, 
    test.s 
){
  
  ###########################################################################
  # 1. Output paths
  ###########################################################################
  
  new_finemapping_dir <- gsub("_iter_100", "",finemapping_dir)
  # output_dir <- paste0("~/scratch/B_cell_roadmap_2/mCOGS/mCOGS_input/", trait_name)
  output_dir <- paste0("/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/mCOGS/mCOGS_input/", trait_name)
  ifelse(!dir.exists(output_dir), dir.create(output_dir, recursive = T), FALSE)
  output_file.finemapping = paste0(output.finemapped,"/",  trait_name, '_', new_finemapping_dir, "_finemapping.tsv")
  output_file.cogs <- paste0(output_dir, "/", trait_name, '_' ,new_finemapping_dir, "_gwas_4mCOGS.tsv")
  
  
  ###########################################################################
  # 2. Process each LD block
  ###########################################################################
  
  count <- 0
  
  for(file_num in 1:length(files_results)){
    
    ldidnum <- stringr::str_split(basename(files_results[file_num]), '_')[[1]][1]
    file_selected <- files_results[file_num]
    
    #########################################################################
    # Apply QC-based LD block selection
    #########################################################################
    
    # Blocks completely excluded from susieR fine-mapping 
    if(ldidnum %in% blocks.blacklist$LDblock){
      
      print(paste0("Don't use block ", ldidnum))
      next
      
      # Use CV1 instead of CV5 if failed QC
    } else if( cv.s == "CV5" & 
               ldidnum %in% blocks.blacklist.to1CV$LDblock ){
      
      print(paste0("Fine mapping not correct for block ", ldidnum, " - using CV1"))
      
      file_selected <- gsub("CV5", "CV1", file_selected)
      file_selected <- gsub("L5", "L1", file_selected)
      
      # Blocks with 1 or 0 CS in CV5
    } else if(!(ldidnum %in% c(blocks.multi$LDblock, blocks.useCV5$LDblock))){
      # In this category we have LDblocks that (A) have >1 CS using L=5 -> use pip of L=5
      # (B) have 1 CS using L=5 but 0 CS using L=1 -> we are using the pip of L=5
      # For all blocks that are not in this two categories, we will use the pips of L=1
      
      print(paste0("Fine mapping with only 1 CV for block ", ldidnum, " - using CV1"))
      
      file_selected <- gsub("CV5", "CV1", file_selected)
      file_selected <- gsub("L5", "L1", file_selected)
    }
    
    #########################################################################
    # Load selected SuSiE result file
    #########################################################################
    
    df <- data.table::fread(file_selected)
    df$ld_idnum <- ldidnum
    
    # I am adding this both columns track the original results from susieR 
    # and flag the SNPs that belong the the identified CS
    df$pip_og <- df$pip
    df$cs_og <- df$cs
    
    #########################################################################
    # Expand multiple credible sets
    #
    # If a block contains several credible sets, create one copy of the
    # block per CS and replace pip with the corresponding CS-specific PIP (alpha).
    # Right now we are doing this only for the selected comb and standard because
    # due to historical reasons, the other epitests don't have the CSXX_pip saved
    #########################################################################
    if(test %in% c("comb2DD", "standard")){
      if(length(unique(df$cs)) > 1){
        
        cs_detected <- setdiff(unique(df$cs), -1)
        
        for(num_cs in 1:length(cs_detected)){
          
          cs_detected.id <- cs_detected[num_cs]
          new.col <- paste0("CS",cs_detected.id, "_pip")
          
          df.f <- df
          df.f$pip <- df.f[,..new.col]
          df.f$cs <- cs_detected.id
          
          if(num_cs == 1){df.merge.cs <- df.f} else{ df.merge.cs <- rbind(df.merge.cs, df.f)}
        }
        df <- df.merge.cs
      }
    }
    #########################################################################
    # Combine all LD blocks
    #########################################################################
    
    if(count == 0){
      full.results <- df
      count <- 1
    } else { 
      full.results <- rbind(full.results, df, fill=T) }
  }
  
  ###########################################################################
  # 3. Save full fine-mapping results for fine-mapping assessment
  ###########################################################################
  
  to.save <- full.results
  to.save$celltype <- cell_type.s
  to.save$epitest <- test.s
  to.save$CV <- cv.s
  message("Saving finemppaing in: ", output_file.finemapping)
  data.table::fwrite(to.save, file=output_file.finemapping, sep="\t" )
  
  
  ###########################################################################
  # 4. Create mCOGS input
  #
  # Convert SuSiE output into:
  #   chr | pos | ld | cred_set | ppi
  #
  # All SNPs from retained LD blocks are preserved.
  ###########################################################################
  full.results_cs <- full.results
  full.results_cs <- full.results_cs[,c("pid", "CHR", "BP", "ld", "pip", "cs")]
  names(full.results_cs) <- c("pid", "chr", "pos", "ld", "ppi", "cred_set")
  full.results.f <- data.frame(full.results_cs)
  
  # Convert SuSiE "not in CS" label (-1) into CS0
  full.results.f$cred_set <- gsub("-1", "0", full.results.f$cred_set)
  full.results.f$cred_set <- paste0("pip_set_", full.results.f$cred_set)
  
  ###########################################################################
  # Add ABF results for LD blocks not included in epigenetic fine-mapping
  ###########################################################################
  abf.f <- abf[!(abf$ld_id %in% full.results_cs$ld),]
  
  abf.f <- abf.f[, c("pid", "chr", "pos", "ld_id", "ABF.PP", "ABF.PP_CS")]
  names(abf.f) <- c("pid", "chr", "pos", "ld", "ppi", "cred_set")
  abf.f$cred_set <- "only.abf"
  
  ###########################################################################
  # Merge SuSiE and ABF results
  ###########################################################################
  abf.merge.total <- rbind(abf.f, full.results.f)
  abf.merge.total$pid <- NULL
  abf.merge.total <- abf.merge.total[,  c("chr", "pos", "ld",  "cred_set","ppi")]
  
  ###########################################################################
  # Write final mCOGS input file
  ###########################################################################
  data.table::fwrite(abf.merge.total, file = output_file.cogs, sep='\t')
  message("Saving mCOGS input file in: ", output_file.cogs)
  return(to.save)
}


###########################################################################
# Variables and paths 
###########################################################################
output.finemapped <- paste0("/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/susie/analysis/summarise_finemapping/", trait)
path_results <- paste0("/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/susie/results/", trait, "/finemapping")
abf.results <- paste0("/gpfs/scratch/bsc08/bsc008246/Data/GWAS/ABF_finemapping_GRCh38/", trait, "/", trait, "_allBlocks_ABF_finemapping.tsv")
ldblocks.qc.path = "/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/susie/analysis/checks/select_LDblocks/"
coo.list <- list("MM_ICR" = "COO4", "pedBALL_ICR" = "COO1", "pedBALLer_ICR" = "COO1", "pedBALLhyper_ICR" = "COO1","BALL_Alyssa" = "COO1", "pedBALL_Alyssa" = "COO1", "CLL_ILH" = "COO2", "MZL_ILH" = "COO2", "FL_ILH" = "COO3", "HL_ILH" = "COO3", "DLBCL_ILH" = "COO3", "LYMPHOMA_UKBB" = "COO5")

# Create output directory
dir.create(output.finemapped, showWarnings = F, recursive = T)

###########################################################################
# Select the proper fine-mapping results files
###########################################################################
print(paste0("Processing ", trait))
coo.selected = coo.list[[trait]]
print(coo.selected)

# Read list of fine-mapping results files
dir_results <- list.files(path_results, full.names = F)

# Filter results of the selected epitest
# We are reading intially L=5 (and the script automatically will find the corresponding L=1 files to complement)
tofind = paste0(coo.selected, "_", selected.comb, "_")

if(selected.comb == "comb2DD"){ # process also standard fine-mapping (but only once)
  dir_results <- grep(paste0(tofind, "|CV5_cov90"), dir_results, value=T)
} else{
  dir_results <- grep(tofind, dir_results, value=T)
}

dir_results <- grep("CV5", dir_results, value=T)
print(dir_results)


###########################################################################
# Import ABF results
###########################################################################
abf <- data.table::fread(abf.results)


###########################################################################
# Iterate by results files
# This will correspond to two files:
# (1) epigenetic fine-mapping
# (2) standard fine-mapping
###########################################################################
for(n_file in 1:length(dir_results)){
  
  dir_name <- dir_results[n_file]
  print(dir_name)
  
  # Metadata from file names: CV, epitest fullname...
  cv <- stringr::str_split(dir_name, '_')[[1]][1]
  cell_type <- stringr::str_split(dir_name, '_')[[1]][2]
  test <- stringr::str_split(dir_name, '_')[[1]][3]
  
  if(is.na(test)){
    cell_type <- "standard"
    test <- "standard"
  }
  
  testname = gsub("_iter_100", "", dir_name)
  message("Testname: ", testname)
  
  # Reformat directory name for downstream mCOGS file naming
  dir_name_new <- gsub(paste0("_", cell_type, "_"),"_", dir_name)
  dir_name_new <- paste0(cell_type, "_", dir_name_new) # import to have this order in filename for running COGS
  
  
  ###########################################################################
  # Read LD blocks information
  ###########################################################################
  blocks.info <- data.table::fread(paste0(ldblocks.qc.path, trait, "_", testname, "_LDblocks_selection.tsv"))
  message("LD blocks file: ",paste0(ldblocks.qc.path, trait, "_", testname, "_LDblocks_selection.tsv"))
  
  # Blocks blacklist: we have to use ABF due to high lambda values
  blocks.blacklist <- blocks.info[blocks.info$exclude == "BAD",]
  
  # Blocks to use 1 CV because there are converging problems in CV 5
  blocks.blacklist.to1CV <- blocks.info[blocks.info$use_cv1 == "yes",]
  
  # Blocks with more than one CS
  blocks.multi <- blocks.info[blocks.info$nCS_CV5 > 1 & blocks.info$use_cv1 == "no",]
  
  # Ldblocks that have 1 CS using L=5 but not using L=1 -> use the pip of L=5
  # CV5 is preferred in these cases because it identifies a signal that would otherwise be missed.
  blocks.useCV5 <- blocks.info[blocks.info$nCS_CV5 == 1 & blocks.info$nCS_CV1 == 0,]
  # Remove blocks already excluded or forced to CV1
  blocks.useCV5 <- blocks.useCV5[!(blocks.useCV5$LDblock %in% c(blocks.blacklist$LDblock, blocks.blacklist.to1CV$LDblock))]
  
  ###########################################################################
  # Retrieve all fine-mapping result files for this analysis
  # Each file is one LD block
  ###########################################################################
  files_finemapping <- list.files(paste0(path_results, "/", dir_name), full.names = T)
  

  ###########################################################################
  # Format results
  ###########################################################################
  result.block <- format_results(finemapping_dir = dir_name_new, 
                                 trait_name = trait, 
                                 files_results = files_finemapping,
                                 cv.s = cv, 
                                 cell_type.s = cell_type, 
                                 test.s = test)
  
}


print("Finished")




