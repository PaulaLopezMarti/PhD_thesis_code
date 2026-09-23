#!/usr/bin/env Rscript
# ==============================================================================
# Format GWAS summary statistics for PolyFUN
#
# This script:
#   1. Loads GWAS summary statistics.
#   2. Standardizes column names to match PolyFUN conventions.
#   3. Adds minor allele frequency (MAF) information from the
#      1000 Genomes EUR reference panel.
#   4. Harmonizes alleles against the reference panel using
#      bigsnpr::snp_match(): with (strand_flip = T) and without strand flipping (strand_flip = F).
#   5. Corrects effect directions when allele orientation is reversed.
#   6. Exports cleaned summary statistics to report QC matching.
#
# Two output files are generated:
#
#   <trait>_cleansumstats_withflip.tsv
#       Summary statistics after allele harmonization allowing strand flips.
#
#   <trait>_cleansumstats_NOflip.tsv
#       Summary statistics after allele harmonization without allowing strand flips.
#
# To be run in MN5.
# Usage:
#   module load R/4.3.2
#   Rscript 03_PolyFun_format_sumstats_for_munge.R <trait_id>
#
# Example:
#   Rscript 03_PolyFun_format_sumstats_for_munge.R BALL_Alyssa
#
# Original script: 
# /home/plopez/scratch/B_cell_roadmap_2/polyFUN/scripts/checks/02_Check_sumstats_flip_strands.R
# ==============================================================================

args = commandArgs(trailingOnly=TRUE)
options(scipen=999)

# ==============================================================================
# 0. Load required packages
# ==============================================================================
library("data.table")
library("bigsnpr")


# ==============================================================================
# 1. Global configuration
# ==============================================================================
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  trait <- args[1]
}

print(paste0("Selected trait: ", trait ))


# Hardcoded paths and variables
gwas_info_file = "/gpfs/scratch/bsc08/bsc008246/Data/GWAS/SummaryStatistics_AllGWAS_complete.tsv"
path_ss <- "/gpfs/scratch/bsc08/bsc008246/Data/GWAS/summary_statistics_nodup"
outputdir <- paste0("/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/polyFUN/checks/sumstats/")
af.file <- "/gpfs/projects/bsc08/shared_projects/IJC_3Dchromatin/genomes/1000genomes/1000genomes_GRCh38/EUR_AF_PID_nodup/EUR.chrALL.GRCh38.AlleleFrequency.bed"
path.plinkfiles <- "/gpfs/projects/bsc08/shared_projects/IJC_3Dchromatin/genomes/1000genomes/1000genomes_GRCh38/plink_PID_EUR_nodup/"
chr_list <- as.character(1:22)

# outputdir <- paste0("~/scratch/B_cell_roadmap_2/polyFUN/checks/sumstats/")
# path_ss <- "~/scratch/Data/GWAS/summary_statistics_nodup"

# ==============================================================================
# 2. Load GWAS data
# ==============================================================================
print("Load GWAS...")
ss.file <- list.files(path_ss, pattern = paste0("^",trait,'_'), full.names = T)
ss <- data.table::fread(ss.file, sep='\t')
ss$CHR <- as.character(ss$CHR)

# ==============================================================================
# 3. Standardizes column names to match PolyFUN conventions
# ==============================================================================
# Remove redundant columns and rename columns as read by munge_polyfun_sumstats
if(c("EFF") %in% names(ss)){ss$EFF <- NULL}
if(c("EFFECT") %in% names(ss)){ss$EFFECT <- NULL}
if(c("EFFECT_SIZE") %in% names(ss)){ss$EFFECT_SIZE <- NULL}
if(c("OR") %in% names(ss)){ss$OR <- NULL}
if(c("CI") %in% names(ss)){ss$CI <- NULL}
names(ss)[names(ss) == 'effect_allele'] <- 'A1'
names(ss)[names(ss) == 'reference_allele'] <- 'A2'
ss$A1 = toupper(ss$A1)
ss$A2 = toupper(ss$A2)
if(c("beta") %in% names(ss)){names(ss)[names(ss) == 'beta'] <- 'BETA'}
names(ss)[names(ss) == 'z'] <- 'Z'
names(ss)[names(ss) == 'SNP'] <- 'SNP_og'
names(ss)[names(ss) == 'pid'] <- 'SNP'

# Num cases - controls
if(c("NUM_CASE") %in% names(ss)){names(ss)[names(ss) == 'NUM_CASE'] <- 'N_cases'}
if(c("NUM_CONTROL") %in% names(ss)){names(ss)[names(ss) == 'NUM_CONTROL'] <- 'N_controls'}
if(trait == "HL_ILH"){
  ss$N_controls <-  NULL
  ss$N_cases <-  NULL
}
# Remove missing values in SE
if("se" %in% names(ss)){ 
  names(ss)[names(ss) == 'se'] <- 'SE'
  ss <- ss[!is.na(ss$SE),]
}
# provisional, remove duplicated SNP in Alyssa dataset
if(trait %in% c("BALL_Alyssa", "pedBALL_Alyssa")){
  ss <- ss[ss$SNP != "8_144398231",]
}

# ==============================================================================
# 4. Add MAF from 1000 Genomes
# ==============================================================================
print("Load MAF...")
# following mCOGS code
maf_df <- data.table::fread(af.file)
maf_df$chr <- as.character(maf_df$chr)
maf_df[,maf:=ifelse(maf>0.5,1-maf,maf)]
maf_df$pid <- paste0(maf_df$chr,'_',maf_df$pos) # positional ID

# merge with summary statistics
ss.merge <- merge(ss, maf_df, by.y=c("pid"), by.x=c("SNP"), all.x=T)
ss.merge[,c("pos", "chr")] <- NULL
names(ss.merge)[names(ss.merge) == 'maf'] <- 'MAF'

# ==============================================================================
# 5. Check consistency with bigsnpr
# ==============================================================================
print("start iterating by chromosome...")

for(chr in chr_list){
  
  print(chr)
  
  ss.f <- ss.merge[ss.merge$CHR == chr,]
  
  if("BETA" %in% names(ss.f)){ # for those datasets that have BETA column present (I know this could be optimized to avoid this if statement)
    ss.f2 <- ss.f[,c("CHR","BP","SNP","A2","A1","P", "BETA","SE","Z","MAF","SNP_og")]
    names(ss.f2) = c('chr','pos','id','a0','a1',"P", "beta","SE","Z","MAF","SNP_og")
    ss.f2$chr <- as.character(ss.f2$chr)
    
    # read bim file
    bim.file <- paste0(path.plinkfiles, "ALL.chr", chr, ".phase3.genotypes.PIDs.EUR.bim")
    bimsnps.check = data.table::fread(bim.file)
    names(bimsnps.check) = c('chr','id','GD','pos','a0','a1')
    bimsnps.check$chr <- as.character(bimsnps.check$chr)
    
    # Match both data
    # return_flip_and_rev is set to TRUE to have reported each variant
    # Run with strand flipping and without
    result = snp_match(ss.f2, bimsnps.check, return_flip_and_rev = TRUE)
    result_noflip = snp_match(ss.f2, bimsnps.check, return_flip_and_rev = TRUE, strand_flip = F)
    
    # Change sign of Z for reversed SNPs
    result[result[c('_REV_')] == TRUE,c('Z')] = -result[result[c('_REV_')] == TRUE,c('Z')]
    result_noflip[result_noflip[c('_REV_')] == TRUE,c('Z')] = -result_noflip[result_noflip[c('_REV_')] == TRUE,c('Z')]
    
    # Format column names
    cols_names <- c("SNP","SNP_og","CHR","BP","P","A2","A1","BETA","SE","Z","MAF", "_FLIP_", "_REV_")
    
    result.f <- result[,c("id.ss", "SNP_og", "chr","pos","P","a0","a1","beta","SE","Z","MAF", "_FLIP_", "_REV_" )]
    names(result.f) <- cols_names
    
    result_noflip.f <- result_noflip[,c("id.ss", "SNP_og", "chr","pos","P","a0","a1","beta","SE","Z","MAF", "_FLIP_", "_REV_" )]
    names(result_noflip.f) <- cols_names
    
    # Select columns to conserve
    ss.f.tomerge <- ss.f[,-c("SNP_og","CHR","BP","P","A2","A1","BETA","SE","Z","MAF")]
    
    
  } else { # for those datasets that do not have BETA column
    ss.f2 <- ss.f[,c("CHR","BP","SNP","A2","A1","P", "Z","MAF","SNP_og")]
    names(ss.f2) = c('chr','pos','id','a0','a1',"P", "Z","MAF","SNP_og")
    ss.f2$chr <- as.character(ss.f2$chr)
    ss.f2$beta <- 1 # we are adding temporary fake beta values because it is required by snp_match
    
    # read bim file
    bim.file <- paste0(path.plinkfiles, "ALL.chr", chr, ".phase3.genotypes.PIDs.EUR.bim")
    bimsnps.check = data.table::fread(bim.file)
    names(bimsnps.check) = c('chr','id','GD','pos','a0','a1')
    bimsnps.check$chr <- as.character(bimsnps.check$chr)
    
    # Match both data
    result = snp_match(ss.f2, bimsnps.check, return_flip_and_rev = TRUE)
    result$beta <- NULL # remove fake beta values
    
    result_noflip = snp_match(ss.f2, bimsnps.check, return_flip_and_rev = TRUE, strand_flip = F)
    result_noflip$beta <- NULL # remove fake beta values
    
    # Change sign of Z for reversed SNPs
    result[result[c('_REV_')] == TRUE,c('Z')] = -result[result[c('_REV_')] == TRUE,c('Z')]
    result_noflip[result_noflip[c('_REV_')] == TRUE,c('Z')] = -result_noflip[result_noflip[c('_REV_')] == TRUE,c('Z')]
    
    # Format column names
    cols_names <- c("SNP","SNP_og","CHR","BP","P","A2","A1","Z","MAF", "_FLIP_", "_REV_")
    
    result.f <- result[,c("id.ss", "SNP_og", "chr","pos","P","a0","a1","Z","MAF", "_FLIP_", "_REV_" )]
    names(result.f) <- cols_names
    
    result_noflip.f <- result_noflip[,c("id.ss", "SNP_og", "chr","pos","P","a0","a1","Z","MAF", "_FLIP_", "_REV_" )]
    names(result_noflip.f) <- cols_names
    
    # Select columns to conserve
    ss.f.tomerge <- ss.f[,-c("SNP_og","CHR","BP","P","A2","A1","Z","MAF")]
  }
  
  if(ncol(ss.f.tomerge) >= 2){
    
    # Merge with original summary statistics
    # We are recovering those SNPs that failed the snp_match() to consider them for the report
    result.f <- merge(ss.f.tomerge, result.f, by="SNP", all.x = T)
    result_noflip.f <- merge(ss.f.tomerge, result_noflip.f, by="SNP", all.x = T)
    
    # Identify unmatched SNPs
    missing_idx <- is.na(result.f$SNP_og) 
    missing_idx_noflip <- is.na(result_noflip.f$SNP_og)
    
    # recover original GWAS information for the unmatched SNPs
    filled_rows <- merge(result.f[missing_idx, .(SNP)], ss.f, by = "SNP", all.x = TRUE)
    filled_rows_noflip <- merge(result_noflip.f[missing_idx_noflip, .(SNP)], ss.f, by = "SNP", all.x = TRUE)
    
    # Replace NA columns with info from the summary statistics
    cols_to_fill <- setdiff(names(result.f), c("SNP","coded_af", "_FLIP_", "_REV_"))  # all except SNP
    for (col in cols_to_fill) { 
      result.f[missing_idx, (col) := filled_rows[[col]]] 
      }
    
    cols_to_fill_noflip <- setdiff(names(result_noflip.f), c("SNP","coded_af", "_FLIP_", "_REV_"))  # all except SNP
    for (col in cols_to_fill_noflip) { 
      result_noflip.f[missing_idx_noflip, (col) := filled_rows_noflip[[col]]]  
      }
  }
  
  # Merge all chromosomes in one file
  if(chr == '1'){
    result.all <- result.f
    result.all_noflip <- result_noflip.f
  } else { 
    result.all <- rbind(result.all, result.f)
    result.all_noflip <- rbind(result.all_noflip, result_noflip.f)
    }
  
}

print(paste0("Saving results at... ", paste0(outputdir, trait, '_cleansumstats_withflip.tsv')))
data.table::fwrite(result.all, file=paste0(outputdir, trait, '_cleansumstats_withflip.tsv'), sep='\t')

print(paste0("Saving results at... ", paste0(outputdir, trait, '_cleansumstats_NOflip.tsv')))
data.table::fwrite(result.all_noflip, file=paste0(outputdir, trait, '_cleansumstats_NOflip.tsv'), sep='\t')

