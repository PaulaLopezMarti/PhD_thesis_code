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
#      bigsnpr::snp_match() and remove strand flips.
#   5. Corrects effect directions when allele orientation is reversed.
#   6. Exports cleaned summary statistics to report QC matching.
#
# Output:
#   <trait>_polyFUNformattedClean.tsv: Summary statistics after allele harmonization.
#
# To be run in MN5.
# Usage:
#   module load R/4.3.2
#   Rscript 03C_PolyFun_format_sumstats_for_munge.R <trait_id>
#
# Example:
#   Rscript 03C_PolyFun_format_sumstats_for_munge.R BALL_Alyssa
#
# Original script: 
# /gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/polyFUN/scripts/03_PolyFun_format_sumstats_for_munge.R
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
outputdir <- paste0("/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/polyFUN/data/sumstats/")
path.plinkfiles <- "/gpfs/projects/bsc08/shared_projects/IJC_3Dchromatin/genomes/1000genomes/1000genomes_GRCh38/plink_PID_EUR_nodup/"
af.file <- "/gpfs/projects/bsc08/shared_projects/IJC_3Dchromatin/genomes/1000genomes/1000genomes_GRCh38/EUR_AF_PID_nodup/EUR.chrALL.GRCh38.AlleleFrequency.bed"
snps.toflip.path <- paste0("/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/polyFUN/checks/sumstats/", trait, "_cleansumstats_withflip.tsv")
chr_list <- as.character(1:22)

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


# ==============================================================================
# 4. Add Allele Frequency from 1000 Genomes
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
# 5. Import list of SNPs that should be strand flipped and remove them
# ==============================================================================
# We are using the output from 03A script
snps.toflip <- data.table::fread(snps.toflip.path)
snps.toflip.pid <- snps.toflip[snps.toflip$`_FLIP_` == T,]$SNP
print(paste0("Total SNPs to flip discarded: ", length(snps.toflip.pid)))
ss.merge <- ss.merge[!(ss.merge$SNP %in% snps.toflip.pid),]

# ==============================================================================
# 6. Check consistency with bigsnpr
# ==============================================================================
print("start iterating by chromosome...")

for(chr in chr_list){
  
  print(chr)
  
  ss.f <- ss.merge[ss.merge$CHR == chr,]
  
  if("BETA" %in% names(ss.f)){  # for those datasets that have BETA column present (I know this could be optimized to avoid this if statement)
    
    ss.f2 <- ss.f[,c("CHR","BP","SNP","A2","A1","P", "BETA","SE","Z","MAF","SNP_og")]
    names(ss.f2) = c('chr','pos','id','a0','a1',"P", "beta","SE","Z","MAF","SNP_og")
    ss.f2$chr <- as.character(ss.f2$chr)
    
    # read bim file
    bim.file <- paste0(path.plinkfiles, "ALL.chr", chr, ".phase3.genotypes.PIDs.EUR.bim")
    bimsnps.check = data.table::fread(bim.file)
    names(bimsnps.check) = c('chr','id','GD','pos','a0','a1')
    bimsnps.check$chr <- as.character(bimsnps.check$chr)
    
    # Match both data
    # we are setting strand_flip = F because we have already remove the variants that required strand flipping
    result = snp_match(ss.f2, bimsnps.check, strand_flip = F, return_flip_and_rev = TRUE)
    
    # Change sign of Z for reversed SNPs
    result[result[c('_REV_')] == TRUE,c('Z')] = -result[result[c('_REV_')] == TRUE,c('Z')]
    
    # Format column names
    cols_names <- c("SNP","SNP_og","CHR","BP","P","A2","A1","BETA","SE","Z","MAF")
    result.f <- result[,c("id.ss", "SNP_og", "chr","pos","P","a0","a1","beta","SE","Z","MAF" )]
    names(result.f) <- cols_names
    
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
    result = snp_match(ss.f2, bimsnps.check, strand_flip = F, return_flip_and_rev = TRUE)
    result$beta <- NULL # remove fake beta values
    
    # Change sign of Z for reversed SNPs
    result[result[c('_REV_')] == TRUE,c('Z')] = -result[result[c('_REV_')] == TRUE,c('Z')]
    
    # Format column names
    cols_names <- c("SNP","SNP_og","CHR","BP","P","A2","A1","Z","MAF")
    result.f <- result[,c("id.ss", "SNP_og", "chr","pos","P","a0","a1","Z","MAF" )]
    names(result.f) <- cols_names
    
    # Select columns to conserve
    ss.f.tomerge <- ss.f[,-c("SNP_og","CHR","BP","P","A2","A1","Z","MAF")]
  }
  
  # Merge all chromosomes in one file
  if(ncol(ss.f.tomerge) >= 2){
    result.f <- merge(result.f, ss.f.tomerge, by="SNP")
  }
  
  if(chr == '1'){
    result.all <- result.f
  } else { result.all <- rbind(result.all, result.f)}
  
}

# Save output
print(paste0("Saving results at... ", paste0(outputdir, trait, '_polyFUNformattedClean.tsv')))
data.table::fwrite(result.all, file=paste0(outputdir, trait, '_polyFUNformattedClean.tsv'), sep='\t')

# Things to fix for final pipeline:
# It is quite redundant with 03A: take advantage of what it is done in that script so we avoid repeating the same analysis twice