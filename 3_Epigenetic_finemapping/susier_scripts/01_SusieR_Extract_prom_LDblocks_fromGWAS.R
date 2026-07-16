#!/usr/bin/env Rscript
# ==============================================================================
# Select LD blocks to perform fine-mapping
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
# /home/plopez/scratch/B_cell_roadmap_2/susie/scripts/01_SusieR_Extract_prom_LDblocks_fromGWAS.R
# ==============================================================================

# ==============================================================================
# 0. Load required packages
# ==============================================================================
suppressMessages(library(GenomicRanges))


# ==============================================================================
# 1. Global configuration
# ==============================================================================
args = commandArgs(trailingOnly=TRUE)

# 1. Arguments ----
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>=1) {
  trait <- args[1]
}

# MAF threshold
min.maf <- 0.01 

# P-val cutoff to select LD blocks for fine-mapping
pval.cutoff <- 5e-5

# Harcoded paths and variables ----
# MN5
ld_blocks_path  <- "/gpfs/scratch/bsc08/bsc008246/mCOGS/db/LD_blocks/LD_blocks_GRCh38_EUR.bed"
gwas_info_file = "/gpfs/scratch/bsc08/bsc008246/Data/GWAS/SummaryStatistics_AllGWAS_complete.tsv"
path_ss <- "/gpfs/scratch/bsc08/bsc008246/Data/GWAS/summary_statistics_nodup"
maf_path <- "/gpfs/projects/bsc08/shared_projects/IJC_3Dchromatin/genomes/1000genomes/1000genomes_GRCh38/EUR_AF_PID_nodup/EUR.chrALL.GRCh38.AlleleFrequency.bed"
outdir <- paste0("/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/susie/data/", trait)
plink.commands.file <- paste0("/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/susie/jobs/02_SusieR_Generate_plink/02_SusieR_Generate_plink_", trait, "_commands.sh")
plink.job.file <- paste0("/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/susie/jobs/02_SusieR_Generate_plink/02_SusieR_Generate_plink_", trait, ".cmd")

jobs_dir_mn5 <- "/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/susie/jobs/02_SusieR_Generate_plink/"
logs_dir_mn5 <- "/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/susie/logs/02_SusieR_Generate_plink/"

# Local (for debugging)
# ld_blocks_path  <- "~/scratch/mCOGS/db/LD_blocks/LD_blocks_GRCh38_EUR.bed"
# gwas_info_file = "~/scratch/Data/GWAS/SummaryStatistics_AllGWAS_complete.tsv"
# path_ss <- "~/scratch/Data/GWAS/summary_statistics_nodup"
# maf_path <- "~/IJC/genomes/1000genomes/1000genomes_GRCh38/EUR_AF_PID_nodup/EUR.chrALL.GRCh38.AlleleFrequency.bed"
# outdir <- "~/scratch/Data/GWAS/SusieR_finemapping_GRCh38/pedBALL_ICR"
# plink.commands.file <- paste0("/gpfs/scratch/bsc08/bsc008246/Data/GWAS/SusieR_finemapping_GRCh38/jobs/", trait, "_plink_commands.sh")

# MHC coordinates ----
# The original coord are chr6:28510120-33480577. I expanded them to capture the whole LD blocks containing the NCBI coords
# because the SNPs remaining in those blocks were being problematic to compute their ppi (ppi == NA):
mhcStart <- 28390030
mhcEnd <- 33935919

# Create directories ----
ifelse(!dir.exists(outdir), dir.create(outdir), FALSE)
ifelse(!dir.exists(jobs_dir_mn5), dir.create(jobs_dir_mn5), FALSE)
ifelse(!dir.exists(logs_dir_mn5), dir.create(logs_dir_mn5), FALSE)
ifelse(!dir.exists(paste0(outdir, "/SNPs_list")), dir.create(paste0(outdir, "/SNPs_list")), FALSE)
ifelse(!dir.exists(paste0(outdir, "/plink_files")), dir.create(paste0(outdir, "/plink_files")), FALSE)
ifelse(!dir.exists(paste0(outdir, "/plots")), dir.create(paste0(outdir, "/plots")), FALSE)
ifelse(!dir.exists(paste0(outdir, "/plots/diagnostic_plots")), dir.create(paste0(outdir, "/plots/diagnostic_plots")), FALSE)
ifelse(!dir.exists(paste0(outdir, "/plots/comparison_Raw.Clean")), dir.create(paste0(outdir, "/plots/comparison_Raw.Clean")), FALSE)
ifelse(!dir.exists(paste0(outdir, "/blocks_ss")), dir.create(paste0(outdir, "/blocks_ss")), FALSE)

# set working dir
setwd(outdir)

# ==============================================================================
# 2. Read LD blocks information
# ==============================================================================
ld = data.table::fread(ld_blocks_path, sep='\t')
ld$chr <- gsub("chr","",ld$chr)
ld$start <- ld$start + 1 # 0-based bed
ld$ld_id <- paste0(ld$chr,'_', ld$start,'_',ld$stop)
ld$width <- ld$stop - ld$start
ld$ld_order <- 1:dim(ld)[1]

# Remove weird chromosomes
chr_list <- as.character(1:23)
ld$chr <- gsub("X","23",ld$chr)
ld <- ld[ld$chr %in% chr_list,]
ld.gr <- makeGRangesFromDataFrame(ld, keep.extra.columns = T)

# ==============================================================================
# 3. Load GWAS metadata
# ==============================================================================
gwas_info <- data.table::fread(gwas_info_file)
study.type <- tolower(gwas_info[gwas_info$trait_id == trait,]$study_type)
ncases <- gwas_info[gwas_info$trait_id == trait,]$N_cases
ncontrols <- gwas_info[gwas_info$trait_id == trait,]$N_controls
study.size <- ncases + ncontrols
prop <- ncases/study.size

# ==============================================================================
# 4. Load GWAS summary stats
# ==============================================================================
print("Load GWAS...")
ss.file <- list.files(path_ss, pattern = paste0("^",trait,'_'), full.names = T)
ss <- data.table::fread(ss.file, sep='\t')
ss$CHR <- as.character(ss$CHR)
ss.gr <- makeGRangesFromDataFrame(ss, keep.extra.columns = T, start.field = "BP", end.field = "BP")

# Merge with LD blocks
ss.merge <- as.data.frame(mergeByOverlaps(ss.gr, ld.gr))
colnames(ss.merge)[colnames(ss.merge) == 'ss.gr.seqnames'] <- 'chr'
colnames(ss.merge)[colnames(ss.merge) == 'ss.gr.start'] <- 'start'
colnames(ss.merge)[colnames(ss.merge) == 'ss.gr.end'] <- 'end'
remove_cols <- colnames(ss.merge)[grepl("ss\\.gr|ld\\.gr", colnames(ss.merge))]
ss.merge[,remove_cols] <- NULL

# ==============================================================================
# 5. Load MAF from 1KG EUR
# ==============================================================================
print("Load MAF data...")
maf_df <- data.table::fread(maf_path)
maf_df$chr <- as.character(maf_df$chr)


# ==============================================================================
# 6. Select LD blocks for fine-mapping
# ==============================================================================
print("Start to iterate throught the LD blocks...")
count = 0
total_commands = 0
cat("", file=plink.commands.file)

# Iterate through all LD blocks ----
for(i in 1:dim(ld)[1]){ 
  
  # LD block info
  ld_block_selected <- ld[i,]$ld_order
  chrom <- as.character(ld[i,]$chr)
  start_pos <- ld[i,]$start
  stop_pos <- ld[i,]$stop
  
  # Filter GWAS to the LD block region
  ss.f <- ss.merge[(ss.merge$chr == chrom) & (ss.merge$start >= start_pos) & (ss.merge$start <= stop_pos),]
  
  # check if LD block is in the MHC region; otherwise, it will continue looping
  is.mhc <- any((chrom=='6' & dplyr::between(c(start_pos, stop_pos),mhcStart,mhcEnd)))
  
  if(!is.mhc){ 
    
    if(dim(ss.f)[1] == 0){
      
      print(paste0("No SNPs found in LD block ", ld_block_selected))
      print(paste0("LDblock", ld_block_selected, " chr", chrom, ":", start_pos, '-', stop_pos))
      
    } else {
      
      # Check if there is any SNP over the p-value threhsold
      if(any(ss.f$P < pval.cutoff)){
        
        print("Found relevant SNPs!")
        print(paste0("LDblock", ld_block_selected, " chr", chrom, ":", start_pos, '-', stop_pos))
        
        # Include MAF information
        maf.f <- maf_df[(maf_df$chr == chrom) & (maf_df$pos >= start_pos) & (maf_df$pos <= stop_pos),]
        maf.f[,maf:=ifelse(maf>0.5,1-maf,maf)]
        maf.f <- maf.f[!is.na(maf.f$maf),] # remove snps not in 1KG
        maf.f <- maf.f[maf.f$maf > 0,]
        
        # Filter out SNPs by MAF
        maf.f <- maf.f[maf > min.maf,] 
        maf.f$pid <- paste0(maf.f$chr,'_',maf.f$pos) # add positional ID to match with sum stats
        
        # Merge sums stats with MAF
        ss.merge2 <- merge(ss.f, maf.f[,c("pid", "maf")], by="pid")
        
        # If the LD block still contain promising SNP:
        if(any(ss.merge2$P < pval.cutoff)){
          
          # Write full table
          table.file <- paste0("blocks_ss/LDblock", ld_block_selected, "_chr", chrom, ".", start_pos, '.', stop_pos, ".tsv" )
          data.table::fwrite(ss.merge2, file=table.file, sep='\t')
          
          # Write SNPs list
          snps.file <- paste0("SNPs_list/LDblock", ld_block_selected, "_chr", chrom, ".", start_pos, '.', stop_pos, "_snps.txt" )
          snps.file.mn5 <- paste0(outdir, "/", snps.file)
          write.table(ss.merge2$pid, file=snps.file, col.names=F, row.names=F, quote=F)
          
          # Write plink command 
          total_commands <- total_commands + 2
          plink.file <- paste0(outdir, "/plink_files/LDblock", ld_block_selected, "_chr", chrom, ".", start_pos, '.', stop_pos, "_1KG" )
          plink_command <- paste0("/gpfs/scratch/bsc08/bsc008246/Programs/plink/plink --bfile /gpfs/projects/bsc08/shared_projects/IJC_3Dchromatin/genomes/1000genomes/1000genomes_GRCh38/plink_PID_EUR_nodup/ALL.chr",chrom,".phase3.genotypes.PIDs.EUR ",
                                  "--keep-allele-order --r square ",
                                  "--extract ", snps.file.mn5, " ",
                                  "--out ", plink.file, "\n")
          cat(plink_command, file=plink.commands.file, append = T)
          plink_command <- paste0("/gpfs/scratch/bsc08/bsc008246/Programs/plink/plink --bfile /gpfs/projects/bsc08/shared_projects/IJC_3Dchromatin/genomes/1000genomes/1000genomes_GRCh38/plink_PID_EUR_nodup/ALL.chr",chrom,".phase3.genotypes.PIDs.EUR ",
                                  "--keep-allele-order --make-bed ",
                                  "--extract ", snps.file.mn5, " ",
                                  "--out ", plink.file, "\n")
          cat(plink_command, file=plink.commands.file, append = T)
        }
      }
    }
  }
}



# ==============================================================================
# 7. Prepare SLURM JOB
# ==============================================================================
job <- paste0("#!/bin/bash\n",
              "#SBATCH --job-name=02_SusieR_Generate_plink_", trait, "\n",
              "#SBATCH -D ", logs_dir_mn5,"\n",
              "#SBATCH --output=",logs_dir_mn5, "02_SusieR_Generate_plink_", trait, "_%j.out\n",
              "#SBATCH --error=",logs_dir_mn5, "02_SusieR_Generate_plink_", trait, "_%j.err\n",
              "#SBATCH --cpus-per-task=",4,"\n",
              "#SBATCH --ntasks=",total_commands,"\n",
              "#SBATCH --time=2:00:00\n",
              "#SBATCH --qos=gp_bscls\n",
              "#SBATCH --account=bsc08\n\n",
              
              "export GREASY_LOGFILE=",logs_dir_mn5, "02_SusieR_Generate_plink_", trait, ".log\n\n",
              
              "/apps/GPP/GREASY/2.2.4.2/INTEL/IMPI/bin/greasy ", jobs_dir_mn5, "02_SusieR_Generate_plink_", trait, "_commands.sh\n"
)


cat("", file = plink.job.file)
cat(job, file = plink.job.file, append=TRUE)

system2(command = 'sbatch', arg = plink.job.file)

# Run in local:
# gwas_info_file = "~/scratch/Data/GWAS/SummaryStatistics_AllGWAS_complete.tsv"
# gwas_info <- data.table::fread(gwas_info_file)
# gwas_selected <- gwas_info[gwas_info$selected_cogs == 'yes' & gwas_info$trait_class == "BCELL_CANCER"]$trait_id
# gwas_selected <- grep("UKBB", gwas_selected, value = T, invert = T)
# gwas_selected <- c(gwas_selected, "LYMPHOMA_UKBB")
# 
# for(gwas_id in gwas_selected){
#   cat(paste0("/gpfs/apps/MN5/GPP/R/4.3.2/INTEL/bin/Rscript /gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/susie/scripts/01_SusieR_Extract_prom_LDblocks_fromGWAS.R ",
#              gwas_id, "\n"))
# }
