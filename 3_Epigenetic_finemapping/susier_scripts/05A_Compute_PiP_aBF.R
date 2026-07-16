#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Original script:
# /gpfs/scratch/bsc08/bsc008246/Programs/mCOGS/my_scripts/clean/03_Compute_PiP_aBF.R

##################################################

## Usage:
## module load R/4.3.2-cairo (use cairo to print pdf)
## /gpfs/apps/MN5/GPP/R/4.3.2-cairo/INTEL/bin/Rscript /gpfs/scratch/bsc08/bsc008246/Programs/mCOGS/my_scripts/clean/03_Compute_PiP_aBF.R HL_ILH

######################
## Prepare commands ----
# gwas_info_file = "~/scratch/Data/GWAS/SummaryStatistics_AllGWAS_complete.tsv"
# gwas_info <- data.table::fread(gwas_info_file)
# gwas_selected <- gwas_info[gwas_info$selected_cogs == 'yes' & gwas_info$trait_class == "BCELL_CANCER"]$trait_id
# gwas_selected <- grep("UKBB", gwas_selected, value = T, invert = T)
# gwas_selected <- c(gwas_selected, "LYMPHOMA_UKBB")
# 
# for(trait in gwas_selected){
#   command <- paste0("/gpfs/apps/MN5/GPP/R/4.3.2-cairo/INTEL/bin/Rscript /gpfs/scratch/bsc08/bsc008246/Programs/mCOGS/my_scripts/clean/03_Compute_PiP_aBF.R ", trait, "\n")
#   cat(command)
# }
######################

## 1. Arguments ----
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>=1) {
  trait <- args[1]
}


### Paula function ----
# Function to extract credible sets based on a function from Guille and Chris
extract_cs <- function(locus.df, threshold = 0.95, ppi_column = "PIP") {
  locus.df <- locus.df[order(locus.df[,..ppi_column], decreasing = T),]
  locus.df[,"cum_prob"] <- cumsum(locus.df[,..ppi_column])
  
  # From Guille and Chris function:
  w=which(cumsum(locus.df[,..ppi_column])>=threshold)[1]
  
  name_credset <- paste0(ppi_column, '_CS')
  locus.df[,name_credset] <- "no"
  if(!is.na(w)){
    locus.df[ 1:w ,][,name_credset] <- 'yes'
  }
  return(locus.df)
}


## 2. Load libraries -----
library(GenomicRanges)
library(tidyverse)
library(cowplot)
library(coloc)

ld_blocks_path  <- "/gpfs/scratch/bsc08/bsc008246/mCOGS/db/LD_blocks/LD_blocks_GRCh38_EUR.bed"
gwas_info_file = "/gpfs/scratch/bsc08/bsc008246/Data/GWAS/SummaryStatistics_AllGWAS_complete.tsv"
# ld_blocks_path  <- "~/scratch/mCOGS/db/LD_blocks/LD_blocks_GRCh38_EUR.bed"
outputdir <- "/gpfs/scratch/bsc08/bsc008246/Data/GWAS/ABF_finemapping_GRCh38"
outputdir <- paste0(outputdir, "/", trait)
ifelse(!dir.exists(outputdir), dir.create(outputdir), FALSE)
ifelse(!dir.exists(paste0(outputdir, "/plots")), dir.create(paste0(outputdir, "/plots")), FALSE)

# Read LD blocks ----
ld = data.table::fread(ld_blocks_path, sep='\t')
ld$chr <- gsub("chr","",ld$chr)
ld$start <- ld$start + 1 # 0-based bed
ld$ld_id <- paste0(ld$chr,'_', ld$start,'_',ld$stop)
ld$width <- ld$stop - ld$start
ld$ld_order <- 1:dim(ld)[1]

## remove weird chromosomes
chr_list <- as.character(1:23)
ld$chr <- gsub("X","23",ld$chr)
ld <- ld[ld$chr %in% chr_list,]
ld.gr <- makeGRangesFromDataFrame(ld, keep.extra.columns = T)

# Load GWAS data ----
gwas_info <- data.table::fread(gwas_info_file)
study.type <- tolower(gwas_info[gwas_info$trait_id == trait,]$study_type)
ncases <- gwas_info[gwas_info$trait_id == trait,]$N_cases
ncontrols <- gwas_info[gwas_info$trait_id == trait,]$N_controls
study.size <- ncases + ncontrols
prop <- ncases/study.size

# Load GWAS summary statistics 
ss.file <- list.files("/gpfs/scratch/bsc08/bsc008246/Data/GWAS/summary_statistics_nodup", pattern = paste0("^",trait,'_'), full.names = T)
ss <- data.table::fread(ss.file, sep='\t')
ss <- ss[,c("pid", "CHR", "BP", "P", "SNP")]
ss$CHR <- as.character(ss$CHR)
ss.gr <- makeGRangesFromDataFrame(ss, keep.extra.columns = T, start.field = "BP", end.field = "BP")

## Merge with LD blocks
ss.merge <- as.data.frame(mergeByOverlaps(ss.gr, ld.gr))
ss.merge <- ss.merge[,c("ss.gr.seqnames","ss.gr.start","ss.gr.end","pid","P","SNP","ld_id","ld_order")]
names(ss.merge) <- c("chr","start","end","pid","P","SNP","ld_id","ld_order")


## load MAF
maf_path <- "/gpfs/projects/bsc08/shared_projects/IJC_3Dchromatin/genomes/1000genomes/1000genomes_GRCh38/EUR_AF_PID_nodup/EUR.chrALL.GRCh38.AlleleFrequency.bed"
# maf_path <- "~/IJC/genomes/1000genomes/1000genomes_GRCh38/EUR_AF/EUR.chrALL.GRCh38.AlleleFrequency.bed"
maf_df <- data.table::fread(maf_path)
maf_df$chr <- as.character(maf_df$chr)
min.maf <- 0.01

print("Start to iterate through the LD blocks...")
pdf(paste0(outputdir,"/plots/", trait, "_promLDblocks_ABF_finemapping.pdf"), width=10, height=14)

count = 0
## Compute PPi for all the LD blocks of a given GWAS ----
for(i in 1:dim(ld)[1]){
  # LD block info
  ld_block_selected <- ld[i,]$ld_order
  chrom <- as.character(ld[i,]$chr)
  start_pos <- ld[i,]$start
  stop_pos <- ld[i,]$stop
  print(paste0("LDblock", ld_block_selected, " chr", chrom, ":", start_pos, '-', stop_pos))
  # Filter GWAS to the LD block region
  ss.f <- ss.merge[(ss.merge$chr == chrom) & (ss.merge$start >= start_pos) & (ss.merge$start <= stop_pos),]
  
  # If there are no SNPs in that block:
  if(dim(ss.f)[1] == 0){
    print(paste0("No SNPs found in LD block ", ld_block_selected))
  } else {
    # Clean MAF data
    maf.f <- maf_df[(maf_df$chr == chrom) & (maf_df$pos >= start_pos) & (maf_df$pos <= stop_pos),]
    maf.f[,maf:=ifelse(maf>0.5,1-maf,maf)]
    maf.f <- maf.f[maf>min.maf,]
    maf.f$pid <- paste0(maf.f$chr,'_',maf.f$pos)
    
    ss.merge2 <- merge(ss.f, maf.f[,c("pid", "maf")], by="pid")
    if(dim(ss.merge2)[1] == 0){
      print(paste0("No matched SNPs with MAF in 1KG found in LD block ", ld_block_selected))
    } else {

      # Perform fine-mapping with ABF function from coloc package
      if(study.type == 'cc'){
        ds <- list(ss.merge2$P, ss.merge2$maf, ss.merge2$pid, study.size, study.type, prop, ss.merge2$start)
        names(ds) <- c("pvalues", "MAF", "snp", "N", "type", "s", "position")
      } else {
        ds <- list(ss.merge2$P, ss.merge2$maf, ss.merge2$pid, study.size, study.type, ss.merge2$start)
        names(ds) <- c("pvalues", "MAF", "snp", "N", "type", "position")
      }
      
      ss.abf <- finemap.abf(ds)
      
      # Note the last line in this data.frame does not correspond to a SNP, but to the null model, that no SNP is causal.
      # I will remove this column
      ss.abf <- ss.abf[!is.na(ss.abf$position),]
      ss.abf <- data.table::as.data.table(ss.abf)
      ss.abf <- extract_cs(ss.abf, ppi_column = "SNP.PP")
      
      # Save LD block
      ss.abf.tomerge <- ss.abf[, c("snp","SNP.PP","SNP.PP_CS")]
      names(ss.abf.tomerge) <- c("snp","ABF.PP","ABF.PP_CS")
      ss.merge3 <- merge(ss.merge2, ss.abf.tomerge, by.x="pid", by.y="snp")
      ss.merge3$cum_prob <- NULL
      ss.merge3$end <- NULL
      names(ss.merge3)[names(ss.merge3) == 'start'] <- 'pos'
      data.table::fwrite(ss.merge3, file=paste0(outputdir,"/", trait, "_LDblock", ld_block_selected, "_chr", chrom, ".", start_pos, '.', stop_pos  ,"_ABF_finemapping.tsv"), sep='\t')
      
      if(count == 0){
        ss.all <- ss.merge3
        count = 1
      } else { ss.all <- rbind(ss.all, ss.merge3)}
      
      # Plot LD block if there is a significant SNP
      if(any(ss.merge2$P < 10^-6)){
        p1 <- ggplot(ss.merge2, aes(start/1000000, -log10(P))) +  geom_point() + ggtitle(" ") + 
          theme_minimal() + geom_hline(yintercept=-log10(5e-8), linetype="dashed", color = "blue") +
          geom_hline(yintercept=-log10(5e-6), linetype="dashed", color = "orange") +
          scale_color_manual(values=c("black", "#bf0000")) + xlab('position (Mb)') +
          theme(legend.position = "none")  + ggtitle(paste0('LDblock', ld_block_selected, " chr", chrom, ":", start_pos, '-', stop_pos))
        p2 <- ggplot(ss.abf, aes(position/1000000, SNP.PP, color=as.factor(SNP.PP_CS))) +  
          geom_point() + 
          theme_minimal() + xlab('position (Mb)') + ylim(0, 1) +
          scale_color_manual(values=c("black", "#bf0000")) +
          theme(legend.position = "none") 
        
        my_plot_list <- list(p1,p2)
        final_plot <- plot_grid(plotlist = my_plot_list, ncol = 1)
        
        # png(paste0(outputdir,"/plots/", trait, "_LDblock", ld_block_selected, " chr", chrom, ":", start_pos, '-', stop_pos  ,"_ABF_finemapping.png"), width=1200, height=1800, res=150)

        print(final_plot)
        
      }
    }
  }
}

dev.off()

print("Finished all blocks...")
print(paste0(outputdir,"/", trait, "_allBlocks_ABF_finemapping.tsv"))
data.table::fwrite(ss.all, file=paste0(outputdir,"/", trait, "_allBlocks_ABF_finemapping.tsv"), sep='\t')


