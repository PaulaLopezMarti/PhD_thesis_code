#!/usr/bin/env Rscript

# Original script
# /home/plopez/scratch/B_cell_roadmap_2/susie/scripts/04_SusieR_run_susier_polyfun.R

suppressMessages(library(argparser))
args = commandArgs(trailingOnly=T)

## Usage:
## module load R/4.3.2
## /gpfs/apps/MN5/GPP/R/4.3.2/INTEL/bin/Rscript /gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/polyFUN/scripts/06_PolyFun_runSusieR_with_polyfun.R HSC

# 1. Arguments ----
p <- arg_parser("Running SusieR + PolyFun", name="Rscript 06_PolyFun_runSusieR_with_polyfun.R")
p <- add_argument(p, arg = "--trait", help = "Trait ID")
p <- add_argument(p, arg = "--celltype", help = "Celltype")
p <- add_argument(p, arg = "--test", help = "PolyFun test name")
p <- add_argument(p, arg = "--cs", help = "CS threshold")
p <- add_argument(p, arg = "--iter", help = "Max num iterations")
p <- add_argument(p, arg = "--savefinemapping", help = "Flag specifying whether to save finemapping results.",
                  flag=TRUE)
p <- add_argument(p, arg = "--ppiplots", help = "Flag specifying whether to save posterior probabilities plots.",
                  flag=TRUE)


opts = parse_args(p, args)

trait = opts[["trait"]]
celltype= opts[["celltype"]]
polyfun_test=opts[["test"]]
cs = opts[["cs"]]
n_iter = as.numeric(opts[["iter"]])
diagnostic_plots = opts[["dplots"]]
compute_lambda= opts[["lambda"]]
save_data= opts[["savedata"]]
save_results=opts[["savefinemapping"]]
ppi_plots = opts[["ppiplots"]]

# diagnostic_plots <- "no"
# compute_lambda <- "yes"
# save_data <- "no"
# save_results <- 'yes'
# ppi_plots <- 'yes'


#############################################
## START HERE ###############################
library(susieR)
library(ggplot2)
library(Rfast)
library(cowplot)

# Variables to debug in local:
# trait <- "HL_ILH"
# polyfun_test <- 'partial1'
# celltype <- "HSC"
# wd <- paste0("~/scratch/Data/GWAS/SusieR_finemapping_GRCh38/", trait, "/")
# gwas_info_file = "~/scratch/Data/GWAS/SummaryStatistics_AllGWAS_complete.tsv"
# clean_sumstats_path = paste0("~/scratch/B_cell_roadmap_2/polyFUN/data/sumstats/", trait, "_polyFUNformattedClean.tsv")
# inputdir <- paste0("~/scratch/B_cell_roadmap_2/susie/data/", "HL_ILH")
# outdir <- paste0("~/scratch/B_cell_roadmap_2/susie/results/", trait)
# polyfun_dir_path = paste0("~/scratch/B_cell_roadmap_2/polyFUN/results/", trait)

# 1. Common variables ----
gwas_info_file = "/gpfs/scratch/bsc08/bsc008246/Data/GWAS/SummaryStatistics_AllGWAS_complete.tsv"
outdir <- paste0("/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/susie/results/", trait)
inputdir <- paste0("/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/susie/data/", trait)
polyfun_dir_path = paste0("/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/polyFUN/results/", trait)

if(cs == "95"){ coverage.value <- 0.95} else if(cs == "90"){ coverage.value <- 0.90}
L_values_list <- c(1,5)
outdir_name <- paste0(celltype, '_', polyfun_test, "_iter_", n_iter, "_cov", cs)

# 2. Load Data ----
## Import GWAS data
files.ss <- list.files(paste0(inputdir, "/blocks_ss/clean"), full.names = T)

## Gwas information
gwas_info <- data.table::fread(gwas_info_file)
study.type <- tolower(gwas_info[gwas_info$trait_id == trait,]$study_type)
ncases <- gwas_info[gwas_info$trait_id == trait,]$N_cases
ncontrols <- gwas_info[gwas_info$trait_id == trait,]$N_controls
study.size <- ncases + ncontrols
prop <- ncases/study.size

## Create output directories
dir.create(paste0(outdir, "/finemapping/CV1_", outdir_name), showWarnings = FALSE, recursive = T)
dir.create(paste0(outdir, "/finemapping/CV5_", outdir_name), showWarnings = FALSE, recursive = T)
dir.create(paste0(outdir, "/plots/finemapping_polyfun/"), showWarnings = FALSE, recursive = T)

# ifelse(!dir.exists(output_polyfun), dir.create(output_polyfun, recursive = T), FALSE) # create output dir in polyFun directory
# ifelse(!dir.exists("results_finemapping"), dir.create("results_finemapping", recursive = T), FALSE) # create output dir in Data/GWAS directory
# ifelse(!dir.exists(paste0("plots/comparison_SusieR.PolyFun/", outdir_name)), dir.create(paste0("plots/comparison_SusieR.PolyFun/", outdir_name), recursive = T), FALSE) # create output dir in Data/GWAS directory
# ifelse(!dir.exists("diagnostic_tables"), dir.create("diagnostic_tables", recursive = T), FALSE) # create output dir in Data/GWAS directory

## Import polyfun results
outdir_name_polyfun <- paste0(trait, '_', celltype, '_', polyfun_test)
polyfun_files_path <- paste0(polyfun_dir_path, "/", outdir_name_polyfun)
polyfun_files <- list.files(polyfun_files_path, pattern="snpvar_ridge_constrained.gz", full.names = T)


# 3. Iterate through all LD blocks ----
print(paste0("Starting to iterate throuhg ", length(files.ss), " LDblocks..."))



for(L_value in L_values_list){
  if(ppi_plots == TRUE){
    pdf(paste0(outdir, "/plots/finemapping_polyfun/Allblocks_", trait, "_", outdir_name, "_susieR_cs", cs, "_L", L_value, ".pdf"), width=14, height=9)
  }
  for(i in 1:length(files.ss)){
  
    print(paste0("LDblock: ",i, " ",files.ss[i]))
    ss.merge.clean <- data.table::fread(files.ss[i])
  
    coord_block <- strsplit(x = gsub(".*_", "", basename(files.ss[i])), split = '\\.')
    chr_block <- coord_block[[1]][1]
    start_block <- as.numeric(coord_block[[1]][2])
    end_block <- as.numeric(coord_block[[1]][3])
    ld_id <- paste0(gsub("chr","",chr_block), '_', start_block, '_', end_block)
  
    # Ld matrix
    ld.mat.f <- as.matrix(data.table::fread(paste0(inputdir,"/LDmatrices_clean/", gsub(".f.tsv",".ld",basename(files.ss[i])))))
    rownames(ld.mat.f) <- ss.merge.clean$pid
    colnames(ld.mat.f) <- ss.merge.clean$pid
    
    # Import polyfun results
    polyfun.file.chr <- grep(paste0("\\.", gsub("chr","",chr_block), "\\."), polyfun_files, value=T)
    print(polyfun_files)
    print(polyfun.file.chr)
    print(polyfun_files_path)
    polyfun.df <- data.table::fread(polyfun.file.chr)
    result.poly <- merge(ss.merge.clean, polyfun.df[,c("SNP", "SNPVAR")], all.x=T, by.x="pid", by.y='SNP')
    
    # filter
    result.poly.f <- result.poly[!is.na(result.poly$SNPVAR),]
    ld.mat.poly.f <- ld.mat.f[result.poly.f$pid, result.poly.f$pid]
    
    # Create empty lists for storing each type of plot
    p1_list <- list()
    p2_list <- list()
    p3_list <- list()
    
  
    # RUN SUSIER
    set.seed(103)
    # fitted_rss2 = susie_rss(z = result.poly.f$Z, R = ld.mat.poly.f, n = study.size, L=L_value,
    #                         prior_weights = result.poly.f$SNPVAR, refine=T)
    fitted_rss2 = susie_rss(z = result.poly.f$Z, R = ld.mat.poly.f, n = study.size, L=L_value,
                            prior_weights = result.poly.f$SNPVAR, 
                            max_iter = n_iter, coverage = coverage.value)

    # check plot
    locus1 <- result.poly.f
    locus1$variable = as.numeric(rownames(locus1))
    locus1_susie = merge(x=locus1, y=summary(fitted_rss2)$vars, by="variable", all.x=TRUE)
    
    # read susier results 
    susier <- data.table::fread(paste0(outdir, "/finemapping/CV", L_value, "_cov", cs, "/",gsub(".f.tsv","",basename(files.ss[i])), "_susieR_cs", cs, "_L", L_value,".tsv"))
    susier <- susier[,c("pid", "pip", "cs")]
    names(susier) <- c("pid", "pip_susier", "cs_susier")
    locus1_susie <- merge(locus1_susie, susier, by="pid")
    
    if(save_results == TRUE){
      if(length(unique(locus1_susie$cs)) <= 1){
        output_file_standard_susieR <- paste0(outdir, "/finemapping/CV", L_value, "_", outdir_name, "/",gsub(".f.tsv","",basename(files.ss[i])), "_susieR_cs", cs, "_L", L_value,".tsv")
        print(paste0("Saving results in ", output_file_standard_susieR))
        save_locus <- locus1_susie
        save_locus$ld <- ld_id
        names(save_locus)[names(save_locus) == 'variable_prob'] <- 'pip'
        save_locus$variable <- NULL
        data.table::fwrite(save_locus, file= output_file_standard_susieR, sep='\t')
        
        # If more than 1 CS is detected
      } else if(length(unique(locus1_susie$cs)) > 1){
        print(paste0(trait, " ", ld_id, " multifinemapping"))
        locus1_susie.multi <- locus1_susie
        cs_detected <- setdiff(unique(locus1_susie$cs), -1)
        for(cs_detected.id in cs_detected){
          ppi.snps <- fitted_rss2$alpha[cs_detected.id,]
          names(ppi.snps) <- 1:length(ppi.snps)
          ppi.snps <- data.frame(ppi.snps)
          ppi.snps$variable <- as.numeric(rownames(ppi.snps))
          names(ppi.snps) <- c(paste0("CS",cs_detected.id, "_pip"), "variable")
          locus1_susie.multi = merge(x=locus1_susie.multi, y=ppi.snps, by="variable", all.x=TRUE)
        }
        
        output_file_standard_susieR <- paste0(outdir, "/finemapping/CV", L_value, "_", outdir_name, "/",gsub(".f.tsv","",basename(files.ss[i])), "_susieR_cs", cs, "_L", L_value,".tsv")
        print(paste0("Saving results in ", output_file_standard_susieR))
        save_locus <- locus1_susie.multi
        save_locus$ld <- ld_id
        names(save_locus)[names(save_locus) == 'variable_prob'] <- 'pip'
        save_locus$variable <- NULL
        data.table::fwrite(save_locus, file= output_file_standard_susieR, sep='\t')
      }
    }
    
    if(ppi_plots == TRUE){
      # Standard susieR
      p1 <- ggplot(locus1_susie, aes(x=BP, y=-log10(P))) +
        geom_point(alpha=0.8, size=1.3) +
        geom_point(data=subset(locus1_susie, locus1_susie$cs_susier!=-1), aes(color=as.factor(cs_susier)), size=2) +
        geom_hline(yintercept=-log10(5e-8), linetype="dashed", color = "blue") +
        geom_hline(yintercept=-log10(5e-6), linetype="dashed", color = "orange") +
        theme_bw() + ggtitle(paste0("Standard susieR ", gsub(".f.tsv","",basename(files.ss[i]))," - CV", L_value)) +
        theme(axis.title.x=element_blank())

      p2 <- ggplot(locus1_susie, aes(x=BP, y=Z)) +
        geom_point(alpha=0.8, size=1.3) +
        geom_hline(yintercept=0, linetype="dashed", color = "red") +
        geom_point(data=subset(locus1_susie, locus1_susie$cs_susier!=-1), aes(color=as.factor(cs_susier)), size=2) +
        theme_bw()  + ggtitle(paste0("L = ", L_value)) +
        theme(axis.title.x=element_blank())

      p3 <- ggplot(locus1_susie, aes(x=BP, y=pip_susier)) +
        geom_point(alpha=0.8, size=1.3) + ylim(0,1) +
        geom_point(data=subset(locus1_susie, locus1_susie$cs_susier!=-1), aes(color=as.factor(cs_susier)), size=1.4) +
        # geom_label_repel(data=subset(locus1_susie, locus1_susie$cs!=-1), aes(label=SNP), size=4, max.overlaps = 20) +
        theme_bw() + ggtitle(paste0("L = ", L_value)) +
        theme(axis.title.x=element_blank())
      
      
      # Polyfun + susier
      p4 <- ggplot(locus1_susie, aes(x=BP, y=-log10(P))) +
        geom_point(alpha=0.8, size=1.3) +
        geom_point(data=subset(locus1_susie, locus1_susie$cs!=-1), aes(color=as.factor(cs)), size=2) +
        geom_hline(yintercept=-log10(5e-8), linetype="dashed", color = "blue") +
        geom_hline(yintercept=-log10(5e-6), linetype="dashed", color = "orange") +
        theme_bw() + ggtitle(paste0( "Polyfun ",outdir_name)) +
        theme(axis.title.x=element_blank())
      
      p5 <- ggplot(locus1_susie, aes(x=BP, y=Z)) +
        geom_point(alpha=0.8, size=1.3) +
        geom_hline(yintercept=0, linetype="dashed", color = "red") +
        geom_point(data=subset(locus1_susie, locus1_susie$cs!=-1), aes(color=as.factor(cs)), size=2) +
        theme_bw() + ggtitle(paste0("L = ", L_value)) +
        theme(axis.title.x=element_blank())
      
      p6 <- ggplot(locus1_susie, aes(x=BP, y=variable_prob)) +
        geom_point(alpha=0.8, size=1.3) + ylim(0,1) +
        geom_point(data=subset(locus1_susie, locus1_susie$cs!=-1), aes(color=as.factor(cs)), size=1.4) +
        # geom_label_repel(data=subset(locus1_susie, locus1_susie$cs!=-1), aes(label=SNP), size=4, max.overlaps = 20) +
        theme_bw() + ggtitle(paste0("L = ", L_value)) +
        theme(axis.title.x=element_blank())
      
      my_plot_list <- list(p1,p4,p2,p5,p3,p6)
      final_plot <- plot_grid(plotlist = my_plot_list, ncol = 2, align=T)
      print(final_plot)
    }
  }
  if(ppi_plots == TRUE){  dev.off() }
}




