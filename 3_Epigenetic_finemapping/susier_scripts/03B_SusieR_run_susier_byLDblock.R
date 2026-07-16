# Original script:
# /home/plopez/scratch/B_cell_roadmap_2/susie/scripts/03_SusieR_run_susier_byblock_v2.R

#!/usr/bin/env Rscript
suppressMessages(library(argparser))
args = commandArgs(trailingOnly=T)

## Usage:
## module load R/4.3.2
## /gpfs/apps/MN5/GPP/R/4.3.2/INTEL/bin/Rscript /gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/polyFUN/scripts/06_PolyFun_runSusieR_with_polyfun.R HSC

# 1. Arguments ----
p <- arg_parser("Running SusieR + PolyFun", name="Rscript 06_PolyFun_runSusieR_with_polyfun.R")
p <- add_argument(p, arg = "--trait", help = "Trait ID")
p <- add_argument(p, arg = "--block", help = "LDblock ID")
p <- add_argument(p, arg = "--cs", help = "CS threshold")
p <- add_argument(p, arg = "--iter", help = "Max num iterations")
p <- add_argument(p, arg = "--dplots", help = "Flag specifying whether to generate diagnostic plots (observed vs. expected Z-score).",
                  flag=TRUE)
p <- add_argument(p, arg = "--lambda", help = "Flag specifying whether to compute lambda values.",
                  flag=TRUE)
p <- add_argument(p, arg = "--savedata", help = "Flag specifying whether to save LD blocks clean data.",
                  flag=TRUE)
p <- add_argument(p, arg = "--savefinemapping", help = "Flag specifying whether to save finemapping results.",
                  flag=TRUE)
p <- add_argument(p, arg = "--ppiplots", help = "Flag specifying whether to save posterior probabilities plots.",
                  flag=TRUE)


opts = parse_args(p, args)

trait = opts[["trait"]]
LDblock_selected = opts[["block"]]
cs = opts[["cs"]]
n_iter = as.numeric(opts[["iter"]])
celltype= opts[["celltype"]]
polyfun_test=opts[["test"]]
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
# trait <- "pedBALL_ICR"
# # polyfun_test <- 'partial1'
# # celltype <- "HSC"
# # wd <- paste0("~/scratch/Data/GWAS/SusieR_finemapping_GRCh38/", trait, "/")
# gwas_info_file = "~/scratch/Data/GWAS/SummaryStatistics_AllGWAS_complete.tsv"
# clean_sumstats_path = paste0("~/scratch/B_cell_roadmap_2/polyFUN/data/sumstats/", trait, "_polyFUNformattedClean.tsv")
# inputdir <- paste0("~/scratch/B_cell_roadmap_2/susie/data/", trait)
# outdir <- paste0("~/scratch/B_cell_roadmap_2/susie/results/", trait)
# iter = 100
# cs = "90"
# LDblock_selected = "LDblock177"

# diagnostic_plots <- "no"
# compute_lambda <- "no"
# save_data <- "no"
# save_results <- 'no'
# ppi_plots <- 'no'

# 1. Common variables ----
gwas_info_file = "/gpfs/scratch/bsc08/bsc008246/Data/GWAS/SummaryStatistics_AllGWAS_complete.tsv"
clean_sumstats_path = paste0("/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/polyFUN/data/sumstats/", trait, "_polyFUNformattedClean.tsv")
polyfun_dir_path = "/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/polyFUN/results/"
outdir <- paste0("/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/susie/results/", trait)
inputdir <- paste0("/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/susie/data/", trait)

# cs <- "95"
L_values_list <- c(1,5)
if(cs == "95"){ coverage.value <- 0.95} else if(cs == "90"){ coverage.value <- 0.90}

# 2. Load Data ----
## Import GWAS data
files.ss <- list.files(paste0(inputdir, "/blocks_ss"), full.names = T)
# files.ss <- grep("clean|paintor", files.ss, value=T, invert=T)
files.ss <- grep(paste0(LDblock_selected, "_"), files.ss, value = T) # select only the LD block of interest

## Gwas information
gwas_info <- data.table::fread(gwas_info_file)
study.type <- tolower(gwas_info[gwas_info$trait_id == trait,]$study_type)
ncases <- gwas_info[gwas_info$trait_id == trait,]$N_cases
ncontrols <- gwas_info[gwas_info$trait_id == trait,]$N_controls
study.size <- ncases + ncontrols
prop <- ncases/study.size

## Create output directories
dir.create(paste0(outdir, "/plots/diagnostic_plots"), showWarnings = FALSE, recursive = T)
dir.create(paste0(outdir, "/diagnostic_tables"), showWarnings = FALSE, recursive = T)
dir.create(paste0(outdir, "/warnings"), showWarnings = FALSE, recursive = T)
dir.create(paste0(outdir, "/finemapping/CV1_cov", cs,"/"), showWarnings = FALSE, recursive = T)
dir.create(paste0(outdir, "/finemapping/CV5_cov", cs,"/"), showWarnings = FALSE, recursive = T)
dir.create(paste0(inputdir, "/blocks_ss/clean/"), showWarnings = FALSE, recursive = T)
dir.create(paste0(inputdir, "/LDmatrices_clean"), showWarnings = FALSE, recursive = T)
dir.create(paste0(outdir, "/plots/finemapping/"), showWarnings = FALSE, recursive = T)
dir.create(paste0(outdir, "/plots/finemapping_cov", cs,"/"), showWarnings = FALSE, recursive = T)

# ifelse(!dir.exists(output_polyfun), dir.create(output_polyfun, recursive = T), FALSE) # create output dir in polyFun directory
# ifelse(!dir.exists("results_finemapping"), dir.create("results_finemapping", recursive = T), FALSE) # create output dir in Data/GWAS directory
# ifelse(!dir.exists(paste0("plots/comparison_SusieR.PolyFun/", outdir_name)), dir.create(paste0("plots/comparison_SusieR.PolyFun/", outdir_name), recursive = T), FALSE) # create output dir in Data/GWAS directory
# ifelse(!dir.exists("diagnostic_tables"), dir.create("diagnostic_tables", recursive = T), FALSE) # create output dir in Data/GWAS directory

clean.ss <- data.table::fread(clean_sumstats_path)
names(clean.ss)[names(clean.ss) == 'SNP'] <- 'pid'




# 3. Iterate through all LD blocks ----
print(paste0("Starting to iterate throuhg ", length(files.ss), " LDblocks..."))


for(i in 1:length(files.ss)){

  print(paste0("LDblock: ",i, " ",files.ss[i]))
  ss.merge2 <- data.table::fread(files.ss[i])

  coord_block <- strsplit(x = gsub(".*_", "", basename(files.ss[i])), split = '\\.')
  chr_block <- coord_block[[1]][1]
  start_block <- as.numeric(coord_block[[1]][2])
  end_block <- as.numeric(coord_block[[1]][3])
  ld_id <- paste0(gsub("chr","",chr_block), '_', start_block, '_', end_block)
  
  # Ld matrix
  LDblock <- gsub("_.*","",basename(files.ss[i]))
  ld.file <- grep(paste0(LDblock,'_'), list.files(paste0(inputdir, "/plink_files"), pattern='.ld', full.names = T), value = T)
  bim.file <- grep(paste0(LDblock,'_'), list.files(paste0(inputdir, "/plink_files"), pattern='.bim', full.names = T), value = T)


  
  ld = data.table::fread(ld.file,sep="\t",header=F)
  ld.mat = matrix(as.vector(data.matrix(ld)), nrow=nrow(ss.merge2), ncol=nrow(ss.merge2))
  rownames(ld.mat) <- ss.merge2$pid
  colnames(ld.mat) <- ss.merge2$pid

  ## Check allele consistency with bigsnpr ----
  ss.merge.clean <- merge(ss.merge2[,c("pid")], clean.ss, by=c("pid"))

  # filter matrix
  ld.mat.f <- ld.mat[ss.merge.clean$pid, ss.merge.clean$pid]

  if(compute_lambda == TRUE){
    print("Computing lambda value....")
    lambda.raw = estimate_s_rss(ss.merge2$z, ld.mat, n=study.size)
    lambda.f = estimate_s_rss(ss.merge.clean$Z, ld.mat.f, n=study.size)

    lambda.df <- c(LDblock, chr_block, start_block, end_block, lambda.raw, lambda.f)
    names(lambda.df) <-  c("LDblock", "chr_block", "start_block", "end_block", "lambda_raw", "lambda_clean")
    if(i == 1){
      lambda.all <- lambda.df
    } else {
      lambda.all <- rbind(lambda.all, lambda.df)
      data.table::fwrite(lambda.all, file=paste0(outdir, "/diagnostic_tables/lambda_values_tmp.tsv"), sep='\t', row.names = F, col.names = T)
      }

  }

  if(diagnostic_plots == TRUE){
    print("Start diagnostic plots...")
    # Diagnostic plot ----
    # compute alpha value ?? could be useful but it's slow for now
    condz.raw = kriging_rss(ss.merge2$z, ld.mat, n=study.size)
    condz.f = kriging_rss(ss.merge.clean$Z, ld.mat.f, n=study.size)
    p1 <- condz.raw$plot + ggtitle(paste0("Raw SNPs (", dim(ss.merge2)[1], ")"))
    p2 <- condz.f$plot + ggtitle(paste0("Clean SNPs (", dim(ss.merge.clean)[1], ")"))

    my_plot_list <- list(p1,p2)
    final_plot <- plot_grid(plotlist = my_plot_list, ncol = 2, align=T)

    pdf(paste0(outdir, "/plots/diagnostic_plots/", LDblock, "_", ld_id, "_diagnostic_plot.pdf"), width=12, height=6)
    print(final_plot)
    dev.off()
    print(paste0("Finished diagnostic plots, please check: ", paste0("plots/diagnostic_plots/", LDblock, "_", ld_id, "_diagnostic_plot.pdf")))
  }

  if(save_data == TRUE){
    # Save clean data:
    data.table::fwrite(ss.merge.clean, file=paste0(inputdir, "/blocks_ss/clean/", gsub(".tsv",".f.tsv",basename(files.ss[i]))), sep='\t')
    
    #LD matrix
    data.table::fwrite(ld.mat.f, file=paste0(inputdir,"/LDmatrices_clean/", gsub(".tsv",".ld",basename(files.ss[i]))), sep=' ', col.names = F, row.names = F)
  }

  # # RUN SUSIER
  # polyfun.file.chr <- grep(paste0("\\.", gsub("chr","",chr_block), "\\."),polyfun_files, value=T)
  # polyfun.df <- data.table::fread(polyfun.file.chr)
  # result.poly <- merge(result, polyfun.df[,c("SNP", "SNPVAR")], all.x=T, by.x="id", by.y='SNP')
  
  # if(lambda.f > 0.5){  next }
  
  # Create empty lists for storing each type of plot
  p1_list <- list()
  p2_list <- list()
  p3_list <- list()
  
  if(save_results == TRUE){
  
    for(L_value in L_values_list){
      print(L_value)
  
      # RUN SUSIER
      set.seed(103)
      print(dim(ss.merge.clean))
      print(dim(ld.mat.f))
      # Save warning for plots not converging after 100 iterations
      log_file <- paste0(paste0(outdir, "/warnings/CV", L_value, "_cov", cs,"_",gsub(".tsv","",basename(files.ss[i])), "_warnings.txt"))
      log_con <- file(log_file, open = "wt")
      # fitted_rss2 = susie_rss(z = ss.merge.clean$Z, R = ld.mat.f, n = study.size, L=L_value,
      #                         max_iter = n_iter, coverage = coverage.value )
  
      result.log <- withCallingHandlers(
        tryCatch({
          fitted_rss2 = susie_rss(z = ss.merge.clean$Z, R = ld.mat.f, n = study.size, L=L_value,
                              max_iter = n_iter, coverage = coverage.value )
          fitted_rss2
        }, error = function(e) {
          writeLines(paste("ERROR:", conditionMessage(e)), log_con)
          NULL
        }),
        warning = function(w) {
          writeLines(paste("WARNING:", conditionMessage(w)), log_con)
          invokeRestart("muffleWarning")  # suppress from console
        }
      )
      close(log_con)
      
      # check plot
      locus1 <- ss.merge.clean
      locus1$variable = as.numeric(rownames(locus1))
      locus1_susie = merge(x=locus1, y=summary(fitted_rss2)$vars, by="variable", all.x=TRUE)
      
      if(save_results == TRUE){
        if(length(unique(locus1_susie$cs)) <= 1){
          output_file_standard_susieR <- paste0(outdir, "/finemapping/CV", L_value, "_cov", cs,"/",gsub(".tsv","",basename(files.ss[i])), "_susieR_cs", cs, "_L", L_value,".tsv")
          save_locus <- locus1_susie
          save_locus$ld <- ld_id
          names(save_locus)[names(save_locus) == 'variable_prob'] <- 'pip'
          save_locus$variable <- NULL
          # save_locus <- save_locus[,c("pid", "CHR", "BP", "ld", "cs", "variable_prob")]
          # names(save_locus) <- c("chr", "pos", "ld", "cred_set", "ppi")
          data.table::fwrite(save_locus, file= output_file_standard_susieR, sep='\t')
          
          # If more than 1 CS is detected
          } else if(length(unique(locus1_susie$cs)) > 1){
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
            
            output_file_standard_susieR <- paste0(outdir, "/finemapping/CV", L_value, "_cov", cs,"/",gsub(".tsv","",basename(files.ss[i])), "_susieR_cs", cs, "_L", L_value,".tsv")
            save_locus <- locus1_susie.multi
            save_locus$ld <- ld_id
            names(save_locus)[names(save_locus) == 'variable_prob'] <- 'pip'
            save_locus$variable <- NULL
            data.table::fwrite(save_locus, file= output_file_standard_susieR, sep='\t')
          }
          
      }
  
      if(ppi_plots == TRUE){
        p1 <- ggplot(locus1_susie, aes(x=BP, y=-log10(P))) +
          geom_point(alpha=0.8, size=1.3) +
          geom_point(data=subset(locus1_susie, locus1_susie$cs!=-1), aes(color=as.factor(cs)), size=2) +
          geom_hline(yintercept=-log10(5e-8), linetype="dashed", color = "blue") +
          geom_hline(yintercept=-log10(5e-6), linetype="dashed", color = "orange") +
          theme_bw() + ggtitle(paste0(gsub(".tsv","",basename(files.ss[i]))," - CV", L_value)) +
          theme(axis.title.x=element_blank())
  
        p2 <- ggplot(locus1_susie, aes(x=BP, y=Z)) +
          geom_point(alpha=0.8, size=1.3) +
          geom_hline(yintercept=0, linetype="dashed", color = "red") +
          geom_point(data=subset(locus1_susie, locus1_susie$cs!=-1), aes(color=as.factor(cs)), size=2) +
          theme_bw() + ggtitle(paste0(LDblock," Raw - L = ", L_value)) +
          theme(axis.title.x=element_blank())
  
        p3 <- ggplot(locus1_susie, aes(x=BP, y=variable_prob)) +
          geom_point(alpha=0.8, size=1.3) + ylim(0,1) +
          geom_point(data=subset(locus1_susie, locus1_susie$cs!=-1), aes(color=as.factor(cs)), size=1.4) +
          # geom_label_repel(data=subset(locus1_susie, locus1_susie$cs!=-1), aes(label=SNP), size=4, max.overlaps = 20) +
          theme_bw() + ggtitle(paste0("L = ", L_value)) +
          theme(axis.title.x=element_blank())
        
        # Store plots
        p1_list[[length(p1_list) + 1]] <- p1
        p2_list[[length(p2_list) + 1]] <- p2
        p3_list[[length(p3_list) + 1]] <- p3
      }
    }
  
    if(ppi_plots == TRUE){
      my_plots <- c(p1_list, p2_list, p3_list)
      final_plot <- plot_grid(plotlist = my_plots, ncol = 2, align=T)
      # pdf(paste0(outdir, "/plots/finemapping/",gsub(".tsv","",basename(files.ss[i])), "_susieR_cs", cs,".pdf"), width=14, height=9)
      # Save plots
      pdf(paste0(outdir, "/plots/finemapping_cov", cs,"/", LDblock, "_", ld_id ,"_", trait, "_finemapping_susieR_cs", cs,".pdf"), width=14, height=9)
      print(final_plot)
      dev.off()
    # dev.off()
    }
  }
}



if(compute_lambda == TRUE){
  lambda.all <- rbind(lambda.all, c("LDblock", "chr_block", "start_block", "end_block", "lambda_raw", "lambda_clean"))
  data.table::fwrite(lambda.all, file=paste0(outdir, "/diagnostic_tables/lambda_values_", LDblock_selected, ".tsv"), sep='\t', row.names = F, col.names = T)
}
