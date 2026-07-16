# Original script:
# /home/plopez/scratch/B_cell_roadmap_2/susie/scripts/03.1_SusieR_run_susier_PREPARECOMMANDS.R
# Prepare commands to run SusieR ----

gwas_info_file = "~/scratch/Data/GWAS/SummaryStatistics_AllGWAS_complete.tsv"
gwas_info <- data.table::fread(gwas_info_file)
gwas_selected <- gwas_info[gwas_info$selected_cogs == 'yes' & gwas_info$trait_class == "BCELL_CANCER"]$trait_id
gwas_selected <- grep("UKBB", gwas_selected, value = T, invert = T)
gwas_selected <- c(gwas_selected, "LYMPHOMA_UKBB")

jobs_dir_mn5 <- "/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/susie/jobs/03_runSusieR/"
logs_dir_mn5 <- "/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/susie/logs/03_runSusieR/"
jobs_dir_local <- "~/scratch/B_cell_roadmap_2/susie/jobs/03_runSusieR"
logs_dir_local <- "~/scratch/B_cell_roadmap_2/susie/logs/03_runSusieR"
ifelse(!dir.exists(jobs_dir_local), dir.create(jobs_dir_local, recursive = T), FALSE)
ifelse(!dir.exists(logs_dir_local), dir.create(logs_dir_local, recursive = T), FALSE)



## BY LD BLOCK V2 ----
gwas_info_file = "~/scratch/Data/GWAS/SummaryStatistics_AllGWAS_complete.tsv"
gwas_info <- data.table::fread(gwas_info_file)
gwas_selected <- gwas_info[gwas_info$selected_cogs == 'yes' & gwas_info$trait_class == "BCELL_CANCER"]$trait_id
gwas_selected <- grep("UKBB", gwas_selected, value = T, invert = T)
gwas_selected <- c(gwas_selected, "LYMPHOMA_UKBB")

# blocks_lambda <- data.table::fread("~/scratch/B_cell_roadmap_2/susie/analysis/checks/lambda_allblocks_withoutflip.tsv")


jobs_dir_mn5 <- "/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/susie/jobs/03_runSusieR/"
logs_dir_mn5 <- "/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/susie/logs/03_runSusieR/"
jobs_dir_local <- "~/scratch/B_cell_roadmap_2/susie/jobs/03_runSusieR"
logs_dir_local <- "~/scratch/B_cell_roadmap_2/susie/logs/03_runSusieR"
ifelse(!dir.exists(jobs_dir_local), dir.create(jobs_dir_local, recursive = T), FALSE)
ifelse(!dir.exists(logs_dir_local), dir.create(logs_dir_local, recursive = T), FALSE)

cs <- "90"
iter <- "100"

for(trait in gwas_selected){
  # Create one job for each trait
  file_basename <- paste0("03_runSusieR_", trait, "_cs", cs, "_iter", iter)
  output_commands_file <- paste0(jobs_dir_local, "/", file_basename ,"_commands.sh")
  cat("",file=output_commands_file)
  
  inputdir <- paste0("~/scratch/B_cell_roadmap_2/susie/data/", trait)
  files.ss <- list.files(paste0(inputdir, "/blocks_ss"), full.names = T)
  files.ss <- grep("clean", files.ss, value = T, invert = T)
  
  for(file_name in files.ss){
    block <- gsub("_.*", "", basename(file_name))
    # if(dim(blocks_lambda[blocks_lambda$LDblock == block & blocks_lambda$traitid == trait,])[1] > 0){
    #   if(blocks_lambda[blocks_lambda$LDblock == block & blocks_lambda$traitid == trait,]$warning == "BAD"){
    #   next
    # }}
    command = paste0("/gpfs/apps/MN5/GPP/R/4.3.2/INTEL/bin/Rscript /gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/susie/scripts/03_SusieR_run_susier_byblock_v2.R ",
                     "--trait ", trait, " ",
                     "--block ", block, " ",
                     "--cs ", cs, " ",
                     "--iter ", iter, " ",
                     "--savefinemapping --ppiplots --savedata\n")
    cat(command, file=output_commands_file, append = T)
  }
  if(trait %in% c("CLL_ILH", "FL_ILH", "pedBALL_ICR", "DLBCL_ILH", "MM_ICR", "LYMPHOMA_UKBB")){
    time <- "10"
  } else if(trait %in% c("HL_ILH", "BALL_Alyssa", "pedBALL_Alyssa")){ 
    time <- "5"
  } else { time <- "5"}
  
  
  # Print slurm command ----
  slurm_command <- paste0(jobs_dir_local, "/", file_basename,".cmd")
  fileConn <- file(slurm_command, open="wt")
  
  job <- paste0("#!/bin/bash\n",
                "#SBATCH --job-name=", file_basename, "\n",
                "#SBATCH -D ", logs_dir_mn5,"\n",
                "#SBATCH --output=",logs_dir_mn5, file_basename,"_%j.out\n",
                "#SBATCH --error=",logs_dir_mn5, file_basename, "_%j.err\n",
                "#SBATCH --cpus-per-task=",18,"\n",
                "#SBATCH --ntasks=",round(length(files.ss)/6),"\n",
                "#SBATCH --time=", "2", ":00:00\n",
                "#SBATCH --qos=gp_bscls\n",
                "#SBATCH --account=bsc08\n\n",
                
                "export GREASY_LOGFILE=",logs_dir_mn5, file_basename, ".log\n\n",
                
                
                "/apps/GPP/GREASY/2.2.4.2/INTEL/IMPI/bin/greasy ",jobs_dir_mn5,  file_basename, "_commands.sh\n"
  )
  
  cat(job, file = fileConn)
  close(fileConn)
  
  cat(paste0("sbatch ", file_basename,".cmd", "\n"))
}


# for(trait in gwas_selected){
#   cat(paste0("sbatch 03_runSusieR_", trait, "_cs90_iter100.cmd\n"))
# }


## ONLY DIAGNOSTIC PLOTS ----
cs <- "90"
iter <- "100"

for(trait in gwas_selected){
  # Create one job for each trait
  file_basename <- paste0("03_runSusieR_", trait, "_diagnostic")
  output_commands_file <- paste0(jobs_dir_local, "/", file_basename ,"_commands.sh")
  cat("",file=output_commands_file)
  
  inputdir <- paste0("~/scratch/B_cell_roadmap_2/susie/data/", trait)
  files.ss <- list.files(paste0(inputdir, "/blocks_ss"), full.names = T)
  files.ss <- grep("clean", files.ss, value = T, invert = T)
  
  for(file_name in files.ss){
    block <- gsub("_.*", "", basename(file_name))
    # if(dim(blocks_lambda[blocks_lambda$LDblock == block & blocks_lambda$traitid == trait,])[1] > 0){
    #   if(blocks_lambda[blocks_lambda$LDblock == block & blocks_lambda$traitid == trait,]$warning == "BAD"){
    #   next
    # }}
    command = paste0("/gpfs/apps/MN5/GPP/R/4.3.2/INTEL/bin/Rscript /gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/susie/scripts/03_SusieR_run_susier_byblock_v2.R ",
                     "--trait ", trait, " ",
                     "--block ", block, " ",
                     "--cs ", cs, " ",
                     "--iter ", iter, " ",
                     "--dplots --lambda --savedata\n")
    cat(command, file=output_commands_file, append = T)
  }
  if(trait %in% c("CLL_ILH", "FL_ILH", "pedBALL_ICR", "DLBCL_ILH", "MM_ICR", "LYMPHOMA_UKBB")){
    time <- "48"
  } else if(trait %in% c("HL_ILH", "BALL_Alyssa", "pedBALL_Alyssa")){ 
    time <- "12"
  } else { time <- "10"}
  
  
  # Print slurm command ----
  slurm_command <- paste0(jobs_dir_local, "/", file_basename,".cmd")
  fileConn <- file(slurm_command, open="wt")
  
  job <- paste0("#!/bin/bash\n",
                "#SBATCH --job-name=", file_basename, "\n",
                "#SBATCH -D ", logs_dir_mn5,"\n",
                "#SBATCH --output=",logs_dir_mn5, file_basename,"_%j.out\n",
                "#SBATCH --error=",logs_dir_mn5, file_basename, "_%j.err\n",
                "#SBATCH --cpus-per-task=",18,"\n",
                "#SBATCH --ntasks=",round(length(files.ss)/6),"\n",
                "#SBATCH --time=", "48", ":00:00\n",
                "#SBATCH --qos=gp_bscls\n",
                "#SBATCH --account=bsc08\n\n",
                
                "export GREASY_LOGFILE=",logs_dir_mn5, file_basename, ".log\n\n",
                
                
                "/apps/GPP/GREASY/2.2.4.2/INTEL/IMPI/bin/greasy ",jobs_dir_mn5,  file_basename, "_commands.sh\n"
  )
  
  cat(job, file = fileConn)
  close(fileConn)
  
  cat(paste0("sbatch ", file_basename,".cmd", "\n"))
}











#############################################
# 2) POLYFUN + SUSIER ######################
test_job <- "testCOO"
cs <- "90"
iter <- "100"

gwas_info_file = "~/scratch/Data/GWAS/SummaryStatistics_AllGWAS_complete.tsv"
gwas_info <- data.table::fread(gwas_info_file)
gwas_selected <- gwas_info[gwas_info$selected_cogs == 'yes' & gwas_info$trait_class == "BCELL_CANCER"]$trait_id
gwas_selected <- grep("UKBB", gwas_selected, value = T, invert = T)
gwas_selected <- c(gwas_selected, "LYMPHOMA_UKBB")

jobs_dir_mn5 <- "/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/susie/jobs/03_runSusieR_polyfun/"
logs_dir_mn5 <- "/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/susie/logs/03_runSusieR_polyfun/"
jobs_dir_local <- "~/scratch/B_cell_roadmap_2/susie/jobs/03_runSusieR_polyfun"
logs_dir_local <- "~/scratch/B_cell_roadmap_2/susie/logs/03_runSusieR_polyfun"
ifelse(!dir.exists(jobs_dir_local), dir.create(jobs_dir_local, recursive = T), FALSE)
ifelse(!dir.exists(logs_dir_local), dir.create(logs_dir_local, recursive = T), FALSE)

# variables
annots_list <- list("comb1A" = c("AEnh", "APromoter", "ATACpeaks", "PoEnh", "PrEnh"),
                    "comb1B" = c("AEnh", "APromoter", "ATACpeaks", "PoEnh", "PrEnh", "CDS"),
                    "comb1C" = c("AEnh", "APromoter", "ATACpeaks", "PIRs","PoEnh", "PrEnh"),
                    "comb1D" = c("AEnh", "APromoter", "ATACpeaks", "PIRs","PoEnh", "PrEnh", "CDS"),
                    "comb1AA" = c("AEnh", "APromoter", "ATACpeaks", "PoEnh", "PrEnh", "COO"),
                    "comb1BB" = c("AEnh", "APromoter", "ATACpeaks", "PoEnh", "PrEnh", "CDS", "COO"),
                    "comb1CC" = c("AEnh", "APromoter", "ATACpeaks", "PIRs","PoEnh", "PrEnh", "COO"),
                    "comb1DD" = c("AEnh", "APromoter", "ATACpeaks", "PIRs","PoEnh", "PrEnh", "CDS", "COO"),
                    
                    "comb2A"= c("AEnhATAC", "APromoterATAC","PoEnhATAC","PrEnhATAC"),
                    "comb2B"= c("AEnhATAC", "APromoterATAC", "PoEnhATAC","PrEnhATAC", "CDS"),
                    "comb2C"= c("AEnhATAC", "APromoterATAC", "PIRs","PoEnhATAC","PrEnhATAC"),
                    "comb2D"= c("AEnhATAC", "APromoterATAC", "PIRs","PoEnhATAC","PrEnhATAC", "CDS"),
                    "comb2AA"= c("AEnhATAC", "APromoterATAC","PoEnhATAC","PrEnhATAC", "COO"),
                    "comb2BB"= c("AEnhATAC", "APromoterATAC", "PoEnhATAC","PrEnhATAC", "CDS", "COO"),
                    "comb2CC"= c("AEnhATAC", "APromoterATAC", "PIRs","PoEnhATAC","PrEnhATAC", "COO"),
                    "comb2DD"= c("AEnhATAC", "APromoterATAC", "PIRs","PoEnhATAC","PrEnhATAC", "CDS", "COO"),
                    
                    "comb3AA"= c("AEnhATAC", "APromoterATAC","PoEnhATAC","PrEnhATAC", "COO"),
                    "comb3BB"= c("AEnhATAC", "APromoterATAC", "PoEnhATAC","PrEnhATAC", "CDS", "COO"),
                    "comb3CC"= c("AEnhATAC", "APromoterATAC", "PIRs","PoEnhATAC","PrEnhATAC", "COO"),
                    "comb3DD"= c("AEnhATAC", "APromoterATAC", "PIRs","PoEnhATAC","PrEnhATAC", "CDS", "COO")
)

celltypes.sort <- c("COOannot")
tests_list <- names(annots_list)
tests_list = c("comb2CC", "comb2DD")

for(trait in gwas_selected){
  # Create one job for each trait
  file_basename <- paste0("03_runSusieR_",test_job, "_", trait, "_cs", cs, "_iter", iter)
  output_commands_file <- paste0(jobs_dir_local, "/", file_basename,"_commands.sh")
  cat("",file=output_commands_file)
  comb_list <- list.files(paste0("~/scratch/B_cell_roadmap_2/polyFUN/results/", trait), pattern="COO")
  
  for(comb.name in comb_list){
    celltype = strsplit(comb.name, "_")[[1]][3]
    test = strsplit(comb.name, "_")[[1]][4]
    if(!(test %in% tests_list)) { next}
    command = paste0("/gpfs/apps/MN5/GPP/R/4.3.2/INTEL/bin/Rscript /gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/susie/scripts/04_SusieR_run_susier_polyfun.R ",
                     "--trait ", trait, " ",
                     "--celltype ", celltype, " ",
                     "--test ", test, " ",
                     "--cs ", cs, " ",
                     "--iter ", iter, " ",
                     "--savefinemapping --ppiplots\n")
    cat(command, file=output_commands_file, append = T)
  }
  
  # Print slurm command ----
  slurm_command <- paste0(jobs_dir_local, "/", file_basename,".cmd")
  fileConn <- file(slurm_command, open="wt")
  
  job <- paste0("#!/bin/bash\n",
                "#SBATCH --job-name=", file_basename, "\n",
                "#SBATCH -D ", logs_dir_mn5,"\n",
                "#SBATCH --output=",logs_dir_mn5, file_basename,"_%j.out\n",
                "#SBATCH --error=",logs_dir_mn5, file_basename,"_%j.err\n",
                "#SBATCH --cpus-per-task=",15,"\n",
                "#SBATCH --ntasks=",length(tests_list),"\n",
                "#SBATCH --time=6:00:00\n",
                "#SBATCH --qos=gp_bscls\n",
                "#SBATCH --account=bsc08\n\n",
                
                "export GREASY_LOGFILE=",logs_dir_mn5, file_basename, ".log\n\n",
                
                "/apps/GPP/GREASY/2.2.4.2/INTEL/IMPI/bin/greasy ",jobs_dir_mn5, file_basename, "_commands.sh\n"
  )
  
  cat(job, file = fileConn)
  close(fileConn)
  
  cat(paste0("sbatch ", file_basename,".cmd\n"))
}

# test_job = "testPIRs"
# for(trait in gwas_selected){
#   # Create one job for each trait
#   file_basename <- paste0("03_runSusieR_",test_job, "_", trait, "_cs", cs, "_iter", iter)
#   cat(paste0("sbatch ", file_basename,".cmd\n"))
# }



## public data
# test_job <- "test4public"
celltypes.sort <- c("E055", "E056", "E059", "E080", "E084", "E085", "E087", "E089", "E090", "E091", "E092", "E097", "E098", "E100", "E109")
annots_list <- list("partial1" = c("AEnh", "APromoter", "ATACpeaks", "PoEnh","PoPromoter", "PrEnh","PrPromoter"),
                    "partial2" = c("AEnh", "APromoter", "PoEnh","PoPromoter", "PrEnh","PrPromoter"),
                    "partial3" = c("AEnhATAC", "APromoterATAC", "PoEnh","PoPromoter", "PrEnh","PrPromoter"),
                    "partial4" = c("AEnh", "APromoter", "PoEnh", "PrEnh","PrPromoter"),
                    "partial5" = c("AEnh", "APromoter", "ATACpeaks", "PoEnh", "PrEnh","PrPromoter"),
                    "partial6" = c("AEnh", "APromoter", "PoEnh", "PrEnh"),
                    "partial7" = c("AEnh", "APromoter", "ATACpeaks", "PoEnh", "PrEnh"),
                    "partial8" = c("AEnhATAC", "APromoterATAC", "PoEnhATAC","PoPromoterATAC", "PrEnhATAC","PrPromoterATAC"),
                    "partial12" = c("AEnh", "AEnhATAC", "APromoter", "APromoterATAC","PoEnh", "PoEnhATAC", 
                                    "PoPromoter","PoPromoterATAC", "PrEnh", "PrEnhATAC","PrPromoter", "PrPromoterATAC"),
                    "partial13" = c("AEnh", "APromoter", "ATACpeaks", "PoEnh", "PrEnh", "CDS"),
                    "partial14"= c("AEnhATAC", "APromoterATAC", "PoEnhATAC","PoPromoterATAC", "PrEnhATAC","PrPromoterATAC", "CDS")
)

# tests_list <- c("full", names(annots_list))

for(trait in gwas_selected){
  # Create one job for each trait
  file_basename <- paste0("03_runSusieR_",test_job, "_", trait, "_cs", cs, "_iter", iter)
  output_commands_file <- paste0(jobs_dir_local, "/", file_basename,"_commands.sh")
  # cat("",file=output_commands_file)
  for(test in tests_list){
    for(celltype in celltypes.sort){
      
      command = paste0("/gpfs/apps/MN5/GPP/R/4.3.2/INTEL/bin/Rscript /gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/susie/scripts/04_SusieR_run_susier_polyfun.R ",
                       "--trait ", trait, " ",
                       "--celltype ", celltype, " ",
                       "--test ", test, " ",
                       "--cs ", cs, " ",
                       "--iter ", iter, " ",
                       "--savefinemapping --ppiplots\n")
      cat(command, file=output_commands_file, append = T)
      
    }
  }
  # 
  # # Print slurm command ----
  # slurm_command <- paste0(jobs_dir_local, "/", file_basename, ".cmd")
  # fileConn <- file(slurm_command, open="wt")
  # 
  # job <- paste0("#!/bin/bash\n",
  #               "#SBATCH --job-name=", file_basename, "\n",
  #               "#SBATCH -D ", logs_dir_mn5,"\n",
  #               "#SBATCH --output=",logs_dir_mn5, file_basename,"_%j.out\n",
  #               "#SBATCH --error=",logs_dir_mn5, file_basename,"_%j.err\n",
  #               "#SBATCH --cpus-per-task=",4,"\n",
  #               "#SBATCH --ntasks=",60,"\n",
  #               "#SBATCH --time=6:00:00\n",
  #               "#SBATCH --qos=gp_bscls\n",
  #               "#SBATCH --account=bsc08\n\n",
  #               
  #               "export GREASY_LOGFILE=",logs_dir_mn5, file_basename, ".log\n\n",
  #               
  #               
  #               "/apps/GPP/GREASY/2.2.4.2/INTEL/IMPI/bin/greasy ",jobs_dir_mn5,  file_basename, "_commands.sh\n"
  # )
  # 
  # cat(job, file = fileConn)
  # close(fileConn)
  
  cat(paste0("sbatch ", file_basename,".cmd\n"))
}



## JOIN PDF FINE-MAPPING PLOTS ----
library(qpdf)
trait = "MZL_ILH"
for(trait in gwas_selected){
  pdf.list = list.files(paste0("~/scratch/B_cell_roadmap_2/susie/results/", trait, "/plots/finemapping_iter_100_cov90/"), full.names = T)
  pdf_combine(
    input = pdf.list[1:2],
    output = paste0("~/scratch/B_cell_roadmap_2/susie/results/",trait, "/plots/", trait, "cov90_allBlocks.pdf")
  )
}
