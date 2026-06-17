# ==============================================================================
# Helper scritp to generate commands to munge sum stats
# Original script:
# /home/plopez/scratch/B_cell_roadmap_2/polyFUN/scripts/04_PolyFun_prepare_commands_mungeSumstats.R
# ==============================================================================

## Run these commands in the terminal
## Remember to: conda activate polyfun
name_file <- "04_PolyFun_mungeSumStats"
gwas_selected <- c("pedBALL_ICR", "pedBALLer_ICR", "pedBALLhyper_ICR","pedBALL_Alyssa","BALL_Alyssa", "LYMPHOMA_UKBB","FL_ILH","HL_ILH","DLBCL_ILH","MZL_ILH","CLL_ILH","MM_ICR")

# Hardcoded paths 
gwas_info_file = "~/scratch/Data/GWAS/SummaryStatistics_AllGWAS_complete.tsv"
output_slurm_file <- paste0("~/scratch/B_cell_roadmap_2/polyFUN/jobs/", name_file, ".cmd")
path_logs <- "/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/polyFUN/logs/"
path_jobs <- "/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/polyFUN/jobs/"
output_commands_file <- paste0("~/scratch/B_cell_roadmap_2/polyFUN/jobs/", name_file, "_commands.sh")

# Read GWAS information (sample size) ----
gwas_info <- data.table::fread(gwas_info_file)

cat("",file=output_commands_file)

# 1. Write commands jobs ----
for(trait in gwas_selected){
  ncases <- gwas_info[gwas_info$trait_id == trait,]$N_cases
  ncontrols <- gwas_info[gwas_info$trait_id == trait,]$N_controls
  study.size <- ncases + ncontrols
  
  sumstats_file <- paste0("/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/polyFUN/data/sumstats/", trait, "_polyFUNformattedClean.tsv")
  output_file <- paste0("/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/polyFUN/data/sumstats/", trait, "_munged.parquet")
  
  #check if there is info about case/control size in the sum stats
  path_local <- paste0("~/scratch/B_cell_roadmap_2/polyFUN/data/sumstats/", trait, "_polyFUNformattedClean.tsv") # change this
  header_file <- read.table(path_local, head = TRUE, nrows = 1, sep = "\t")[- 1, ]
  
  if("N_cases" %in% names(header_file)){
    command_munge <- paste0("python /gpfs/scratch/bsc08/bsc008246/Programs/polyfun/munge_polyfun_sumstats.py ",
                            "--sumstats ", sumstats_file, " ",
                            # "--n ", study.size, " ",
                            "--out ", output_file, " ",
                            "--min-info 0 ",
                            "--min-maf 0.01 ",
                            "\n")
    cat(command_munge, file=output_commands_file,append = T)
  } else {
    command_munge <- paste0("python /gpfs/scratch/bsc08/bsc008246/Programs/polyfun/munge_polyfun_sumstats.py ",
                            "--sumstats ", sumstats_file, " ",
                            "--n ", study.size, " ",
                            "--out ", output_file, " ",
                            "--min-info 0 ",
                            "--min-maf 0.01 ",
                            "\n")
    cat(command_munge, file=output_commands_file,append = T)
  }
}

# Write slurm command ----
job <- paste0("#!/bin/bash\n",
              "#SBATCH --job-name='04_PolyFun_mungeSumStats'\n",
              "#SBATCH -D ", path_logs,"\n",
              "#SBATCH --output=", path_logs, name_file, "_%j.out\n",
              "#SBATCH --error=", path_logs, name_file, "_%j.err\n",
              "#SBATCH --cpus-per-task=", 5,"\n",
              "#SBATCH --ntasks=", length(gwas_selected), "\n",
              "#SBATCH --time=1:00:00\n",
              "#SBATCH --qos=gp_debug\n",
              "#SBATCH --account=bsc08\n\n",
              
              "export GREASY_LOGFILE=", path_logs, name_file, ".log\n\n",
              
              "module load R/4.3.2\n\n",
              
              "/apps/GPP/GREASY/2.2.4.2/INTEL/IMPI/bin/greasy ", path_jobs, name_file, "_commands.sh", "\n"
)

## Print slurm job
fileConn <- file(output_slurm_file, open="wt")  
cat(job, file = output_slurm_file, append=TRUE)
close(fileConn)



