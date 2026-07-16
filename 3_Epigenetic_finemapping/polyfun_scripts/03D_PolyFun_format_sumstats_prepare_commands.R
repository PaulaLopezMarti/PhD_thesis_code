# ==============================================================================
# Helper script to generate commands to match summary statistics to reference panel
# Original script:
# /home/plopez/scratch/B_cell_roadmap_2/polyFUN/scripts/03.1_PolyFun_format_sumstats_for_munge_PREPARECOMMANDS.R
# ==============================================================================

## Run these commands in the terminal
## Remember to: conda activate polyfun
gwas_selected <- c("pedBALL_ICR", "pedBALLer_ICR", "pedBALLhyper_ICR","pedBALL_Alyssa","BALL_Alyssa", "LYMPHOMA_UKBB","FL_ILH","HL_ILH","DLBCL_ILH","MZL_ILH","CLL_ILH","MM_ICR")

# Variables and hardcoded paths ----
path_logs <- "/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/polyFUN/logs/"
path_jobs <- "/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/polyFUN/jobs/"


##########################################################
## 1) CHECK FLIP STRANDS AND MATCHING QUALITY CONTROL ----
##########################################################
name_dir <- "CHECK_sumstats_flip"

dir.create(paste0("~/scratch/B_cell_roadmap_2/polyFUN/jobs/", name_dir), showWarnings = FALSE, recursive = T)
dir.create(paste0("~/scratch/B_cell_roadmap_2/polyFUN/logs/", name_dir), showWarnings = FALSE, recursive = T)

path_logs <- paste0("/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/polyFUN/logs/", name_dir, "/")
path_jobs <- paste0("/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/polyFUN/jobs/", name_dir, "/")


# 1. Write commands jobs ----
for(trait in gwas_selected){
  output_slurm_file <- paste0("~/scratch/B_cell_roadmap_2/polyFUN/jobs/", name_dir, "/", name_dir, "_", trait, ".cmd")
  
  output_commands_file <- paste0("~/scratch/B_cell_roadmap_2/polyFUN/jobs/", name_dir, "/", name_dir, "_", trait, "_commands.sh")
  cat("",file=output_commands_file)
  
  # Command to clean sumstats
  command_clean = paste0("/gpfs/apps/MN5/GPP/R/4.3.2/INTEL/bin/Rscript /gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/polyFUN/scripts/checks/03A_PolyFun_check_sumstats_strand_flip.R ",
                         trait, "\n")
  cat(command_clean, file=output_commands_file,append = T)
  
  
  # Write slurm command ----
  
  # Prepare slurm
  job <- paste0("#!/bin/bash\n",
                "#SBATCH --job-name=", name_dir, "_", trait, "\n",
                "#SBATCH -D ", path_logs,"\n",
                "#SBATCH --output=", path_logs, name_dir, "_", trait, "_%j.out\n",
                "#SBATCH --error=", path_logs, name_dir, "_", trait, "_%j.err\n",
                "#SBATCH --cpus-per-task=", 20,"\n",
                "#SBATCH --ntasks=", 1, "\n",
                "#SBATCH --time=1:00:00\n",
                "#SBATCH --qos=gp_bscls\n",
                "#SBATCH --account=bsc08\n\n",
                
                "export GREASY_LOGFILE=", path_logs, name_dir, "_", trait, ".log\n\n",
                
                "module load R/4.3.2\n\n",
                
                "/apps/GPP/GREASY/2.2.4.2/INTEL/IMPI/bin/greasy ", path_jobs, name_dir, "_", trait, "_commands.sh", "\n"
  )
  
  ## Print slurm job
  fileConn <- file(output_slurm_file, open="wt")  
  cat(job, file = output_slurm_file, append=TRUE)
  close(fileConn)
  
  cat(paste0("sbatch ", name_dir, "_", trait, ".cmd\n"))
}

# copy and paste these commands in a command file and run as sbatch job


##########################################################
## 2) FINAL HARMONIZATION ----
##########################################################
name_file <- "03_PolyFun_format_sumstats"
output_slurm_file <- paste0("~/scratch/B_cell_roadmap_2/polyFUN/jobs/", name_file, ".cmd")
output_commands_file <- paste0("~/scratch/B_cell_roadmap_2/polyFUN/jobs/", name_file, "_commands.sh")
cat("",file=output_commands_file)

for(trait in gwas_selected){
  # Command to clean sumstats
  command_clean = paste0("/gpfs/apps/MN5/GPP/R/4.3.2/INTEL/bin/Rscript /gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/polyFUN/scripts/03C_PolyFun_format_sumstats_for_munge.R ",
                         trait, "\n")
  cat(command_clean, file=output_commands_file,append = T)
}

# Write slurm command ----
job <- paste0("#!/bin/bash\n",
              "#SBATCH --job-name=", name_file, "\n",
              "#SBATCH -D ", path_logs,"\n",
              "#SBATCH --output=", path_logs, name_file, "_%j.out\n",
              "#SBATCH --error=", path_logs, name_file, "_%j.err\n",
              "#SBATCH --cpus-per-task=", 20,"\n",
              "#SBATCH --ntasks=", length(gwas_selected), "\n",
              "#SBATCH --time=1:00:00\n",
              "#SBATCH --qos=gp_bscls\n",
              "#SBATCH --account=bsc08\n\n",
              
              "export GREASY_LOGFILE=", path_logs, name_file, ".log\n\n",
              
              "module load R/4.3.2\n\n",
              
              "/apps/GPP/GREASY/2.2.4.2/INTEL/IMPI/bin/greasy ", path_jobs, name_file, "_commands.sh", "\n"
)


## Print slurm job
fileConn <- file(output_slurm_file, open="wt")  
cat(job, file = output_slurm_file, append=TRUE)
close(fileConn)


# copy and paste these commands in a command file and run as sbatch job






