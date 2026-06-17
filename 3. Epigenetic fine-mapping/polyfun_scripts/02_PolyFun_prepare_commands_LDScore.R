# ==============================================================================
# Generate commands to compute annotation LD-scores and regression weights
# Original script:
# /home/plopez/scratch/B_cell_roadmap_2/polyFUN/scripts/02_PolyFun_prepare_commands_LDScore.R
# ==============================================================================

# Variables and hardcoded paths ----
annot.dir <- "COOannot"
chr_list <- as.character(1:22)
output_dir_ld <- paste0("/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/polyFUN/data/annotations/", annot.dir)
output_commands_file <- paste0("/home/plopez/scratch/B_cell_roadmap_2/polyFUN/jobs/02_PolyFun_compute_LDScores_1KG_", annot.dir, "_commands.sh")

cat("", file=output_commands_file)

# Output files need to end with this extension '.l2.ldscore.parquet' in order to run polyFun
# Files need to be stored in the same directory as annotations

# Commands to compute LD-scores for annotations
for(chr in chr_list){
  plink_file <- paste0("/gpfs/projects/bsc08/shared_projects/IJC_3Dchromatin/genomes/1000genomes/1000genomes_GRCh38/plink_PID_EUR_nodup/ALL.chr", chr, ".phase3.genotypes.PIDs.EUR")
  parquet_file <- paste0("/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/polyFUN/data/annotations/", annot.dir, "/",annot.dir,".", chr, ".annot.parquet")
  output_file <- paste0(output_dir_ld, '/', annot.dir,'.',chr ,'.l2.ldscore.parquet')
  command_ldscores <- paste0("python /gpfs/scratch/bsc08/bsc008246/Programs/polyfun/compute_ldscores.py ",
                             "--bfile ", plink_file, " ",
                             "--annot ", parquet_file, " ",
                             "--out ", output_file, " ",
                             "--ld-wind-kb 1000", 
                             "\n")
  # using --ld-wind-kb 1Mb because the plink file I'm using don't have CM data (https://github.com/omerwe/polyfun/issues/103)
  cat(command_ldscores, file=output_commands_file, append = T)
}

# Commands to compute weights (without annotation parameter)
output_dir_w <- paste0("/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/polyFUN/data/annotations/", annot.dir)
for(chr in chr_list){
  plink_file <- paste0("/gpfs/projects/bsc08/shared_projects/IJC_3Dchromatin/genomes/1000genomes/1000genomes_GRCh38/plink_PID_EUR_nodup/ALL.chr", chr, ".phase3.genotypes.PIDs.EUR")
  parquet_file <- paste0("/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/polyFUN/data/annotations/", annot.dir, "/", annot.dir, ".", chr, ".annot.parquet")
  output_file <- paste0(output_dir_w, '/', annot.dir,'_w.',chr ,'.l2.ldscore.parquet')
  command_ldscores <- paste0("python /gpfs/scratch/bsc08/bsc008246/Programs/polyfun/compute_ldscores.py ",
                             "--bfile ", plink_file, " ",
                             "--out ", output_file, " ",
                             "--ld-wind-kb 1000", 
                             "\n")
  cat(command_ldscores, file=output_commands_file, append = T)
}


## Run these commands in the terminal
## Remember to: conda activate polyfun