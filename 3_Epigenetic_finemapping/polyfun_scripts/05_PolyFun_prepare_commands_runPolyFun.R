# ==============================================================================
# Helper script to generate commands to run PolyFun L2-regularized S-LDSC
# Annotations (.annot) and .l2.M have to be in the same dir as ldscores
# Original script:
# /home/plopez/scratch/B_cell_roadmap_2/polyFUN/scripts/05_PolyFun_prepare_commands_runPolyFun_COO_ogversion.R
# ==============================================================================

# 1. Variables ----
test_job <- "testCOO_definitive"
dir_data <- "COOannot"
name_coo_annotation = "_NoCOOAEnhATAC"

# 2. Hardcoded variables -----

# Annotations combinations tested:
# Annotations have to be in the same order as in the .annot files !!
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
                    "comb2DD"= c("AEnhATAC", "APromoterATAC", "PIRs","PoEnhATAC","PrEnhATAC", "CDS", "COO")
)

coo.list <- list( "COO1" = c("PreProB", "ProB", "PreB"),
                  "COO2" = c("nB", "memB"),
                  "COO3" = c("GCB"),
                  "COO4" = c("PC"),
                  "COO5" = c("nB", "GCB", "memB"),
                  "COO6" = c("GCB","PC") )

traits.list <- list("G1" = c("pedBALL_ICR", "pedBALLer_ICR", "pedBALLhyper_ICR","pedBALL_Alyssa","BALL_Alyssa"),
                    "G2" = c("CLL_ILH", "MZL_ILH"),
                    "G3" = c("FL_ILH", "HL_ILH", "DLBCL_ILH", "LYMPHOMA_UKBB"),
                    "G4" = c("MM_ICR"),
                    "G5" = c("FL_ILH", "HL_ILH", "DLBCL_ILH", "LYMPHOMA_UKBB"),
                    "G6" = c("FL_ILH", "HL_ILH", "DLBCL_ILH", "LYMPHOMA_UKBB", "MM_ICR")
                    )

# 3. Hardcoded paths -----
# LD-scores and weigths files
ld <- paste0("/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/polyFUN/data/annotations/", dir_data, "/", dir_data, ".")
w <- paste0("/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/polyFUN/data/annotations/", dir_data, "/", dir_data, "_w.")
# Job files
jobs_dir_mn5 <- "/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/polyFUN/jobs/05_runPolyFun/"
logs_dir_mn5 <- "/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/polyFUN/logs/05_runPolyFun/"
jobs_dir_local <- "~/scratch/B_cell_roadmap_2/polyFUN/jobs/05_runPolyFun"
logs_dir_local <- "~/scratch/B_cell_roadmap_2/polyFUN/logs/05_runPolyFun"
ifelse(!dir.exists(jobs_dir_local), dir.create(jobs_dir_local, recursive = T), FALSE)
ifelse(!dir.exists(logs_dir_local), dir.create(logs_dir_local, recursive = T), FALSE)

# 4. Read first line of annotation file to get the annotation order ----
file.header <- paste0("~/scratch/B_cell_roadmap_2/polyFUN/data/annotations/", dir_data, "/", dir_data, ".22.annot")
first_line <- readLines(file.header, n=1)
first_line <- strsplit(first_line, "\t")

# 5. Prepare commands ----
# Select combination annotations to run
combis.selected = c("comb2DD")

# Iterate for each trait group (same COO andn therefore, same annotations used)
for(trait.group in names(traits.list)){
  print(trait.group)
  gwas_selected <- traits.list[[trait.group]]
  
  coo.name <- gsub("G", "COO", trait.group)
  coo.celltypes <- coo.list[[coo.name]]
  
  total_commands = 0
  
  # Create one job for each trait group
  output_commands_file <- paste0(jobs_dir_local, "/05_runPolyfun_", trait.group, '_', test_job,"_commands.sh")
  cat("",file=output_commands_file)
  
  for(trait in gwas_selected){
    print(trait)
    
    sumstats <- paste0("/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/polyFUN/data/sumstats/", trait, "_munged.parquet")
    
    for(combi in combis.selected){
      print(combi)
      
      # Take annotation list to use for each trait
      annotations_list = annots_list[[combi]]
      
      annot.vector <- with(expand.grid(coo.celltypes, setdiff(annotations_list, c("CDS", "COO"))),
                           paste0(Var1, "_", Var2, "$"))
      pattern.search <- "Base|%s"
      if("CDS" %in% annotations_list){ pattern.search <- paste0("CDS|", pattern.search)}
      if("COO" %in% annotations_list){ pattern.search <- paste0(coo.name, name_coo_annotation ,"|", pattern.search)}

      annot.toprint <- grep(sprintf(pattern.search, paste0(annot.vector, collapse="|")), first_line[[1]], value = T)

      # Create ouput directory
      outdir_local <- paste0("~/scratch/B_cell_roadmap_2/polyFUN/results/",trait, "/", trait, '_', coo.name, "_", combi)
      ifelse(!dir.exists(outdir_local), dir.create(outdir_local, recursive = T), FALSE)
      
      # PolyFun command
      output_prefix <- paste0("/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/polyFUN/results/", trait, 
                              "/", trait, '_', coo.name, "_", combi,
                              "/", trait, '_', coo.name, "_", combi)
      
      command_polyfun <- paste0("python /gpfs/scratch/bsc08/bsc008246/Programs/polyfun/polyfun.py ",
                                "--compute-h2-L2 ",
                                "--no-partitions ",
                                "--anno ", paste(annot.toprint, collapse=','), " ",
                                "--output-prefix ",output_prefix, " ",
                                "--sumstats ", sumstats," ",
                                "--ref-ld-chr ", ld, " ",
                                "--w-ld-chr ", w, "\n")
      
      cat(command_polyfun, file=output_commands_file, append = T)
      total_commands = total_commands + 1
    
    }
  }
  
  # Print slurm command ----
  slurm_command <- paste0(jobs_dir_local, "/05_runPolyfun_", trait.group, "_", test_job,".cmd")
  fileConn <- file(slurm_command, open="wt")
  
  job <- paste0("#!/bin/bash\n",
                "#SBATCH --job-name=05_runPolyfun_", trait.group, "_", test_job, "\n",
                "#SBATCH -D ", logs_dir_mn5,"\n",
                "#SBATCH --output=",logs_dir_mn5, "05_runPolyfun_", trait.group, "_", test_job,"_%j.out\n",
                "#SBATCH --error=",logs_dir_mn5, "05_runPolyfun_", trait.group, "_", test_job,"_%j.err\n",
                "#SBATCH --cpus-per-task=",20,"\n",
                "#SBATCH --ntasks=",15,"\n",
                "#SBATCH --time=24:00:00\n",
                "#SBATCH --qos=gp_bscls\n",
                "#SBATCH --account=bsc08\n\n",
                
                "export GREASY_LOGFILE=",logs_dir_mn5, "05_runPolyfun_", trait.group, "_", test_job, ".log\n\n",
                
                
                "/apps/GPP/GREASY/2.2.4.2/INTEL/IMPI/bin/greasy ",jobs_dir_mn5,  "05_runPolyfun_", trait.group, '_', test_job, "_commands.sh\n"
  )

  cat(job, file = fileConn)
  close(fileConn)
  # system2(command="sbatch", arg = slurm_command )
  # copy and paste these commands in a command file and run as sbatch job
  cat(paste0("sbatch 05_runPolyfun_", trait.group, "_", test_job,".cmd\n"))
}


# To prepare sbatch commands
# copy and paste these commands in a command file and run as sbatch job
# remember to conda activate polyfun
for(trait.group in names(traits.list)){
  cat(paste0("sbatch 05_runPolyfun_", trait.group, "_", test_job,".cmd\n"))
}

