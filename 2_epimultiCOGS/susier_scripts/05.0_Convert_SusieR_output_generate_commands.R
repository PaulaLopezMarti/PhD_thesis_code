gwas_selected <- c("pedBALL_ICR", "pedBALLer_ICR", "pedBALLhyper_ICR","pedBALL_Alyssa","BALL_Alyssa", "LYMPHOMA_UKBB","FL_ILH","HL_ILH","DLBCL_ILH","MZL_ILH","CLL_ILH","MM_ICR")
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
tests_list <- names(annots_list)
tests_list <- grep("comb3", tests_list, value = T, invert = T)
# tests_list <- grep("comb2DD", tests_list, value = T, invert = F)

# Write commands
file_basename <- "05B_Convert_SusieR_ouput_for_mCOGS_epitests"
jobs_dir_local <- "~/scratch/B_cell_roadmap_2/susie/jobs/05_convert_susier_for_mCOGS"
output_commands_file <- paste0(jobs_dir_local, "/", file_basename,"_commands.sh")
cat("",file=output_commands_file)

for(trait in gwas_selected){
  for(annot in tests_list){
    command <- paste0("/gpfs/apps/MN5/GPP/R/4.3.2/INTEL/bin/Rscript /gpfs/scratch/bsc08/bsc008246/PhD_thesis_code/3_Epigenetic_finemapping/susier_scripts/05B_Convert_SusieR_ouput_for_mCOGS.R ", trait, " ", annot, "\n")
    cat(command, file=output_commands_file, append = T)
    }
}

