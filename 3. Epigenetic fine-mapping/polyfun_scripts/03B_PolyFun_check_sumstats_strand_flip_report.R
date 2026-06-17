# ==============================================================================
# Generate report about summary statistics and reference panel matching
#
# To be run in local
# 
# Original script:
# /home/plopez/scratch/B_cell_roadmap_2/polyFUN/scripts/checks/02.2_Check_sumstats_flip_strands_summarise.R
# ==============================================================================


# ==============================================================================
# 1. Global configuration
# ==============================================================================
gwas_info_file = "~/scratch/Data/GWAS/SummaryStatistics_AllGWAS_complete.tsv"
inputdir <- paste0("~/scratch/B_cell_roadmap_2/polyFUN/checks/sumstats/")
bcellmalignancies <- c("pedBALL_ICR", "pedBALLer_ICR", "pedBALLhyper_ICR","pedBALL_Alyssa","BALL_Alyssa", "LYMPHOMA_UKBB","FL_ILH","HL_ILH","DLBCL_ILH","MZL_ILH","CLL_ILH","MM_ICR")
p_sign <- 5e-8
p_prom <- 10^-6

# Select GWAS of interest ----
gwas_data.og <- data.table::fread(gwas_info_file)
gwas_data <- gwas_data.og[,c("trait_id","trait_name_plot_short")]
gwas_data <- unique(gwas_data)
gwas_data = gwas_data[gwas_data$trait_id %in% bcellmalignancies,]

trait <- "pedBALL_ICR"

# Build summary table ----
colnames_summary <- c("trait" )
summary_cs <- data.frame(matrix(ncol = length(colnames_summary), nrow = 0))
names(summary_cs) <- colnames_summary

for(i in 1:length(bcellmalignancies)){
  print(trait)
  trait <- bcellmalignancies[i]
  
  # Read data and filter out non-necessary columns
  flip.path <- paste0(inputdir, trait, '_cleansumstats_withflip.tsv')
  noflip.path <- paste0(inputdir, trait, '_cleansumstats_NOflip.tsv')

  flip <- data.table::fread(flip.path)  
  noflip <- data.table::fread(noflip.path) 
  
  flip <- flip[,c("SNP","SNP_og","CHR","BP","P","Z","MAF","_FLIP_","_REV_")]
  noflip <- noflip[,c("SNP","SNP_og","CHR","BP","P","Z","MAF","_FLIP_","_REV_")]
  
  # Save summary data
  summary_cs[i, "trait"] <- gwas_data[gwas_data$trait_id == trait,]$trait_name_plot_short
  
  # Total number of SNPs
  summary_cs[i, "SNPs"] <- nrow(flip)
  summary_cs[i, "sign"] <- nrow(flip[flip$P < p_sign,])
  summary_cs[i, "prom"] <- nrow(flip[flip$P < p_prom,])
  
  # Number of conserved alleles in 1KG
  summary_cs[i, "SNPs.1kg"] <- nrow(noflip[!is.na(noflip$`_REV_`)])
  summary_cs[i, "sign.1kg"] <- nrow(noflip[!is.na(noflip$`_REV_`) & noflip$P < p_sign,])
  summary_cs[i, "prom.1kg"] <- nrow(noflip[!is.na(noflip$`_REV_`) & noflip$P < p_prom,])
  
  # Number of variants that require strand flip
  summary_cs[i, "SNPs.flip"] <-  nrow(flip[flip$`_FLIP_` == T,])
  summary_cs[i, "sign.flip"] <-  nrow(flip[flip$`_FLIP_` == T & flip$P < p_sign,])
  summary_cs[i, "prom.flip"] <- nrow(flip[flip$`_FLIP_` == T & flip$P < p_prom,])
  
  # Number of ambigous alleles: 
  # number of SNPs loss in flip that are not loss in noflip
  # substracting the aleles that should be flipped
  n_ambigous = nrow(flip[is.na(flip$`_REV_`)]) - nrow(noflip[is.na(noflip$`_REV_`)]) - nrow(flip[flip$`_FLIP_` == T,])
  n_ambigous.sign = nrow(flip[is.na(flip$`_REV_`) & flip$P < p_sign]) - nrow(noflip[is.na(noflip$`_REV_`) & noflip$P < p_sign]) - nrow(flip[flip$`_FLIP_` == T & flip$P < p_sign,])
  n_ambigous.prom = nrow(flip[is.na(flip$`_REV_`) & flip$P < p_prom]) - nrow(noflip[is.na(noflip$`_REV_`) & noflip$P < p_prom]) - nrow(flip[flip$`_FLIP_` == T  & flip$P < p_prom,])
  
  summary_cs[i, "SNPs.ambigous"] <- n_ambigous
  summary_cs[i, "sign.ambigous"] <- n_ambigous.sign
  summary_cs[i, "prom.ambigous"] <- n_ambigous.prom
  
  # Number of variants that require allele swap
  summary_cs[i, "SNPs.rev"] <-  nrow(noflip[noflip$`_REV_` == T,])
  summary_cs[i, "sign.rev"] <-  nrow(noflip[noflip$`_REV_` == T & noflip$P < p_sign,])
  summary_cs[i, "prom.rev"] <- nrow(noflip[noflip$`_REV_` == T & noflip$P < p_prom,])
  
  # Number of variants removing strand flip and ambigous alleles
  summary_cs[i, "SNPs.withoutambigous"] <- nrow(flip[!is.na(flip$`_REV_`)])
  summary_cs[i, "sign.withoutambigous"] <- nrow(flip[!is.na(flip$`_REV_`) & flip$P < p_sign,])
  summary_cs[i, "prom.withoutambigous"] <- nrow(flip[!is.na(flip$`_REV_`) & flip$P < p_prom,])
  
  # Number of variants removing strand flip
  summary_cs[i, "SNPs.withambigous"] <- nrow(noflip[noflip$`_FLIP_` == F,]) - nrow(flip[flip$`_FLIP_` == T,])
  summary_cs[i, "sign.withambigous"] <- nrow(noflip[noflip$`_FLIP_` == F & noflip$P < p_sign,]) - nrow(flip[flip$`_FLIP_` == T & flip$P < p_sign,])
  summary_cs[i, "prom.withambigous"] <- nrow(noflip[noflip$`_FLIP_` == F & noflip$P < p_prom,]) - nrow(flip[flip$`_FLIP_` == T & flip$P < p_prom,])
  
  # BY MAF
  # summary_cs[i, "SNPs.maf05"] <-  nrow(noflip[noflip$MAF < 0.05,])
  # summary_cs[i, "sign.maf05"] <-  nrow(noflip[noflip$MAF < 0.05 & noflip$P < p_sign,])
  # summary_cs[i, "prom.maf05"] <- nrow(noflip[noflip$MAF < 0.05 & noflip$P < p_prom,])
  # 
  # summary_cs[i, "SNPs.maf01"] <-  nrow(noflip[noflip$MAF < 0.01,])
  # summary_cs[i, "sign.maf01"] <-  nrow(noflip[noflip$MAF < 0.01 & noflip$P < p_sign,])
  # summary_cs[i, "prom.maf01"] <- nrow(noflip[noflip$MAF < 0.01 & noflip$P < p_prom,])
  # 
  # summary_cs[i, "SNPs.maf0"] <-  nrow(noflip[noflip$MAF == 0,])
  # summary_cs[i, "sign.maf0"] <-  nrow(noflip[noflip$MAF == 0  & noflip$P < p_sign,])
  # summary_cs[i, "prom.maf0"] <- nrow(noflip[noflip$MAF == 0 & noflip$P < p_prom,])
}


data.table::fwrite(summary_cs, file="~/scratch/B_cell_roadmap_2/polyFUN/checks/alltraits_summary_flip_ambigous_v2.tsv", sep="\t")

# Format report table
summary_long <- summary_cs %>%
  pivot_longer(
    cols = -trait,
    names_to = "var",
    values_to = "value"
  ) %>%
  mutate(
    var = str_replace(var, "^SNPs$", "SNPs.all"),
    metric = str_extract(var, "^[^.]+"),
    category = str_extract(var, "(?<=\\.).*"),
    category = ifelse(is.na(category), "all", category)
  ) %>%
  select(trait, metric, category, value)

summary_wide <- summary_long %>%
  pivot_wider(
    names_from = category,
    values_from = value
  )

names(summary_wide) <- c("Trait", "SNPs selected", "Number SNPs", "Included in 1KG", "Flip strands", "Ambigous", "Swap alleles",
                         "Total SNPs conserved (removing ambigous)",
                         "Total SNPs conserved (with ambigous, remove strand flip)")
summary_wide$`SNPs selected` = gsub("SNPs", "all", summary_wide$`SNPs selected`)

data.table::fwrite(summary_wide, file="~/scratch/B_cell_roadmap_2/polyFUN/checks/Report_alltraits_summary_reference_matching.tsv", sep="\t")
