# ==============================================================================
# Generate report about summary statistics and reference panel matching
#
# To be run in local
# 
# Original script:
# /home/plopez/scratch/B_cell_roadmap_2/polyFUN/scripts/checks/02.2_Check_sumstats_flip_strands_summarise.R
# ==============================================================================

library(tidyverse)

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
  
  # snps to flip
  snps.toflip = flip[flip$`_FLIP_` == T,]$SNP
  
  # Conserved alleles in 1KG: these are all SNPs that are not NA in noflip + SNPs that should be flipped 
  conserved.snps = noflip[!is.na(noflip$`_REV_`) | noflip$SNP %in% snps.toflip,]
  snps.filter.ambigous = flip[!is.na(flip$`_REV_`) & flip$`_FLIP_` == F,] # remove strand flip alleles
  snps.nofilter.ambigous = noflip[!is.na(noflip$`_REV_`)] # here strand flip alleles are NA because they could not be matched to 1KG
  ambigous.alleles = snps.nofilter.ambigous[!(snps.nofilter.ambigous$SNP %in% snps.filter.ambigous$SNP),]
  
  # Save summary data
  summary_cs[i, "trait"] <- gwas_data[gwas_data$trait_id == trait,]$trait_name_plot_short
  
  # Total number of SNPs
  summary_cs[i, "SNPs"] <- nrow(flip)
  summary_cs[i, "sign"] <- nrow(flip[flip$P < p_sign,])
  summary_cs[i, "prom"] <- nrow(flip[flip$P < p_prom,])
  
  # Number of conserved alleles in 1KG: 
  summary_cs[i, "SNPs.1kg"] <- nrow(conserved.snps)
  summary_cs[i, "sign.1kg"] <- nrow(conserved.snps[conserved.snps$P < p_sign,])
  summary_cs[i, "prom.1kg"] <- nrow(conserved.snps[conserved.snps$P < p_prom,])
  
  # Number of variants that require strand flip
  summary_cs[i, "SNPs.flip"] <- n_distinct(snps.toflip)
  summary_cs[i, "sign.flip"] <- nrow(flip[flip$`_FLIP_` == T & flip$P < p_sign,])
  summary_cs[i, "prom.flip"] <- nrow(flip[flip$`_FLIP_` == T & flip$P < p_prom,])
  
  # Number of ambigous alleles: 
  # NO: number of SNPs loss in flip that are not loss in noflip substracting the aleles that should be flipped
  # n_ambigous = nrow(flip[is.na(flip$`_REV_`)]) - nrow(noflip[is.na(noflip$`_REV_`)]) - nrow(flip[flip$`_FLIP_` == T,])
  # n_ambigous.sign = nrow(flip[is.na(flip$`_REV_`) & flip$P < p_sign]) - nrow(noflip[is.na(noflip$`_REV_`) & noflip$P < p_sign]) - nrow(flip[flip$`_FLIP_` == T & flip$P < p_sign,])
  # n_ambigous.prom = nrow(flip[is.na(flip$`_REV_`) & flip$P < p_prom]) - nrow(noflip[is.na(noflip$`_REV_`) & noflip$P < p_prom]) - nrow(flip[flip$`_FLIP_` == T  & flip$P < p_prom,])
  n_ambigous = nrow(ambigous.alleles)
  n_ambigous.sign = nrow(ambigous.alleles[ ambigous.alleles$P < p_sign,])
  n_ambigous.prom = nrow(ambigous.alleles[ ambigous.alleles$P < p_prom,])
  
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
  remove.flipstrands = noflip[noflip$`_FLIP_` == F,] # no flip have already removed the flip strands because they don't match with 1KG
  summary_cs[i, "SNPs.withambigous"] <- nrow(remove.flipstrands)
  summary_cs[i, "sign.withambigous"] <- nrow(remove.flipstrands[ remove.flipstrands$P < p_sign,]) 
  summary_cs[i, "prom.withambigous"] <- nrow(remove.flipstrands[ remove.flipstrands$P < p_prom,])
  
  # BY MAF
  summary_cs[i, "SNPs.maf05"] <-  nrow(remove.flipstrands[remove.flipstrands$MAF < 0.05,])
  summary_cs[i, "sign.maf05"] <-  nrow(remove.flipstrands[remove.flipstrands$MAF < 0.05 & remove.flipstrands$P < p_sign,])
  summary_cs[i, "prom.maf05"] <- nrow(remove.flipstrands[remove.flipstrands$MAF < 0.05 & remove.flipstrands$P < p_prom,])
  
  summary_cs[i, "SNPs.maf01"] <-  nrow(remove.flipstrands[remove.flipstrands$MAF < 0.01,])
  summary_cs[i, "sign.maf01"] <-  nrow(remove.flipstrands[remove.flipstrands$MAF < 0.01 & remove.flipstrands$P < p_sign,])
  summary_cs[i, "prom.maf01"] <- nrow(remove.flipstrands[remove.flipstrands$MAF < 0.01 & remove.flipstrands$P < p_prom,])
  
  summary_cs[i, "SNPs.maf0"] <-  nrow(remove.flipstrands[remove.flipstrands$MAF == 0,])
  summary_cs[i, "sign.maf0"] <-  nrow(remove.flipstrands[remove.flipstrands$MAF == 0  & remove.flipstrands$P < p_sign,])
  summary_cs[i, "prom.maf0"] <- nrow(remove.flipstrands[remove.flipstrands$MAF == 0 & remove.flipstrands$P < p_prom,])
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

summary_wide$withoutambigous = NULL
names(summary_wide) <- c("Trait", "SNPs.selected", "total_SNPs", 
                         "In_1KG", "Strand_flip", "Ambiguous", "Allele_swap",
                         "Filter_strand_flips",
                         "MAF_0.05", "MAF_0.01", "MAF_0")

# change SNPs selected labels to be more meaningful
summary_wide$SNPs.selected = gsub("SNPs", "all", summary_wide$SNPs.selected)
summary_wide$SNPs.selected <- gsub("prom", "<1e-6", summary_wide$SNPs.selected)
summary_wide$SNPs.selected <- gsub("sign", "<5e-8", summary_wide$SNPs.selected)
summary_wide$SNPs.selected <- factor(summary_wide$SNPs.selected, levels = c("all","<1e-6", "<5e-8" ))
summary_wide$Trait <- factor(summary_wide$Trait, levels = traits.short)
summary_wide <- summary_wide[order(summary_wide$SNPs.selected, summary_wide$Trait),]

data.table::fwrite(summary_wide, file="~/scratch/B_cell_roadmap_2/polyFUN/checks/Report_alltraits_summary_reference_matching.tsv", sep="\t")
