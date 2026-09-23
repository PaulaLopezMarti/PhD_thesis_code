# ============================================================
# Merge fine-mapping QC all traits all epigenetic tests
# Original script:
# /home/plopez/scratch/B_cell_roadmap_2/susie/scripts/checks/02.2_Check_bad_finemapping_MN5_mergetables.R
# ============================================================

library(tidyverse)
library(data.table)

# ============================================================
# 1. Global configuration & Variables
# ============================================================
input.dir = "~/scratch/B_cell_roadmap_2/susie/analysis/checks/select_LDblocks/"
output.dir.alltraits = "~/scratch/B_cell_roadmap_2/susie/analysis/checks/LDblocks_alltraits/"

bcellmalignancies <- c("pedBALL_ICR", "pedBALLer_ICR", "pedBALLhyper_ICR","pedBALL_Alyssa","BALL_Alyssa", "LYMPHOMA_UKBB",
                       "FL_ILH","HL_ILH","DLBCL_ILH","MZL_ILH","CLL_ILH","MM_ICR")
coo.list <- list("MM_ICR" = "COO4", "pedBALL_ICR" = "COO1", "pedBALLer_ICR" = "COO1", "pedBALLhyper_ICR" = "COO1","BALL_Alyssa" = "COO1", "pedBALL_Alyssa" = "COO1", "CLL_ILH" = "COO2", "MZL_ILH" = "COO2", "FL_ILH" = "COO3", "HL_ILH" = "COO3", "DLBCL_ILH" = "COO3", "LYMPHOMA_UKBB" = "COO5")
annots_list <- list(
  "comb1A" = 1,  "comb1B" = 3,  "comb1C" = 5,  "comb1D" = 7,
  "comb1AA" = 2, "comb1BB" = 4, "comb1CC" = 6, "comb1DD" = 8,
  "comb2A" = 9,  "comb2B" = 11, "comb2C" = 13, "comb2D" = 15,
  "comb2AA" = 10, "comb2BB" = 12, "comb2CC" = 14, "comb2DD" = 16,
  "standard" = "S"
)
annots_vec <- unlist(annots_list)

# Import GWAS metadata
gwas_data.og <- data.table::fread("/home/plopez/scratch/Data/GWAS/SummaryStatistics_AllGWAS_complete.tsv")
gwas_data <- gwas_data.og[,c("trait_id","trait_name_plot_short", "N_cases")]
gwas_data <- unique(gwas_data)
gwas_data = gwas_data[gwas_data$trait_id %in% bcellmalignancies,]
traits.short <- c( "MM", "CLL", "LYM",  "MZL","FL",  "DLBCL","HL", "pB-ALL", "pB-ALL HD", "pB-ALL ER", "B-ALL 2", "pB-ALL 2")


# ============================================================
# 2. Summarise LD block QC results for each trait
# Aggregates LD block QC across traits
# and produces a global summary table stratified by p-value bins
# ============================================================

for(i in 1:length(bcellmalignancies)){
  
  trait <- bcellmalignancies[i]
  print(trait)
  coo.selected = coo.list[[trait]]
  trait.short = gwas_data[gwas_data$trait_id == trait,]$trait_name_plot_short
  
  # ============================================================
  # 2.1 Loop over epigenetic annotation combinations
  # ============================================================
  
  for(x in 1:length(names(annots_vec))){
    
    comb.selected <- names(annots_vec[x])
    comb.id <- annots_list[[comb.selected]]
    
    # Rename for standard fine-mapping
    if(comb.selected == "standard"){ 
      df <- data.table::fread(paste0(input.dir, trait, "_CV5_cov90_LDblocks_selection.tsv"))
    } else {
      df <- data.table::fread(paste0(input.dir, trait, "_CV5_", coo.selected, "_",comb.selected , "_cov90_LDblocks_selection.tsv"))
    }
    

    # ============================================================
    # 2.2 Define QC categories
    # ============================================================
    
    ## Each LD block is assigned to:
    ##  - p-value category
    ##  - QC category based on lambda, convergence, and CS signals
    
    df <- df %>%
      mutate(
        pval_group = case_when(
          min_pval < 5e-8 ~ "<5e-8",
          min_pval < 1e-6 ~ "1e-6 - 5e-8",
          TRUE ~ ">1e-6"
        ),
        
        QC = case_when(
          lambda_over_0.4 == "BAD" ~ "Lambda > 0.4",
          lambda_over_0.4 == "ok" & convIBSS_CV5 == "BAD" ~ "No converge",
          # lambda_over_0.4 == "ok" & convIBSS_CV5 == "ok" & nCS_CV5 >0 ~ "Has CS",
          lambda_over_0.4 == "ok" & convIBSS_CV5 == "ok" & pass5e5_CV5 == "yes" ~ "Good CV5 finemapping",
          lambda_over_0.4 == "ok" & convIBSS_CV5 == "ok" & pass5e5_CV5 == "no" ~ "Bad CV5 finemapping",
          lambda_over_0.4 == "ok" & convIBSS_CV5 == "ok" & nCS_CV5 == 0 ~ "Not CS"
        )
      )
    
    # ============================================================
    # 2.3 Summarise QC per p-value group
    # ============================================================
    
    df3 <- df %>% group_by(pval_group, QC) %>% summarise(nBlocks = n_distinct(LDblock),
                                                         .groups = "drop_last")
    
    ## Total LD blocks and blocks with at least one CS
    df3.total = df %>% group_by(pval_group) %>% summarise(total_LDblocks = n_distinct(LDblock),
                                                          CV5_has_CS = n_distinct(LDblock[nCS_CV5 > 0 & !(QC %in% c("Lambda > 0.4", "No converge"))])
                                                          )
    
    # Reshape QC table to wide format
    df4 = pivot_wider(df3, names_from = QC, values_from = "nBlocks", values_fill=0)
    df4 = merge(df3.total, df4, by="pval_group")
    df4$trait = trait.short
    
    
    # ============================================================
    # 2.4 Ensure consistent column structure
    # ============================================================
    
    cols.tosave = c("trait","total_LDblocks", "pval_group","Lambda > 0.4", "No converge", "CV5_has_CS","Good CV5 finemapping", "Bad CV5 finemapping", "Not CS" )
    cols.selected = intersect(names(df4), cols.tosave)
    df4 = df4[,cols.selected]
    
    ## Add missing columns if absent in a given subset
    if(length(setdiff(cols.tosave, names(df4) ) > 0)){
      for(col.missing in setdiff(cols.tosave, names(df4) )){
        df4[[col.missing]] <- 0
      }
    }
    
    # Add for MZL GWAS the row of sign. snps (there is no sign snp outside MHC)
    if(!("<5e-8" %in% unique(df4$pval_group))){
      empty_row <- df4 %>% slice(1) %>%    mutate(
          pval_group = "<5e-8",
          trait = trait,
          across(-c(pval_group, trait), ~ 0)
        )
      
      df4 <- bind_rows(df4, empty_row)
    }
    
    # Order p-value groups
    # df4$pval_group <- factor(df4$pval_group , levels = c("<5e-8", "<1e-6","<5e-5"))
    # df4 <- df4[order(df4$pval_group), ]
    
    # ============================================================
    # 2.5 Add total row (sum across p-value groups)
    # ============================================================
    
    df.total <- df4 %>% bind_rows( df4 %>%
          summarise(
            across(where(is.numeric), sum),
            pval_group = "Total",
            trait = first(trait)
          )
      )

    df.total$epitest = comb.id
    df.total <- df.total[, c("trait", "epitest", setdiff(cols.tosave, c("trait", "epitest")))]
    df.total$epitest <- as.character(df.total$epitest)
    
    
    # ============================================================
    # 2.6 Merge epitests
    # ============================================================
    
    if(x == 1){ df.trait <- df.total} else{ df.trait <- bind_rows(df.trait, df.total )}
    
  }
  
  # Merge traits
  if(i == 1){ df.all <- df.trait} else{ df.all <- bind_rows(df.all, df.trait )}
  
}

# Final clean-up
df.all[is.na(df.all$`No converge` ), "No converge"] = 0
df.all[is.na(df.all$`Good CV5 finemapping` ), "Good CV5 finemapping"] = 0
df.all[is.na(df.all$`Bad CV5 finemapping` ), "Bad CV5 finemapping"] = 0
df.all = df.all[,c("trait","epitest","pval_group","total_LDblocks", "Lambda > 0.4", "No converge" , "Not CS", "CV5_has_CS","Good CV5 finemapping", "Bad CV5 finemapping"  )]
df.all$pval_group <- gsub(" ", "", df.all$pval_group)

# ============================================================
# 3. Save outputs
# ============================================================

data.table::fwrite(df.all, paste0(output.dir.alltraits, "LDblocks_finemapping_QC_summary_alltraits_allepitests.tsv"), sep="\t")


# ============================================================
# 3.1 Format publication-ready table
# ============================================================

# keep only standard fine-mapping the epigenetic combination used in the paper
df.filtered <- df.all[df.all$epitest %in% c("S", "16"), ]
data.table::fwrite(df.filtered, paste0(output.dir.alltraits, "LDblocks_finemapping_QC_summary_for_publication.tsv"), sep="\t")


df.filtered <- df.all[df.all$pval_group != "Total", ]
df.filtered$trait <- factor(df.filtered$trait, levels=traits.short)
df.filtered$pval_group <- factor(df.filtered$pval_group , levels = c("<5e-8", "1e-6-5e-8",">1e-6"))
df.filtered <- df.filtered[order(df.filtered$trait, df.filtered$pval_group),]

# Standard
df.filtered2 <- df.filtered[df.filtered$epitest %in% c("S"), ]
df.filtered2$epitest <- NULL
data.table::fwrite(df.filtered2, paste0(output.dir.alltraits, "LDblocks_finemapping_QC_summary_for_publication_Standard.tsv"), sep="\t")

# Epigenetic
df.filtered3 <- df.filtered[df.filtered$epitest %in% c("16"), ]
df.filtered3$epitest <- NULL
data.table::fwrite(df.filtered3, paste0(output.dir.alltraits, "LDblocks_finemapping_QC_summary_for_publication_Epigenetic.tsv"), sep="\t")



# ============================================================
# 3.2 Format table for thesis
# ============================================================

# keep only total SNPs for all epitests
df.filtered <- df.all[df.all$pval_group == "Total", ]
df.filtered$epitest <- factor(df.filtered$epitest, levels = c(1:16, "S"))
df.filtered$trait <- factor(df.filtered$trait, levels=traits.short)
df.filtered <- df.filtered[order(df.filtered$trait, df.filtered$epitest),]

data.table::fwrite(df.filtered, paste0(output.dir.alltraits, "LDblocks_finemapping_QC_summary_for_thesis.tsv"), sep="\t")

