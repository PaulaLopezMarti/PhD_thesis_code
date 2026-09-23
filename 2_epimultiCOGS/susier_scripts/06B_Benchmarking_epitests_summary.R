library(tidyverse)
# library(cowplot)
library(ggplot2)

###############################################################################
# 1. Global settings
###############################################################################
bcellmalignancies <- c( "MM_ICR", "CLL_ILH", "pedBALL_ICR", 
                        "DLBCL_ILH","FL_ILH", "HL_ILH", "MZL_ILH", "LYMPHOMA_UKBB",
                        "pedBALLer_ICR", "pedBALLhyper_ICR",
                        "BALL_Alyssa","pedBALL_Alyssa")

coo.list <- list("MM_ICR" = "COO4", "pedBALL_ICR" = "COO1", "pedBALLer_ICR" = "COO1", "pedBALLhyper_ICR" = "COO1","BALL_Alyssa" = "COO1", "pedBALL_Alyssa" = "COO1", "CLL_ILH" = "COO2", "MZL_ILH" = "COO2", "FL_ILH" = "COO3", "HL_ILH" = "COO3", "DLBCL_ILH" = "COO3", "LYMPHOMA_UKBB" = "COO5")

comb.order <- c("comb1A",  "comb1B" , "comb1C" , "comb1D" , "comb1AA" ,"comb1BB" ,"comb1CC", 
                "comb1DD", "comb2A" , "comb2B" , "comb2C" , "comb2D",  "comb2AA" ,"comb2BB" ,"comb2CC" ,"comb2DD","standard")
trait = "DLBCL_ILH"

tables.summary.path <- "~/scratch/B_cell_roadmap_2/susie/analysis/tables_summary_finemapping/"
outputdir.path <- "~/scratch/B_cell_roadmap_2/susie/analysis/benchmarking_epitests/"

annots_list <- list("comb1A" = 1,
                    "comb1B" = 3,
                    "comb1C" = 5,
                    "comb1D" = 7,
                    "comb1AA" = 2,
                    "comb1BB" = 4,
                    "comb1CC" = 6,
                    "comb1DD" = 8,
                    
                    "comb2A"= 9,
                    "comb2B"= 11,
                    "comb2C"= 13,
                    "comb2D"= 15,
                    "comb2AA"= 10,
                    "comb2BB"= 12,
                    "comb2CC"= 14,
                    "comb2DD"= 16,
                    "standard" = "S"
)

annots_vec <- unlist(annots_list)



# gwas metadata
gwas_data.og <- data.table::fread("/home/plopez/scratch/Data/GWAS/SummaryStatistics_AllGWAS_complete.tsv")
gwas_data <- gwas_data.og[,c("trait_id","trait_name_plot_short")]
gwas_data <- unique(gwas_data)
gwas_data = gwas_data[gwas_data$trait_id %in% bcellmalignancies,]
traits.short <- c( "MM", "CLL", "LYM",  "MZL","FL",  "DLBCL","HL", "pB-ALL", "pB-ALL HD", "pB-ALL ER", "B-ALL 2", "pB-ALL 2")




##################################################
## 2. Create summary table all traits
##################################################
for(x in 1:length(bcellmalignancies)){
  trait <- bcellmalignancies[x]
  print(trait)
  coo.selected <- coo.list[trait]
  input_dir <- paste0(tables.summary.path, trait)
  output_dir <- paste0(outputdir.path, trait)
  dir.create(output_dir, showWarnings = F, recursive = T)
  
  # Find list of files to process
  files.sum <- list.files(input_dir)
  files.sum <- grep("comb3", files.sum, value = T, invert = T)
  files.sum <- grep("CV5", files.sum, value = T, invert = F)
  files.sum <- grep(paste0("_", coo.selected, "_|_standard"), files.sum, value=T, invert = F) # Select the correct COO
  
  # Read all epitests summaries
  print(length(files.sum))
  for(i in 1:length(files.sum)){
    # print(i)
    file.selected = files.sum[i]
    df1 <- data.table::fread(paste0(input_dir,"/", file.selected), sep="\t")
    if(i == 1){
      df.all <- df1
      } else{df.all <- rbind(df.all, df1)}
  }
  
  # Merge all traits
  if(x == 1){
    df.alltraits <- df.all
  }else{
    df.alltraits <- rbind(df.alltraits, df.all)}
  
}
 
data.table::fwrite(df.alltraits, file=paste0(outputdir.path, "Table_finemapping_benchmarking_epitests.tsv"), sep="\t")

 
##################################################
## 3. Create summary plots
##################################################
df.alltraits <- data.table::fread(paste0(outputdir.path, "Table_finemapping_benchmarking_epitests.tsv"), sep="\t")
df.alltraits <- merge(df.alltraits, gwas_data)
df.alltraits$annot_id <-  annots_vec[match(df.alltraits$epitest, names(annots_vec))]

##################################################
## 3.1 Number of Credible sets
##################################################
df.f <- df.alltraits[df.alltraits$selectedSNPS == "allSNPs",]
df.f$annot_id <- factor(df.f$annot_id, levels = rev(c(1:16, "S")))
df.f$trait_name_plot_short <- factor(df.f$trait_name_plot_short, levels = traits.short)

df.f$method_group <- ifelse(df.f$epitest == "standard", "standard",
                            ifelse(df.f$epitest == "comb2DD", "selected", "other"))

# one hline value per trait
hline.df <- df.f %>%
  filter(epitest == "standard") %>%
  select(trait_name_plot_short, hline = n_cs)


p1 <- ggplot(df.f, aes(x = annot_id, y = n_cs)) +
  
  geom_point(aes(color=method_group), size = 2.5) +
  
  scale_color_manual(values = c(
    "standard" = "#9b0303",
    "selected"  = "#1f77b4",
    "other"    = "black"
  )) +
  
  geom_hline(
    data = hline.df,
    aes(yintercept = hline),
    inherit.aes = FALSE,
    colour = "#9b0303",
    linetype = "dashed",
    linewidth = 0.6
  ) +
  
  facet_wrap(~trait_name_plot_short, ncol = 2, scales = "free_y") +
  
  # KEY FIX: make x labels readable
  coord_flip() +
  
  labs(x = NULL, y = "Number of credible set variants") +
  
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    
    # improves separation between facets
    panel.spacing = unit(1, "lines"),
    
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    
    # since we flipped coords, this controls epitest labels
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

p1



##################################################
## 3.2 Mean CS size
##################################################

# one hline value per trait
hline.df <- df.f %>%
  filter(epitest == "standard") %>%
  select(trait_name_plot_short, hline = mean_snps.inCS)


p2 <- ggplot(df.f, aes(x = annot_id, y = mean_snps.inCS)) +
  
  geom_point(aes(color=method_group), size = 2.5) +
  
  scale_color_manual(values = c(
    "standard" = "#9b0303",
    "selected"  = "#1f77b4",
    "other"    = "black"
  )) +
  
  geom_hline(
    data = hline.df,
    aes(yintercept = hline),
    inherit.aes = FALSE,
    colour = "#9b0303",
    linetype = "dashed",
    linewidth = 0.6
  ) +
  
  facet_wrap(~trait_name_plot_short, ncol = 1, scales = "free_y") +
  
  # KEY FIX: make x labels readable
  coord_flip() +
  
  labs(x = NULL, y = "Number of credible set variants") +
  
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    
    # improves separation between facets
    panel.spacing = unit(1, "lines"),
    
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    
    # since we flipped coords, this controls epitest labels
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

p2




##################################################
## 3.3 Number of high-confidence variants
##################################################

# one hline value per trait
hline.df <- df.f %>%
  filter(epitest == "standard") %>%
  select(trait_name_plot_short, hline = sum_snps_0.5)


p3 <- ggplot(df.f, aes(x = annot_id, y = sum_snps_0.5)) +
  
  geom_point(aes(color=method_group), size = 2.5) +
  
  scale_color_manual(values = c(
    "standard" = "#9b0303",
    "selected"  = "#1f77b4",
    "other"    = "black"
  )) +
  
  geom_hline(
    data = hline.df,
    aes(yintercept = hline),
    inherit.aes = FALSE,
    colour = "#9b0303",
    linetype = "dashed",
    linewidth = 0.6
  ) +
  
  facet_wrap(~trait_name_plot_short, ncol = 1, scales = "free_y") +
  
  # KEY FIX: make x labels readable
  coord_flip() +
  
  labs(x = NULL, y = "Number of credible set variants") +
  
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    
    # improves separation between facets
    panel.spacing = unit(1, "lines"),
    
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    
    # since we flipped coords, this controls epitest labels
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

p3




##################################################
## 3.4 All metrics
##################################################
snps.selected <- "promSNPs" # "allSNPs"
df.f <- df.alltraits[df.alltraits$selectedSNPS == snps.selected,]
df.f$annot_id <- factor(df.f$annot_id, levels = rev(c(1:16, "S")))
df.f$trait_name_plot_short <- factor(df.f$trait_name_plot_short, levels = traits.short)

df.f$method_group <- ifelse(df.f$epitest == "standard", "standard",
                            ifelse(df.f$epitest == "comb2DD", "selected", "other"))


df.long <- df.f %>%
  pivot_longer(
    cols = c(n_ldblocks_withCS, n_cs, mean_snps.inCS, sum_snps_0.5),
    names_to = "metric",
    values_to = "value"
  )
df.long$metric <- gsub("mean_snps.inCS", "Mean CS size", df.long$metric )
df.long$metric <- gsub("n_cs", "Num. CS", df.long$metric )
df.long$metric <- gsub("sum_snps_0.5", "Num. SNPs PIP>0.5", df.long$metric )
df.long$metric <- gsub("n_ldblocks_withCS", "Num. LD blocks with CS", df.long$metric )
df.long$metric <- factor(df.long$metric, levels = c("Num. LD blocks with CS", "Num. CS", "Mean CS size", "Num. SNPs PIP>0.5"))

hline.df <- df.long %>%
  filter(epitest == "standard") %>%
  select(trait_name_plot_short, metric, hline = value)

pp <- ggplot(df.long, aes(x = annot_id, y = value)) +
  
  geom_point(aes(color=method_group), size = 1.8) +
  
  scale_color_manual(values = c(
    "standard" = "#9A4444",
    "selected"  = "#12355B",
    "other"    = "grey60"
  )) +
  
  geom_hline(
    data = hline.df,
    aes(yintercept = hline),
    inherit.aes = FALSE,
    linetype = "dashed",
    colour = "#9A4444"
  ) +
  
  facet_grid(trait_name_plot_short ~ metric, scales = "free_y") +
  
  coord_flip() +
  
  theme_classic(base_size = 11) +
  
  theme(
    panel.border = element_rect( colour = "black", fill = NA,  linewidth = 0.6),    # Add panels borders
    panel.spacing.y =  unit(0.3, "lines"),
    panel.spacing.x =  unit(0.3, "lines"),
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 11),
    legend.position = "none"
  ) +
  
  labs(x = NULL, y = NULL)


#
pdf(paste0(outputdir.path, "Benchmarking_epitests_dotplot_nCS_CSsize_highconfSNPs_", snps.selected,".pdf"), width=13, height = 16)
print(pp)
dev.off()

png(paste0(outputdir.path, "Benchmarking_epitests_dotplot_nCS_CSsize_highconfSNPs_", snps.selected,".png"), width=3200, height = 4600, res = 300)
print(pp)
dev.off()
