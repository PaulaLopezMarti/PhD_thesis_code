# ==============================================================================
# Generate PolyFUN annotation files by intersecting functional genomic features
# with 1000 Genomes Project (EUR) SNPs.
#
# For each chromosome:
#   1. Load the corresponding 1KG BIM file.
#   2. Annotate SNPs with ChromHMM states across all cell types.
#   3. Create additional COO-specific enhancer annotations.
#   4. Annotate coding sequence (CDS) regions.
#   5. Export PolyFUN-compatible .annot and .l2.M files.
#
# Original script:
# /home/plopez/scratch/B_cell_roadmap_2/polyFUN/scripts/01_PolyFun_prepare_annotation_files_1KG_COO_original.R
# ==============================================================================


# ==============================================================================
# 0. Load required packages
# ==============================================================================
library("arrow") # for parquet files
library("data.table")
library("GenomicRanges")


# ==============================================================================
# 1. Global configuration
# ==============================================================================
chr_list <- as.character(1:22)
celltypes.sort <- c("HSC", "PreProB", "ProB", "PreB", "immtransB", "nB", "GCB", "memB", "PC", "CMP", "Mon", "nCD8")
bcells <- c("HSC", "PreProB", "ProB", "PreB", "immtransB", "nB", "GCB", "memB", "PC")
names_cols <- c("SNP","CHR", "BP", "BP1", "A1","A2") # Columns required by PolyFUN annotation files

# Definitions of Cell-of-Origin (COO) groups.
# Each COO is represented by a set of cell types whose enhancer landscape will be contrasted against the non-COO B-cells
coo.list <- list( "COO1" = c("PreProB", "ProB", "PreB"),
                  "COO2" = c("nB", "memB"),
                  "COO3" = c("GCB"),
                  "COO4" = c("PC"),
                  "COO5" = c("nB", "GCB", "memB"),
                  "COO6" = c("GCB","PC")
                  )
# ChromHMM states that will be used for polyfun
selected_states <- c("AEnh", "PoEnh", "PrEnh", "APromoter", "AEnhATAC", "PoEnhATAC", "PrEnhATAC", "APromoterATAC", "ATACpeaks", "PIRs")

# Prefix used for output annotation files
name_test <- "COOannot"


# ==============================================================================
# 2. Input/output paths
# ==============================================================================
# Local paths:
# outputdir <- paste0("~/scratch/B_cell_roadmap_2/polyFUN/data/annotations/", celltype)
# path_plinkfiles <- "~/IJC/genomes/1000genomes/1000genomes_GRCh38/plink_PID_EUR_nodup/"
# cds.file <- "~/scratch/db/genomes/GRCh38/annotations/CDS_GRCh38_104.bed"
# path_chromhmm_states <- "~/scratch/B_cell_roadmap_2/chromHMM/model_GRCh38_by_state"

# MN5 paths:
outputdir <- paste0("/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/polyFUN/data/annotations/", name_test)
path_plinkfiles <- "/gpfs/projects/bsc08/shared_projects/IJC_3Dchromatin/genomes/1000genomes/1000genomes_GRCh38/plink_PID_EUR_nodup/"
cds.file <- "/gpfs/scratch/bsc08/bsc008246/db/genomes/GRCh38/annotations/CDS_GRCh38_104.bed"
path_chromhmm_states <- "/gpfs/scratch/bsc08/bsc008246/B_cell_roadmap_2/chromHMM/model_GRCh38_by_state"
ifelse(!dir.exists(outputdir), dir.create(outputdir, recursive = T), FALSE)

# ==============================================================================
# 3. Load genome-wide coding sequence (CDS) annotation
# ==============================================================================
cds <- data.table::fread(cds.file)
names(cds) <- c("chr", "start", "end")
new_col <- paste0("CDS")
cds[,new_col] <- 1
setkey(cds, chr, start, end)


# ==============================================================================
# 4. Main chromosome processing loop
# ==============================================================================
# Process each chromosome independently. This minimizes memory usage because
# BIM files are the largest inputs.
# ==============================================================================

for(chr in chr_list){
  
  print(chr)
  
  ## 4.1 Read bim files and format them to data table ----
  bim_file <- paste0(path_plinkfiles, "ALL.chr", chr, ".phase3.genotypes.PIDs.EUR.bim")
  bim <- data.table::fread(bim_file)
  names(bim) <- c("CHR", "SNP", "CM", "BP", "A1", "A2")
  bim$CM <- NULL
  bim$CHR <- as.character(bim$CHR)
  bim$BP1 <- bim$BP #foverlaps() requries start-end coordinates
  bim <- bim[,..names_cols]
  setkey(bim, CHR, BP, BP1)
  
  count <- 0
  
  ## 4.2 Overlap with cell type annotations ----
  # Now, for each cell type and chromatin state, we (1) load the genomic intervals, (2) overlap with SNP positions and create a binary annotation column
  # (3) merge results of all cell types and annotations.
  
  for(celltype in celltypes.sort){ # to modify: no need to add non-B cell annotations, we are not using them

    files_celltype <- list.files(path_chromhmm_states, pattern=paste0('^', celltype, '_'), full.names = T) # chromhmm data

    file_num <- 1

    for(file_num in 1:length(files_celltype)){
      
      file_selected <- files_celltype[file_num]
      annot_file <- strsplit(basename(file_selected), '_')[[1]][2]
      annot_file <- gsub(".bed","", annot_file)
      
      if(annot_file %in% selected_states){ # only use chromhmm data if they are one of the selected annotations
        print(annot_file)
        df <- data.table::fread(file_selected)
        names(df) <- c("chr", "start", "end")
        new_col <- paste0(celltype, '_', annot_file)
        df[,new_col] <- 1
        setkey(df, chr, start, end)
        
        # overlap bim file with annotation
        bim.overlap <- foverlaps(bim, df, by.x = c("CHR", "BP", "BP1"), by.y=c("chr", "start", "end"), nomatch=NA)
        bim.overlap[,c("start", "end")] <- NULL # remove non-necessary columns
        bim.overlap[is.na(get(new_col)), (new_col) := 0] # replace NA with 0 for not matching SNPs
        setcolorder(bim.overlap,c(names_cols, new_col))
        
        # merge all annotations and cell types in bim.final table
        if(count == 0){
          bim.final <- bim.overlap
        } else { bim.final <- merge(bim.final, bim.overlap, by=names_cols)}
        count <- 1
      }
    }
  }
  
  
  ## 4.3 Generate COO-specific enhancer annotation ----
  # For each COO set, we are adding COO specific information: all the B cells enhancers not defined as enhancers in COO
  for(coo in names(coo.list)){
    
    coo.celltypes <- coo.list[[coo]] # COO cell types
    non.coo.celltypes <- setdiff(bcells, coo.celltypes) # non-COO B cell types
  
    # Merge enhancers (AEnh, PoEnh, PrEnh) from COO cell types 
    count.coo = 0
    for(x in 1:length(coo.celltypes)){
      celltype <- coo.celltypes[x]
      print(celltype)
      
      for(annot_selected in c("AEnh","PoEnh", "PrEnh")){
        file_selected <- paste0(path_chromhmm_states, "/", celltype, "_", annot_selected, ".bed")
        annot_file <- strsplit(basename(file_selected), '_')[[1]][2]
        annot_file <- gsub(".bed","", annot_file)
        print(annot_file)
        df <- data.table::fread(file_selected)
        names(df) <- c("chr", "start", "end")
        
        if(count.coo==0){coo.enh = df} else{coo.enh <- rbind(coo.enh,df)}
        count.coo = 1
      }
    }
    
    coo.enh.gr <- makeGRangesFromDataFrame(coo.enh)
    coo.enh.gr <- reduce(coo.enh.gr)
    
    # Get non-COO enhancers: merge AEnhATAC from all remaining B-cells
    for(x in 1:length(non.coo.celltypes)){
      
      celltype <- non.coo.celltypes[x]
      print(celltype)

      file_selected <- paste0(path_chromhmm_states, "/", celltype, "_AEnhATAC.bed")
          
      annot_file <- strsplit(basename(file_selected), '_')[[1]][2]
      annot_file <- gsub(".bed","", annot_file)
      print(annot_file)
      df <- data.table::fread(file_selected)
      names(df) <- c("chr", "start", "end")
      
      if(x==1){non.coo.enh = df} else{non.coo.enh <- rbind(non.coo.enh,df)}
    }
      
    non.coo.enh.gr <- makeGRangesFromDataFrame(non.coo.enh)
    non.coo.enh.gr <- reduce(non.coo.enh.gr)
    
    # Remove enhancers shared with COO
    non.coo.nooverlap.gr <- setdiff(non.coo.enh.gr, coo.enh.gr, ignore.strand=T)
    non.coo.nooverlap.gr <- reduce(non.coo.nooverlap.gr)
    
    non.coo.nooverlap <- data.frame(non.coo.nooverlap.gr)
    non.coo.nooverlap <- non.coo.nooverlap[,c("seqnames", "start", "end")]
    
    new_col <- paste0(coo, '_NoCOOAEnhATAC')
    non.coo.nooverlap[,new_col] <- 1
    non.coo.nooverlap <- data.table(non.coo.nooverlap)
    setkey(non.coo.nooverlap, seqnames, start, end)
            
    # overlap bim file with annotation
    bim.overlap <- foverlaps(bim, non.coo.nooverlap, by.x = c("CHR", "BP", "BP1"), by.y=c("seqnames", "start", "end"), nomatch=NA)
    bim.overlap[,c("start", "end")] <- NULL
    bim.overlap[is.na(get(new_col)), (new_col) := 0] # replace NA with 0
    setcolorder(bim.overlap,c(names_cols, new_col))
    
    # Merge with full annotation SNP table
    bim.final <- merge(bim.final, bim.overlap, by=names_cols)
    }

  
  ## 4.4 Add coding-sequence annotation ----
  bim.overlap <- foverlaps(bim, cds, by.x = c("CHR", "BP", "BP1"), by.y=c("chr", "start", "end"), nomatch=NA)
  bim.overlap[,c("start", "end")] <- NULL
  new_col <- "CDS"
  bim.overlap[is.na(get(new_col)), (new_col) := 0] # replace NA with 0
  setcolorder(bim.overlap,c(names_cols, new_col))
  bim.final <- merge(bim.final, bim.overlap, by=names_cols)
  
  ## 4.5 Add mandatory PolyFUN baseline annotation ----
  # This corresponds to all SNPs
  bim.final$Base <- 1
  bim.final[,c("BP1")] <- NULL
  
  
  
  ## 4.6 Export annotation files ----
  ## Ouput formats: 
  ## (1) *.annot.parquet -> Standard compressed polyfun annotation file
  write_parquet(bim.final, paste0(outputdir, '/', name_test, '.', chr,'.annot.parquet'))
  
  ## (2) *.annot -> same as *.annot.parquet but not compressed. It is not required by the software, so not saving it could save storage space, but it was useful for internal checks. 
  ##            Consider to remove it from the pipeline
  data.table::fwrite(bim.final, file=paste0(outputdir, '/', name_test, '.', chr, '.annot'), sep='\t')
  

  ## (3) *.l2.M -> annotation counts required by PolyFun. Whitespace-delimeted. Each value equals the number of SNPs assigned to the annotation
  l2 <- bim.final[,-c("SNP", "CHR", "BP", "A1", "A2")]
  sum_l2 <- colSums(l2)
  data.table::fwrite(t(data.table(sum_l2)), file=paste0(outputdir, '/', name_test, '.', chr,'.l2.M'), sep=' ', col.names = F)
    
}


# Final table structure:
# SNP metadata columns:
#   SNP, CHR, BP, A1, A2
#
# Binary annotation columns:
#   <celltype>_<state>
#   COO*_NoCOOAEnhATAC
#   CDS
#   Base
#
# Values:
#   1 = SNP overlaps annotation
#   0 = SNP does not overlap annotation