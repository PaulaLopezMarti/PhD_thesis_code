---
editor_options: 
  markdown: 
    wrap: sentence
---

# Chapter 3: Epigenetic fine-mapping

Here I present the code used for the epigenetic and normal fine-mapping and their comparison.

------------------------------------------------------------------------

## Software requirements

-   plink 2.0
-   plink 1.9 *#check if all can be done in one version of plink*
-   [PolyFun](https://github.com/omerwe/polyfun)
-   [multiCOGS](https://github.com/FunctionalGeneControl/multiCOGS)
-   R and R packages:
    -   [SusieR](https://stephenslab.github.io/susieR/)
    -   [bigsnpr](https://github.com/privefl/bigsnpr): To match SNPs with reference panel
    -   [arrow](https://arrow.apache.org/docs/r/): To handle parquet file format
    -   data.table
    -   GenomicRanges

------------------------------------------------------------------------

## Summary of workflow

1.  Data preparation

2.  Compute prior causal probabilities using PolyFUN

    1.  Create your own annotation files

    2.  Compute LD-scores for these annotations

    3.  Match GWAS summary statistics to the reference panel

3.  Fine-mapping with SusieR

------------------------------------------------------------------------

## Data requirements

To run epi-multiCOGS you need three types of input data:

1.  **Functional annotations:** the set of genomic regions with the functional annotations you want to use to perform the epigenetic fine-mapping.
    This can consist on enhancer regions, promoter regions, coding sequence, etc. of your cell type of interest.

    They must be provided in BED format (chr, start, end).

2.  **GWAS summary statistics:** it has to include, for your set of SNPs: SNP ID, chr, BP, p-value, z-score.
    Also need to know the population in which the GWAS is performed and the number of cases and controls.

## 1. Data preparation

### 1.2 Download and prepare reference panel

As we are using summary statistics instead of the genotyped individual, we need to use a reference panel in different steps of the pipeline.
Specifically, it is required for the computation of LD-scores with PolyFun and the computation of LD matrices.
In this case, we will work with European population from 1000 Genomes phase 3.
Even if this option is not recommended due to the small sample size (recommended reference panel of \>3000 individuals, see PolyFun manual), it is sometimes the only publicly available panel.
In case you have access to a bigger reference panel, it is advisable to use it.
In case you have your genotyped individuals, check PolyFun manual to run with this option.

1000 Genomes phase 3 mapped to GRCh38 VCF files can be downloaded from <http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/>

In order to work with 1000 Genomes data and to compute LD matrices, we need to convert VCF files into plink files format (bed, bim and fam).
For this, we can use the plink command --vcf.

```         
for i in {1..22}; do plink --vcf ALL.chr1.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz --make-bed --out ALL.chr${i}.phase3.genotypes; done
```

Moreover, we need to filter by the genetic ancestry of our GWAS data; in this case, European individuals.
In case you are working with a different genetic ancestry GWAS, take this into consideration.
We can filter by individuals using plink *--keep* command.

```         
for i in {1..22}; do plink --bfile chr${i}.phase3.genotypes --keep populations/pop_EUR.txt --make-bed --out ALL.chr${i}.phase3.genotypes.EUR; done
```

Where pop_EUR.txt is a file that contain the IDs of the EUR individuals of 1000 Genomes.

As we are using positional IDs instead of rsID, we will annotate the SNPs ID based on positional ID (chr_pos).
We will conserve only SNPs with the argument *--snps-only* that excludes variants with multi-character alleles codes.
Moreover, we will remove duplicated positional IDs and conserve only the first instance of duplicated SNPs ID with *--rm-dup* force-first.

```         
for i in {1..22}; do plink2 --bfile ALL.chr${i}.phase3.genotypes.EUR --set-all-var-ids @_# --snps-only --make-bed --rm-dup force-first --out ALL.chr${i}.phase3.genotypes.PIDs.EUR; done
```

*Generate AF files!!!*

### 1.3 Approximate independent LD blocks

We used the independent LD blocks computed for european population based on GRCh38 1000 Genomes data.
They can be downloaded from [here](https://github.com/jmacdon/LDblocks_GRCh38).
This corresponds to a total 1361 LDblocks in chromosomes 1 to 22, with an average size of approx 2Mb (minimum: 17kb; max: 28Mb).

------------------------------------------------------------------------

## 2. Compute prior causal probabilities using PolyFUN

Fine-mapping, that consists in two main parts: first, computing prior causal probabilities with PolyFun; and second, performing the multivariate fine-mapping with SusieR.

*Steps 2.1 and 2.2 are based on the tutorial provided in the section [2. Using and creating functional annotations](https://github.com/omerwe/polyfun/wiki/2.-Using-and-creating-functional-annotations) of their Github repository.*

### 2.1 Overlap reference panel SNPs with the functional annotations

To compute LD scores on our custom annotations, PolyFUN requires a binary annotation files indicating, for each SNP, whether it overlaps a given functional annotation.
We evaluated 16 different annotation combination.

Briefly, the tested annotations included:

-   Active promoters and enhancers (active, primed and poised) from the corresponding cell of origin (COO), either overlapping or not overlapping ATAC-seq peaks.

-   Open chromatin regions (ATAC-seq peaks) in the COO.

-   Open active enhancers from non-COO B cells that do not overlap COO enhancers (active, primed or poised).

-   Promoter-interacting regions (PIRs).

-   Coding sequence: CDS regions defined in Ensembl release 104.

*TABLE xxxx*

Annotation files for all the tested annotations were generated using the script `01_PolyFun_prepare_annotation.R`.
For each chromosome:

1.  Load the 1000 Genomes SNPs .bim file.

2.  Annotate SNPs with the selected chromHMM states across all cell types and add COO-specific annotation as well as coding regions.
    It also adds a Baseline annotation (all SNPs) required by polyFUN processing.

3.  Export data.

It generates three main output files:

1.  **\*.annot.parquet** - Standard compressed PolyFUN annotation file containing SNP information (SNP, CHR, BP, A1, A2) and binnary annotation columns indicating whether each SNP overlaps a given annotaiton (1) or not (0).

2.  **\*.annnot** - Uncompressed version of \*.annot.parquet.
    It is not required by PolyFUN, so not saving it could save storage space, but it was useful for internal checks.
    *Consider to remove it from the pipeline*

3.  **\*.l2.M** - Annotation counts required by PolyFUN during LD score computation, containing the number of SNPs assigned to each annotation.
    It has to be whitespace-delimited.

### 2.2 Compute LD-scores for annotations with the reference panel

Once we have generated our annotation files, the next step is to compute LD-scores for these annotations using our reference panel.

This is performed with the PolyFUN script `compute_ldscores.py`, which have to be run separately for each chromosome.
The script takes as input the plink .bed file from the reference panel and the corresponding annotation file generated in Step 2.1.

Specifically, we used the parameter `--ld-wind-kb 1000`, which defines a window size of 1Mb to compute LD-scores.
Under standard LDSC workflows, LD windows are tipically defined using genetic distance (cM).
However, the 1000 Genomes reference files used here do not contain the genetic map coordinates, so we decided to use a physical-distance window with a size of 1Mb as [recommended by PolyFUN authors](https://github.com/omerwe/polyfun/issues/103).

In addition ot annotation-specific LD scores, PolyFUN also requires a set of LD-score weights for regression.
These can be generated by running the same script withoud providing an annotation file (we are adding the *\_w* suffix to the output file name).

The resulting LD-scores files must end with the suffix `.l2.ldscore.parquet`, as required by PolyFUN for downstream processing.
In addition, they must be saved in the same directory as the annotation files, since PolyFUN expects annotation and LD-score outputs to be co-located for automatic file matching and workflow compatibility.

Example of commands for chromosome 22:

``` bash
# LD-scores for custom annotations
python compute_ldscores.py \
--bfile 1000genomes_GRCh38/ALL.chr22.phase3.genotypes.PIDs.EUR \
--annot data/annotations/COOannot.22.annot.parquet \
--ld-wind-kb 1000 \
--out data/annotations/COOannot.22.l2.ldscore.parquet \

# Regression weights (without annotations)
python compute_ldscores.py \
--bfile 1000genomes_GRCh38/ALL.chr22.phase3.genotypes.PIDs.EUR \
--ld-wind-kb 1000 \
--out data/annotations/COOannot_w.22.l2.ldscore.parquet \
```

The script `02_PolyFun_prepare_commands_LDScore.R` generates all commands for computing LD-scores and regression weights for all chromsomes, and print them to a bash script.
Since each chromosome can be processed independently, these commands can be submitted to run in parallel on an HPC cluster.
Runtime typically ranges from approximately 3-15 minutes per chromosome.

### 2.3 Match GWAS summary statistics to the reference panel

We can only perform fine-mapping for variants that are present in both the GWAS summary statisticis and the LD reference panel.
Therefore, before running PolyFUN, it is essential to harmonize the GWAS data with the reference panel and ensure that alleles are consitently defined.

#### Overview

This step performs the following tasks:

1.  Standardizes summary statistics column names to PolyFUN conventions.

2.  Matches GWAS variants to the reference panel.

    1.  Adds allele frequency information from the reference panel.

    2.  Harmonizes allele orientation between the GWAS and the reference panel.

    3.  Identifies and evaluates strand-flipped variants.

3.  Generates the final PolyFUN-compatible summary statistics file.

#### Standardize summary statistics

PolyFUN requires summary statistics in a specific format.
Although the PolyFUN script `munge_polyfun_sumstats.py` can recognize many commonly used column names, some columns still need to be standardized beforehand.
In particular, it is critical that allele columns are correctly defined: A1 correspond to the effect allele and A2 to the non-effect (reference) allele.

#### Match variants to the reference panel

The GWAS summary statistics are matched against the 1000 Genomes EUR reference panel:

-   Variants absent in the reference panel are excluded.

-   Allele frequency (AF) information is added.

-   Alleles are harmonized.
    For this, we used the function `bigsnpr::snp_match()`, that detects and resolves alleles mismatches between the GWAS summary statistics and the reference panel.
    Specifically, we matched variants using their positional identifiers (`chr_position`).

    One common type of mismatch occurs when the **effect and non-effect** **alleles are swapped** between datasets (A1 ↔ A2).
    In these cases, `snp_match()` swap alleles and automatically updates the sign of `BETA` column, to preserve the correct direction of effect.
    However, it does not adjust Z-score direction.
    Therefore, for those variants flagged as reversed (`_REV_ = TRUE`), Z-score was corrected manually by multiplying it by -1.

#### Strand flips and ambiguous SNPs

Another type of allele mismatching are those SNPs that require strand flipping.
This happens when alleles are reported on oppsite DNA strands in the GWAS and the reference panel (A ↔ T and C ↔ G).
This is usually uncommon because most GWAS datasets are aligned to the same reference genome used by 1000 Genomes.
However, when strand flips are present, special care is needed for ambiguous SNPs (A/T; C/G), as it is impossible to determine the correct strand using allele information alone.

To avoid allele mismatching, a common approach is to remove all ambiguous SNPs.
However, this can lead to the loss of a substantial number of variants.
To avoid unnecessary data loss, we first evaluated the number of SNPs requiring strand flips across all GWAS datasets.
Because only a very small number of variants required strand flipping, we decided to remove the few strand-flipped variants while retaining all ambiguous SNPs.

The evaluation of strand flips, as well as other steps of GWAS-reference panel matching, is performed separately in two scripts (`03A_PolyFun_check_sumstats_strand_flip.R` and `03B_PolyFun_check_sumstats_strand_flip_report.R`) and it is described below.

### Scripts

This step requires to use three scripts:

1.  **`03A_PolyFun_check_sumstats_strand_flip.R`**: it generates harmonized summary statistics files that are used to evaluate the impact of the different harmonization and quality control steps.

    This script is executed independently for each GWAS dataset, using the trait name as argument.
    Example:

    ``` bash
    03A_PolyFun_check_sumstats_strand_flip.R BALL
    ```

    The helper script `03D_PolyFun_format_sumstats_prepare_commands.R` was used to generate the corresponding commands and the slurm job files to run on the HPC cluster.

2.  **`03B_PolyFun_check_sumstats_strand_flip_report.R`**: it summarizes the results generated by `03A` across all datasets and creates a report table describing the impact of each harmonization step.
    This script is executed after all datasets have been processed with the `03A` script.

    The report includes, for different p-value thresholds:

    -   Total number of variants in the GWAS.

    -   Number of variants present in the reference panel.

    -   Number of reversed variants (allele swap).

    -   Number of strand-flipped variants.

    -   Number of ambiguous SNPs (A/T and C/G).

    In addition, the report also summarizes the number of variants that would be excluded under different minor allele frequency thresholds (0, 0.01 and 0.05).
    SNPs with very low MAF can lead to unreliable results.
    As the script `munge_polyfun_sumstats.py` used in the next step includes the option of filtering by MAF, we included this information in the report to guide in the decision making of whether a MAF filter should be applied and, if so, which threshold is most appropriate.

    | Trait | SNPs selected | Number SNPs | Included in 1KG | Flip strand | Ambiguous | Swap alleles | TotalSNPs conserved (removing ambiguous) | TotalSNPs conserved (with ambiguous, removing strand flip) |
    |----|----|----|----|----|----|----|----|----|
    |  |  |  |  |  |  |  |  |  |
    |  |  |  |  |  |  |  |  |  |
    |  |  |  |  |  |  |  |  |  |

3.  **`03C_PolyFun_format_sumstats_for_munge.R`**: it performs the final harmonization and formatting of summary statistics, as explained above.
    It also reads output of `03A` to remove variants requiring strand flipping.

    Again, this script is executed independently for each GWAS dataset, and the helper script `03D_PolyFun_format_sumstats_prepare_commands.R` was used to generate the corresponding commands and slurm job files.

    ### Output

    The output file of this step (***trait*****\_polyFUNformattedClean.tsv**) contains the harmonized summary statistics that will be used in all downstream analysis, including polyFUN and fine-mapping.

### 2.4. Munge summary statistics

PolyFUN provides the script `munge_polyfun_sumstats.py`, that creates a munged summary statistics adequate for running PolyFUN, in parquet format.

*More detailed information about this step and the next one (2.4 and 2.5) can be found in the [PolyFun approach 2](https://github.com/omerwe/polyfun/wiki/1.-Computing-prior-causal-probabilities-with-PolyFun#polyfun-approach-2-computing-prior-causal-probabilities-via-an-l2-regularized-extension-of-s-ldsc) section of their GitHub tutorial.*

An example of command used:

``` bash
python munge_polyfun_sumstats.py \ 
--sumstats data/sumstats/BALL_polyFUNformattedClean.tsv \
--n 3473 \
--min-info 0 \
--min-maf 0.01 \
--out data/sumstats/BALL_munged.parquet
```

The argument `-n` indicates the sample size of the GWAS.
For those datasets where the sample size information was available as a column in the summary statistics file, we omitted this parameter.
The argument `--min-info` was set to 0 to not filter by INFO (this step was already performed in the curation of the summary statistics).
We used a minimum MAF of 0.01 (`--min-maf`) to remove rare variants.

The script `04_PolyFun_prepare_commands_mungeSumstats.R` was used to generate the commands to munge all datasets, as well as their correspoding slurm jobs.

### 2.5. Run PolyFun

Now we have the LD-scores from our annotation computed in 2.2 and the GWAS summary statistics in polyfun format from 2.4.
Final step is to run PolyFun using an L2-regularized extension of S-LDSC.
For this, we have to used the polyfun.py script.

Example of command used:

``` bash
python polyfun.py \
--compute-h2-L2 \
--no-partitions \
--sumstats data/sumstats/BALL_munged.parquet \
--ref-ld-chr data/annotations/COOannot/COOannot. \
--w-ld-chr data/annotations/COOannot/COOannot_w. \
--anno annot1,annot2,annot3,Base \
--output-prefix results/BALL \
```

As input, it takes the munged summary statistics created in the previous step (`--sumstats`).
The files passed to `--ref-ld-chr` and `--w-ld-chr` correspond to the LD-scores and weights files computed in step 2.2 (PolyFUN automatically looks for files with that prefix and the *.l2.ldscore.parquet* suffix).
It is also important that the annotations files (.*annot*) and *.l2.M* file are located also in the same directory as the LD-scores and weights.
In our case, we are storing all these files in the *data/annotations* directory.

The argument `--anno` reads a comma-delimited list of annotations to use to compute the per-SNP h2.
Here, we have to provide the annotations that we decided to use for each GWAS dataset.
It is important to include the *Base* annotation and that the annotation list that we provide follows the same order as the in the LD-scores files (otherwise, it crashes).

The script `05_PolyFun_prepare_commands_runPolyFun.R` was used to generate the commands to run PolyFun for all the annotation combinations tested in all datasets.

#### **Output**

As output, we will obtain two files per chromosome.
The one we are interested in is `*.snpvar_ridge_constrained.gz`, which contains the truncated per-SNP heritabilities (***SNPVAR*** column) that can be used as prior causal probabilities for fine-mapping.

------------------------------------------------------------------------

## 3. Fine-mapping with SuSiER

### 3.1 Extract LD blocks with promising SNPs

### 3.2 Compute LD matrices

### 3.3 Run SusieR

### 3.4 Fine-mapping QC

------------------------------------------------------------------------

## 4. Fine-mapping analysis

### 4.1 Fine-mapping evaluation: standard vs. epigenetic fine-mapping

This section will include:

-   Compare number of SNPs by PIP: standard vs. epigenetics

-   Compare CS size: standard vs. epigenetics

-   Overlap with RE / chromHMM from Roadmap

### 4.2 Analysis of disrupted TF motif

This will include how to run motifbreakR in local, compute TF score and generate the plots presentend in the paper/thesis.

### 4.3 B-ALL GWAS: comparison with public data

-   Processing of B-ALL data

-   Comparison with SNPs

-   Cell lines?
    analysis and comparison
