---
editor_options: 
  markdown: 
    wrap: sentence
  chunk_output_type: inline
---

# Epi-multiCOGS

Here I present the code of epi-multiCOGS used for the thesis.
Specifically, this part corresponds to the functional fine-mapping.

------------------------------------------------------------------------

## Software requirements

-   plink 2.0
-   plink 1.9 *#check if all can be done in one version of plink*
-   [PolyFun](https://github.com/omerwe/polyfun)
-   [multiCOGS](https://github.com/FunctionalGeneControl/multiCOGS)
-   R and R packages:
    -   [SusieR](https://stephenslab.github.io/susieR/) (0.14.2)
    -   [bigsnpr](https://github.com/privefl/bigsnpr): To match SNPs with reference panel
    -   [arrow](https://arrow.apache.org/docs/r/): To handle parquet file format
    -   data.table
    -   GenomicRanges
    -   Rfast
    -   cowplot / patchwork

------------------------------------------------------------------------

## Summary of workflow

1.  Data preparation

2.  Compute prior causal probabilities using PolyFUN

    1.  Create your own annotation files

    2.  Compute LD-scores for these annotations

    3.  Match GWAS summary statistics to the reference panel

    4.  Munge summary statistics

    5.  Run PolyFun using an L2-regularized extension of S-LDSC

3.  Fine-mapping with SusieR

    1.  Extract LD blocks with promising SNPs to perform fine-mapping

    2.  Compute LD matrices from the reference panel

    3.  Conduct diagnostics and fine-mapping with susieR

    4.  Quality control of fine-mapping results

    5.  Create final fine-mapping results

    6.  Benchmarking epigenetic combinations

------------------------------------------------------------------------

## Data requirements

To run epi-multiCOGS you need three types of input data:

1.  **Functional annotations:** the set of genomic regions with the functional annotations you want to use to perform the epigenetic fine-mapping.
    This can consist on enhancer regions, promoter regions, coding sequence, etc. of your cell type of interest.

    They must be provided in BED format (chr, start, end).

2.  **GWAS summary statistics:** it has to include, for your set of SNPs: SNP ID, chr, BP, p-value, z-score.
    Also need to know the population in which the GWAS is performed and the number of cases and controls.

## 1. Data preparation

### 1.1 Download and prepare reference panel

As we are using summary statistics instead of the genotyped individuals, we need to use a reference panel in different steps of the pipeline.
Specifically, it is required for the computation of LD-scores with PolyFun and the computation of LD matrices.
In this case, we will work with European population from 1000 Genomes phase 3.
Even if this option is not recommended due to the small sample size (recommended reference panel of \>3000 individuals, see PolyFun manual), it is sometimes the only publicly available panel.
In case you have access to a bigger reference panel, it is advisable to use it.
In case you have your genotyped individuals, check PolyFun manual to run with this option.

1000 Genomes phase 3 mapped to GRCh38 VCF files can be downloaded from <http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/>

In order to work with 1000 Genomes data and to compute LD matrices, we need to convert VCF files into plink files format (bed, bim and fam).
For this, we can use the PLINK command `--vcf`.

``` bash
for i in {1..22}; do plink --vcf ALL.chr1.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz --make-bed --out ALL.chr${i}.phase3.genotypes; done
```

Moreover, we need to filter by the genetic ancestry of our GWAS data; in this case, European individuals.
In case you are working with a different genetic ancestry GWAS, take this into consideration.
We can filter by individuals using `--keep` command.

``` bash
for i in {1..22}; do plink --bfile chr${i}.phase3.genotypes --keep populations/pop_EUR.txt --make-bed --out ALL.chr${i}.phase3.genotypes.EUR; done
```

Where pop_EUR.txt is a file that contain the IDs of the EUR individuals of 1000 Genomes.

As we are using positional IDs instead of rsID, we will annotate the SNPs ID based on positional identifiers (`chr_position`).
We will conserve only SNPs with the argument `--snps-only` that excludes variants with multi-character alleles codes.
Moreover, we will remove duplicated positional IDs and conserve only the first instance of duplicated SNPs ID with `--rm-dup force-first`.

``` bash
for i in {1..22}; do plink2 --bfile ALL.chr${i}.phase3.genotypes.EUR --set-all-var-ids @_# --snps-only --make-bed --rm-dup force-first --out ALL.chr${i}.phase3.genotypes.PIDs.EUR; done
```

We also need to compute allele frequencies for the curated set of 1KG variants.
We can do this with the command `--freq` and then merge all the chromosomes in a single file (*EUR.chrALL.GRCh38.AlleleFrequency.bed*).

``` bash
for i in {1..22}; do plink --bfile /gpfs/projects/bsc08/shared_projects/IJC_3Dchromatin/genomes/1000genomes/1000genomes_GRCh38/plink_PID_EUR_nodup/ALL.chr${i}.phase3.genotypes.PIDs.EUR --freq --out /gpfs/projects/bsc08/shared_projects/IJC_3Dchromatin/genomes/1000genomes/1000genomes_GRCh38/EUR_AF_PID_nodup/EUR.chr${i}.GRCh38.AlleleFrequency; done
```

### 1.2 Approximate independent LD blocks

We used the independent LD blocks computed for european population based on GRCh38 1000 Genomes data.
They can be downloaded from [here](https://github.com/jmacdon/LDblocks_GRCh38).
This corresponds to a total 1361 LDblocks in chromosomes 1 to 22, with an average size of approx 2Mb (minimum: 17kb; max: 28Mb).

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

#### Strand flips and ambiguous SNPs

One common type of mismatch occurs when the **effect and non-effect** **alleles are swapped** between datasets (A1 ↔ A2).
In these cases, `snp_match()` swap alleles and automatically updates the sign of `BETA` column, to preserve the correct direction of effect.
However, it does not adjust Z-score direction.
Therefore, for those variants flagged as reversed (`_REV_ = TRUE`), Z-score was corrected manually by multiplying it by -1.

Another type of allele mismatching are those SNPs that require **strand flipping**.
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

    The report (`Report_alltraits_summary_reference_matching.tsv`) includes, for different p-value thresholds:

    -   Total number of variants in the GWAS.

    -   Number of variants present in the reference panel.

    -   Number of reversed variants (allele swap).

    -   Number of strand-flipped variants.

    -   Number of ambiguous SNPs (A/T and C/G).

    In addition, the report also summarizes the number of variants that would be excluded under different minor allele frequency thresholds (0, 0.01 and 0.05).
    SNPs with very low MAF can lead to unreliable results.
    As the script `munge_polyfun_sumstats.py` used in the next step includes the option of filtering by MAF, we included this information in the report to guide in the decision making of whether a MAF filter should be applied and, if so, which threshold is most appropriate.

3.  **`03C_PolyFun_format_sumstats_for_munge.R`**: it performs the final harmonization and formatting of summary statistics, as explained above.
    It also reads output of `03A` to remove variants requiring strand flipping.

    Again, this script is executed independently for each GWAS dataset, and the helper script `03D_PolyFun_format_sumstats_prepare_commands.R` was used to generate the corresponding commands and slurm job files.

    ### Output

    The types of output files of this step are:

    1.  **Report_alltraits_summary_reference_matching.tsv** QC report

    2.  **\*\_polyFUNformattedClean.tsv** contains the harmonized summary statistics that will be used in all downstream analysis, including polyFUN and fine-mapping.

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
The argument `--min-info` was set to 0 to not filter by INFO (this step was already performed in the data curation of the summary statistics).
We used a MAF threshold of 0.01 (`--min-maf`) to remove rare variants.

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

## 3. Fine-mapping with SuSiER

### 3.1 Extract LD blocks containing association signals

COGS operates on independent LD blocks; therefore, fine-mapping must be performed using the same LD block structure downloaded in **Step 1.2**.

Because fine-mapping is computationally demanding, fine-mapping was restricted to LD blocks showing a minimum association signal.
Specifically, we selected LD blocks containing at least one variant with a p-value of \< 5e-5.

In addition, we excluded LD blocks containing the MHC region (GRCh38 MHC region from NCBI: chr6:28510120-33480577; LDblocks removed: chr6:28390030-33935919).

For this, we used the script `01_SusieR_Extract_prom_LDblocks_fromGWAS.R`, which has to be executed independently for each GWAS dataset using the trait name as argument:

``` bash
Rscript 01_SusieR_Extract_prom_LDblocks_fromGWAS.R BALL
```

The harmonized summary statistics generated for PolyFun would be the ideal input for this step.
However, due to historical reasons, this script uses as input the non-harmonized version of the summary statistics and repeat some curation steps.
This does not affect the validity of the final results, as the actual fine-mapping procedure (Section 3.3) is performed using the fully harmonized summary statistics.

**Processing steps**

The script:

1.  Excludes variants not present in the 1KG reference panel.

2.  Excludes variants with MAF ≤ 0.01

3.  Identifies LD blocks containing at least one variant with p \< 5e-5

4.  Generates the files and commands required to compute LD matrices.

#### Output

For each selected LD block, the script generates:

-   **SNPs_list/LDblockID_snps.txt**: list of positional SNP identifiers that would be used to compute LD matrices.

-   **PLINK Commands and SLURM jobs files to compute LD matrices:** this is explained more in detail in the next step.

### 3.2 Compute LD matrices

For each selected LD block, LD matrices were computed from the 1000 Genomes EUR reference panel using PLINK, with the argument `--r square`, which indicates to compute pairwise correlation (Pearson's *r*) as a square matrix.
We used the `--extract` argument to indicate the list of SNPs to use for the LD matrix calculation (LDblockID_snps.txt, generated in the previous step).

A part from the LD matrices, we also need to generate the bed files resulting for each LD block (`--make-bed`).

Example of commands for one LD block:

``` bash
# Compute LD matrices
plink \
--bfile 1000genomes_GRCh38/ALL.chr18.phase3.genotypes.PIDs.EUR \
--r square \
--keep-allele-order \
--extract SNPs_list/LDblock1230_chr18.26873450.29293548_snps.txt \
--out plink_files/LDblock1230_chr18.26873450.29293548_1KG

# Generate bed files
plink \
--bfile 1000genomes_GRCh38/ALL.chr18.phase3.genotypes.PIDs.EUR \
--keep-allele-order \
--make-bed \
--extract SNPs_list/LDblock1230_chr18.26873450.29293548_snps.txt \
--out plink_files/LDblock1230_chr18.26873450.29293548_1KG
```

The `--keep-allele-order` flag is essential for this step.
It indicates that the original A1/A2 allele order must be preserved.
By default, PLINK would automatically reorder alleles so that the minos alleles become the effect allele.
If this occurs, the LD matrix would not match the orientation used in the summary statistics, leading summary statistics and reference panel inconsistencies and, therefore, wrong fine-mapping results.

### 3.3 Run SusieR

Fine-mapping was performed using the SuSiE-R package and the `susie_rss()` fucntion, which enables Bayesian fine-mapping directly from GWAS summary statistics and an external LD reference panel.

Specifically, we ran `susie_rss()` taking as input the GWAS Z-score and the LD correlation matrices computed with plink.

Example:

``` r
susie_rss(z = ss$Z, 
          R = ld.mat, 
          n = study.size, 
          L = L_value, 
          max_iter = n_iter, 
          coverage = coverage.value )
```

Where:

-   z = GWAS Z-scores

-   R = LD correlation matrix

-   n = GWAS sample size

-   L = maximum number of causal variants allowed within LD block.

-   max_iter = maximum number of IBSS iterations

-   coverage = desired credible set coverage

We set the max_iter = 100 and the coverage to 0.90 for all the fine-mapping tests.
Regarding the maximum number of causal variants, we ran twice: assuming 1 (L=1) or maximum 5 (L=5) causal variants.

#### Fine-mapping diagnostics

Before interpreting fine-mapping results, it is important to assess whether the GWAS summary statistics and the LD reference panel are consistent, as mismatches between them can lead to unreliable posterior probabilities and credible sets.

SuSiE-R provides several **diagnostic tools** for this purpose:

Lambda statistic

The function `estimate_s_rss()` estimates a regularization parameter (`lambda`) that measures the discrepancy between the observed GWAS statistics and the LD matrix.

``` r
estimate_s_rss(z = ss$Z, R = ld.mat, n = study.size)
```

Large values of lambda indicate potential inconsistencies, so ideally lambda value should be close to 0.

**Diagnostic plots**

The function `kriging_rss()` generates observerd vs. expected Z-scores diagnostic plots.
Ideally, variants should lie close to the diagonal, as strong deviations may indicate LD mismatches.

*More detailed information is available in susieR Github about [how to run susie_rss](https://stephenslab.github.io/susieR/articles/finemapping_summary_statistics.html) or [perform diagnostic for fine-mapping with summary statistics](https://stephenslab.github.io/susieR/articles/susierss_diagnostic.html). Here, I aim to document the specific parameters and QC procedures in our study.*

------------------------------------------------------------------------

#### 1. Diagnostics and standard fine-mapping

Standard fine-mapping (without functional priors) and diagnostic analyses were performed using: `03B_SusieR_run_susier_byblock_v2.R`

This script processes one LD block at a time and automatically runs SuSiE-R assuming either:

-   **L=1** causal variants

-   **L=5** causal variants

It can also generates diagnostic plots and compute lambda statistics.

**Script Arguments TO MODIFY**

It has the following parameters:

-   **-t, --trait** Trait ID

-   **-b, --block** LD block ID

-   **-c, --cs** CS threshold

-   **-i, --iter** Max num iterations

-   **-d, --dplots** Flag specifying whether to generate diagnostic plots (observed vs. expected Z-score).

-   **-l, --lambda** Flag specifying whether to compute lambda values.

-   **-s, --savedata** Flag specifying whether to save LD blocks clean data.

-   **--savefinemapping** Flag specifying whether to save finemapping results.

-   **-p, --ppiplots** Flag specifying whether to save posterior probabilities plots.

Specifically, we ran for each LD block in each GWAS dataset the command for the diagnostic plots and the standard fine-mapping, setting up a maximum of 100 IBSS algorithm iterations and a CS coverage of 0.90.

``` bash
# Diagnostics
Rscript 03_SusieR_run_susier_byblock_v2.R \
--trait BALL \
--block LDblock921 \
--cs 90 \
--iter 100 \
--dplots \
--lambda \
--savedata

# Perform standard fine-mapping
Rscript 03_SusieR_run_susier_byblock_v2.R \
--trait BALL \
--block LDblock921 \
--cs 90 \
--iter 100 \
--savefinemapping \
--ppiplots \
--savedata
```

Diagnostics and fine-mapping were run separately because generating diagnostic plots is substantially more computationally intensive than fitting the fine-mapping model itself.

The helper script `03A_SusieR_run_susier_prepare_commands.R` generates all commands and SLURM job files to run diagnostics and standard fine-mapping across all GWAS datasets.

#### 2. PolyFun-informed fine-mapping

To incorporate functional information into fine-mapping, we used the prior causal probabilities computed by PolyFun.
These priors were supplied to `susie_rss()` throught the `prior_weights` arguments.

For this, we used the script `03C_SusieR_run_susier_polyfun.R`, which has the following arguments: *FIX THIS PART AFTER CLEANING THE SCRIPT*

-   **-t, --trait** Trait ID

-   **-c, --celltype** Celltype

-   **--test** PolyFun test name

-   **--cs** CS threshold

-   **-i, --iter** Max num iterations

-   **-s, --savefinemapping** Flag specifying whether to save finemapping results.

-   **-p, --ppiplots** Flag specifying whether to save posterior probabilities plots.

Unlike the script 03B, this script is executed per trait and annotation combination analysed by PolyFun, rather than per LD block.
Again, it also runs automatically for both L=1 and L=5 and stores results separetely.

The script also records LD blocks for which SuSiE-R raise a warning where **IBSS algorithm did not converge after 100 iterations**, as non-convergence may indicate unreliable fine-mapping results.

Example of command:

``` bash
# Perform epigenetic fine-mapping
Rscript 04_SusieR_run_susier_polyfun.R \
--trait BALL \
--celltype COO1 \
--test comb2DD \
--cs 90 \
--iter 100 \
--savefinemapping \
--ppiplots
```

The `--ppiplot` option generates plots comparing the results of standard and functionally informed fine-mapping, including Manhattan plots, Z-score plots and posterior inclusion probabilities (PIP) plots.
These plots will facilitate visual exploration and interpreation of fine-mapping results, as detailed in the next step.

Again, the helper script `03A_SusieR_run_susier_prepare_commands.R` was used to generate the commands and SLURM job files for all epigenetic fine-mapping analysis.

#### Output *PENDING*

-   **results/\<trait\>/finemapping/CV5\_\<epitest\>\_iter_100_cov90/\<LDblockID\>\_susieR_cs90_L5.tsv**: which contain the fine-mapping results for

-   plots....

-   ...

### 3.4 Fine-mapping QC

The diagnostics and fine-mapping outpus should be evaluated with the scripts:

-   `04A_Check_susier_lambda.R`

-   `04B_Check_bad_finemapping_MN5.R`

-   `04C_Check_finemapping_merge_reports.R`

#### **Output**

These scripts generate two main QC reports:

-   **trait_epitest_LDblocks_selection.tsv:** Contains QC results for each LD block within each trait and fine-mapping run.
    This file is used as input for downstream filtering and analysis.

-   **LDblocks_finemapping_QC_summary_alltraits.tsv:** Summarizes, across all traits and epigenetic tests, how many LD blocks pass each QC filter.
    Results are stratified by GWAS p-value thresholds.

    | Trait | pval_group | total_LDblocks | Lambda \> 0.4 | No converge | Not CS | CV5_has_CS | Good CV5 fine-mapping | Bad CV5 fine-mapping |
    |--------|--------|--------|--------|--------|--------|--------|--------|--------|
    | LYM | \<5e-8 | 3 | 1 | 0 | 0 | 2 | 2 | 0 |
    | LYM | 1e-6 - 5e-8 | 12 | 0 | 0 | 7 | 5 | 4 | 1 |
    | LYM | \>1e-6 | 102 | 2 | 0 | 98 | 2 | 2 | 0 |
    | LYM | Total | 115 | 2 | 0 | 105 | 8 | 8 | 0 |

#### QC classification of LD blocks

In the QC report, each LD block is annotated according to its fine-mapping performance and diagnostic metrics:

1.  **High LD/summary statistic inconsistency** (`lambda_over_0.4`): a high lambda value indicates inconsitency between LD matrix and summary statistics.
    These blocks are excluded from SuSie-R fine-mapping and downstream analysis.

2.  **IBSS algorithm non-convergence** (`no_convIBSS_CV1,` `no_convIBSS_CV5`).The IBSS algorithm did not converge after 100 iterations using either L=5 (CV5) or L=1 (CV1) configurations.

3.  **No association signal in credible sets** (`no_pass5e5_CV5`,`no_pass5e5_CV1`): The identified CS do not include any at least one variant with p-value \< 5e-5).
    This suggests suboptimal fine-mapping performance and is often consistent with IBSS non-convergence.

4.  **Model selection** (`use_cv1`): For LD blocks in categories (2) or (3) under CV5 indicates that fine-mapping was not optional and that CV1 results are preferred and should be used instead.

A part from using the QC reports, it is highly recommended a visual inspection of diagnostic plots, as well as the fine-mapping results plots, which can help to detect artifacts.

#### Description of scripts

**`04A_Check_susier_lambda.R`**

-   Read lambda values for all LD blocks and traits and save in a report file.

-   Requires outputs from `03B` (`--diagnostics` step).

-   Generates the commands to launch 04B.

**`04B_Check_bad_finemapping_MN5.R`**

-   Runs independently for each trait and fine-mapping test.

-   Uses lambda report generated by `04A` and performs QC on epigenetic fine-mapping results (outputs from `03C`).

**`04C_Check_finemapping_merge_reports.R`**

-   Merge all tests and all traits fine-mapping QC reports into a single summary.

### 3.5 Generate final fine-mapping results

Final fine-mapping results are generated by integrating QC outputs with the fine-mapping outputs, in order to select the optimal PIP results for each LD block.

**Outputs**

Specifically, we generate two types of ouput files:

-   `summarise_finemapping/*.tsv`: Final fine-mapping results based on SuSiE for LD blocks passing QC.
    It does not include ABF results.
    Used to evaluate fine-mapping performance.

-   `mCOGS_input/*_gwas_4mCOGS.tsv`: full genome-wide summary statistics fincluding all LD blocks (SuSiE + ABF).
    Formatted for input into the mCOGS pipeline.

**Decision rules for final PIP selection *PENDING***

1.  For those were epigenetic fine-mapping was not possible (due to high lambda values), the results from ABF-based results fine mapping was used (coloc function) for COGS. 

    [script: /gpfs/apps/MN5/GPP/R/4.3.2-cairo/INTEL/bin/Rscript /gpfs/scratch/bsc08/bsc008246/Programs/mCOGS/my_scripts/clean/03_Compute_PiP_aBF.R]

2.  Those LD blocks that were flaged in the `use_cv1` column (as fine-mapping with CV5 was not optimal), we take the PIP values from CV1.

3.  For those LD blocks in which using L=5 only 1 CS or 0 CS were found but CV1 found 1 CS, the PIP values from L=1 were used.
    The use of L=1 over L=5 is justified because susieR does not use LD information when assuming a single causal variants and it is a safer option considering that our reference panel is not the best match.

    For those LD blocks in which at least two CS were found, the corresponding PIP for all the SNPs were saved (alpha values from “susie” object) for COGS analysis.

4.  For those LD blocks were fine-mapping with susieR was not performed because they did not contain any signal SNP, ABF results were used (COGS needs genome-wide fine-mapping results).

**Scripts *PENDING***

**05A_Compute_PiP_aBF.R**

Compute ABF

**05B_Convert_SusieR_ouput_for_mCOGS.R**

Select proper results for each LD block.
It also creates

*HELPER SCRIPT PENDING:* `05.0_Convert_SusieR_ouput_generate_commands.R`

### 3.6 Benchmarking annotation combination

`06A_Benchmarking_epitests_generate_tables.R` and `06B_Benchmarking_epitests_summary.R`

These generate a report and plots summarizing, for each functional annotation combination tested:

-   Number of LD blocks with credible sets

-   Number of credible sets

-   Average size of credible sets

-   Number of high-confidence variants (PIP \> 0.5)

Commands generated with the helper script: `06.0_Benchmarking_epitests_generate_commands.R`

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
