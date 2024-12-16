# Date: 16/11/2024 
# query TCGA LUAD data
# Load library
library(tidyverse)
library(TCGAbiolinks)
library(maftools)
library(pheatmap)
library(SummarizedExperiment)
library(DESeq2)

# get a list of projects
gdcprojects <- getGDCprojects()
getProjectSummary('TCGA-LUAD')

# Download gene expression data
# build a query to get gene expression data for entire cohort
query_LUAD_all = GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification",
  sample.type = c("Primary Tumor", "Solid Tissue Normal"),
  access = "open")

output_LUAD_all <- getResults(query_LUAD_all)
# download data - GDCdownload
GDCdownload(query_LUAD_all)

# prepare data
tcga_LUAD_data <- GDCprepare(query_LUAD_all, summarizedExperiment = TRUE)
LUAD_matrix <- assay(tcga_LUAD_data, 'unstranded') # raw data

# Install TCGA clinical data
gene_metadata <- as.data.frame(rowData(tcga_LUAD_data))
coldata <- as.data.frame(colData(tcga_LUAD_data))
colnames(coldata)
# Create code data, Save code: NT and TP
LUAD_code <- coldata |>
  select(barcode, patient, shortLetterCode, days_to_last_follow_up, days_to_death, ajcc_pathologic_stage) 

# mutate LUAD stage variable, stage I-II : early stage, III-IV: late stage
LUAD_code <- LUAD_code |>
  mutate(stage = ifelse(ajcc_pathologic_stage == "Stage I" | ajcc_pathologic_stage == "Stage IA" | ajcc_pathologic_stage == "Stage IB" | ajcc_pathologic_stage == "Stage II" |
                          ajcc_pathologic_stage == "Stage IIA" | ajcc_pathologic_stage == "Stage IIB", "early_stage", "late_stage"))

# metadata survival
LUAD_metadata_survival <- coldata |>
  select(barcode, vital_status, days_to_death, days_to_last_follow_up)

# Check if barcodes in PRAD_code is the same order as in gene expression data (PRAD_matrix)
all(rownames(LUAD_code) == colnames(LUAD_matrix))

# Change shortLetterCode to factor type
LUAD_code$shortLetterCode <- factor(LUAD_code$shortLetterCode)
# Change stage to factor type
LUAD_code$stage <- factor(LUAD_code$stage)

# Create a DGEList object
dds_LUAD <- DESeqDataSetFromMatrix(countData = LUAD_matrix,
                                         colData = LUAD_code,
                                         design = ~ shortLetterCode)
# Create DEG objects between early stage and late stage


# Filter LUAD code data have NA stage
LUAD_code_stage <- LUAD_code |>
  filter(!is.na(stage))

# Filter patients have stage data
LUAD_matrix_stage <- as.data.frame(LUAD_matrix) |>
  select(any_of(LUAD_code_stage$barcode))
df_A <- df_A %>% select(any_of(df_B$x))
all(rownames(LUAD_code_stage) == colnames(LUAD_matrix_stage))
# Filter LUAD matrix to keep only patients that have stage data 
# Create a DGEList object for stage data
dds_LUAD_stage <- DESeqDataSetFromMatrix(countData = LUAD_matrix_stage,
                                   colData = LUAD_code_stage,
                                   design = ~ stage)


# Removing genes with sum total of 10 reads across all samples
keep <- rowSums(counts(dds_LUAD_stage)) >= 10
dds_LUAD_stage <- dds_LUAD_stage[keep,]

# Estimate size factors for normalization
dds_LUAD_stage <- estimateSizeFactors(dds_LUAD_stage)

# normalized counts
normalized_counts_stage <- counts(dds_LUAD_stage, normalized=TRUE)

# Save normalized counts as csv file
normalized_counts_stage <- as.data.frame(normalized_counts_stage) |>
  rownames_to_column(var = "gene_id") |>
  filter(gene_id %in% resLFC_LUAD_filt_stage$gene_id)

# Save normalized counts for stage genes 
normalized_counts_stage <- as.data.frame(no)
# vst: gene count normalization
# calculates a variance stabilizing transformation (VST) from the fitted dispersion-mean relation(s) and then transforms the count data (normalized by division by the size factors or normalization factors
vsd <- vst(dds_LUAD, blind=FALSE) 
LUAD_matrix_vst <- assay(vsd)
LUAD_matrix_vst <- as.data.frame(LUAD_matrix_vst)

# filter gene expression matrix to remain only protein_coding gene
PRAD_matrix_vst_match <- PRAD_matrix_vst_match |>
  rownames_to_column(var = "gene_id")
gene_metadata_protein_coding <- gene_metadata |>
  filter(gene_type == "protein_coding")

# filter gene expression data 
PRAD_matrix_vst_match <- PRAD_matrix_vst_match |>
  filter(gene_id %in% gene_metadata_protein_coding$gene_id)
write.csv(PRAD_matrix_vst_match, '/data/Data/Giang/PRAD_matrix_vst.csv')
# Run DESeq to get DEGs 
dds_LUAD_stage <- DESeq(dds_LUAD_stage)
resLFC_LUAD_stage <- results(dds_LUAD_stage, name = 'stage_late_stage_vs_early_stage')
# 
resLFC_LUAD_stage <- as.data.frame(resLFC_LUAD_stage)
write.csv(resLFC_PRAD_df, '/data/Data/Giang/resLFC_PRAD_df.csv')
# filter only protein coding gene
resLFC_LUAD_stage <- resLFC_LUAD_stage |>
  rownames_to_column(var = "gene_id")
resLFC_LUAD_stage 
# DESeq filter
resLFC_LUAD_filt_stage <- resLFC_LUAD_stage |>
  filter(abs(log2FoldChange) >= 2 & padj <= 0.01)

# Write LUAD_code, normalized LUAD and resLFC in csv
write_csv(LUAD_code, "LUAD_code.csv")
write_csv(normalized_counts_stage, "normalized_counts_stage.csv")
write_csv(resLFC_LUAD_filt_stage, "resLFC_LUAD_filt_stage.csv")

# Write non-filter resLFC bw normal and tumor 
write_csv(resLFC_LUAD, "resLFC_LUAD.csv")
getwd()

# 

