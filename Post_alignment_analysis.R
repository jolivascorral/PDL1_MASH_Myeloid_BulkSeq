###Re-analysis of JOC047 and JOC064 bulkRNAseq experiments 08272026
###Changes in featureCounts to correct/improve for gene annotation and quantification
###Starting post-alignment with STAR and gene counts with featureCounts
###Updated files from featureCounts include combined_counts.txt and combined_counts_summary.txt.

###Pipeline: 
#filter out low-expression/nonprotein-coding/mitochondrial(MT)
#DESeq2 analysis ~experiment + condition (to account for batch variation between exps)
#gene-annotation/ENTREZID
#PCA plots and leading edge genes for PC1 and PC2
#TOP DEGs --> volcano --> heatmap both exps force cluster by genotype
#GSEA enrichment plots --> MsigDB gene lists
#GSVA pathway analysis --> heatmaps/pathway
#Tables for lists of genes included in each gene-set list for GSVA

################################################################################
################################################################################
## UPDATED BULK RNA-SEQ PIPELINE ##
##
####Rebuild analyses after corrected featureCounts####
##
## Analyses:
##   1. JOC047
##   2. JOC064
##   3. Combined JOC047 + JOC064
##
## Filtering:
##   - protein-coding genes only
##   - remove mitochondrial chromosome
##   - remove mt-/MT- gene symbols
##   - expression >= 10 counts in >= 3 samples
##
## DESeq2:
##   JOC047   ~ condition
##   JOC064   ~ condition
##   Combined ~ experiment + condition
##
## All contrasts:
##   CyMt vs WT
################################################################################
################################################################################


################################################################################
##1. LOAD PACKAGES
################################################################################

library(DESeq2)
library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(rtracklayer)
library(SummarizedExperiment)
library(ggplot2)


################################################################################
## 2. DEFINE REFERENCE GTF ####
################################################################################


gtf_file <- file.choose()
gtf_file
################################################################################
## 3. IMPORT JOC047 FEATURECOUNTS FILE
##
## Select the final JOC047 combined_counts.txt

counts_JOC047_raw <- readr::read_tsv(
  file.choose(),
  comment = "#",
  show_col_types = FALSE
)

## 4. IMPORT CORRECTED JOC064 FEATURECOUNTS FILE

counts_JOC064_raw <- readr::read_tsv(
  file.choose(),
  comment = "#",
  show_col_types = FALSE
)

dim(counts_JOC047_raw)
dim(counts_JOC064_raw)

colnames(counts_JOC047_raw)
colnames(counts_JOC064_raw)

################################################################################
## CLEAN SAMPLE NAMES
################################################################################

clean_sample_names <- function(x) {
  
  # Remove leading numbering from featureCounts
  x <- sub("^\\d+_", "", x)
  
  # Keep only WT# or CyMt#
  x <- sub("^(WT[0-9]+).*", "\\1", x)
  x <- sub("^(CyMt[0-9]+).*", "\\1", x)
  
  return(x)
}

# JOC047
colnames(counts_JOC047_raw) <- c(
  colnames(counts_JOC047_raw)[1:6],
  paste0(
    "JOC047_",
    clean_sample_names(colnames(counts_JOC047_raw)[7:ncol(counts_JOC047_raw)])
  )
)

# JOC064
colnames(counts_JOC064_raw) <- c(
  colnames(counts_JOC064_raw)[1:6],
  paste0(
    "JOC064_",
    clean_sample_names(colnames(counts_JOC064_raw)[7:ncol(counts_JOC064_raw)])
  )
)

################################################################################
## 5. CLEAN FEATURECOUNTS SAMPLE COLUMN NAMES
################################################################################

clean_featureCounts_names <- function(count_df) {
  
  annotation_columns <- c(
    "Geneid",
    "Chr",
    "Start",
    "End",
    "Strand",
    "Length"
  )
  
  sample_columns <- setdiff(
    colnames(count_df),
    annotation_columns
  )
  
  clean_names <- basename(sample_columns)
  
  clean_names <- sub(
    "_Aligned\\.sortedByCoord\\.out\\.bam$",
    "",
    clean_names
  )
  
  clean_names <- sub(
    "\\.bam$",
    "",
    clean_names
  )
  
  colnames(count_df)[
    match(sample_columns, colnames(count_df))
  ] <- clean_names
  
  return(count_df)
}


counts_JOC047_raw <- clean_featureCounts_names(
  counts_JOC047_raw
)

counts_JOC064_raw <- clean_featureCounts_names(
  counts_JOC064_raw
)


dim(counts_JOC047_raw)
dim(counts_JOC064_raw)

colnames(counts_JOC047_raw)
colnames(counts_JOC064_raw)

head(counts_JOC047_raw)
head(counts_JOC064_raw)



################################################################################
## 6. IMPORT GTF ####
################################################################################

gtf <- rtracklayer::import(
  gtf_file
)

gtf_df <- as.data.frame(
  gtf
)


################################################################################
## 3. BUILD GENE ANNOTATION ####
################################################################################

gene_anno <- gtf_df %>%
  
  dplyr::filter(
    type == "gene"
  ) %>%
  
  dplyr::transmute(
    
    Geneid = as.character(gene_id),
    
    gene_name = as.character(gene_name),
    
    gene_biotype = as.character(gene_biotype),
    
    chromosome = as.character(seqnames)
    
  ) %>%
  
  dplyr::distinct(
    Geneid,
    .keep_all = TRUE
  )

head(gene_anno)

table(
  gene_anno$gene_biotype,
  useNA = "ifany"
)

unique(gene_anno$chromosome)

################################################################################
## PROTEIN-CODING / NON-MITOCHONDRIAL ANNOTATION####
################################################################################

gene_anno_filtered <- gene_anno %>%
  
  dplyr::filter(
    
    gene_biotype == "protein_coding",
    
    chromosome != "MT",
    
    !is.na(gene_name),
    
    gene_name != "",
    
    !grepl(
      "^mt-",
      gene_name,
      ignore.case = TRUE
    )
    
  )


cat(
  "All annotated genes:",
  nrow(gene_anno),
  "\n"
)

cat(
  "Protein-coding, non-MT genes:",
  nrow(gene_anno_filtered),
  "\n"
)

###############################################################################
## KEEP ONLY GENE ID + COUNTS
################################################################################

counts_JOC047 <- counts_JOC047_raw %>%
  dplyr::select(
    -dplyr::any_of(
      c(
        "Chr",
        "Start",
        "End",
        "Strand",
        "Length"
      )
    )
  )


counts_JOC064 <- counts_JOC064_raw %>%
  dplyr::select(
    -dplyr::any_of(
      c(
        "Chr",
        "Start",
        "End",
        "Strand",
        "Length"
      )
    )
  )

################################################################################
## REORDER SAMPLE COLUMNS
## WT first, then CyMt
################################################################################

counts_JOC047 <- counts_JOC047 %>%
  dplyr::select(
    Geneid,
    JOC047_WT1,
    JOC047_WT2,
    JOC047_WT3,
    JOC047_WT4,
    JOC047_WT5,
    JOC047_CyMt1,
    JOC047_CyMt2,
    JOC047_CyMt3,
    JOC047_CyMt4,
    JOC047_CyMt5
  )

counts_JOC064 <- counts_JOC064 %>%
  dplyr::select(
    Geneid,
    JOC064_WT1,
    JOC064_WT2,
    JOC064_WT3,
    JOC064_CyMt1,
    JOC064_CyMt2,
    JOC064_CyMt3,
    JOC064_CyMt4
  )

################################################################################
## JOIN FILTERED GENE ANNOTATION TO COUNTS
##
## This restricts the count tables to:
##   - protein-coding genes
##   - non-MT chromosome
##   - valid gene symbols
################################################################################

counts_JOC047_annotated <- counts_JOC047 %>%
  dplyr::inner_join(
    gene_anno_filtered %>%
      dplyr::select(
        Geneid,
        gene_name
      ),
    by = "Geneid"
  ) %>%
  dplyr::relocate(
    gene_name,
    .after = Geneid
  )


counts_JOC064_annotated <- counts_JOC064 %>%
  dplyr::inner_join(
    gene_anno_filtered %>%
      dplyr::select(
        Geneid,
        gene_name
      ),
    by = "Geneid"
  ) %>%
  dplyr::relocate(
    gene_name,
    .after = Geneid
  )

################################################################################
## CREATE COMBINED COUNT TABLE
################################################################################

counts_combined_annotated <- counts_JOC047_annotated %>%
  dplyr::inner_join(
    counts_JOC064_annotated,
    by = c(
      "Geneid",
      "gene_name"
    )
  )

################################################################################
## LOW-EXPRESSION FILTER
##
## Keep genes with >= 10 counts in >= 3 samples
################################################################################

filter_low_expression <- function(
    annotated_counts,
    min_count = 10,
    min_samples = 3
) {
  
  count_matrix <- annotated_counts %>%
    dplyr::select(
      -Geneid,
      -gene_name
    ) %>%
    as.matrix()
  
  keep <- rowSums(
    count_matrix >= min_count
  ) >= min_samples
  
  filtered <- annotated_counts[
    keep,
    ,
    drop = FALSE
  ]
  
  message(
    "Genes before low-expression filtering: ",
    nrow(annotated_counts)
  )
  
  message(
    "Genes after low-expression filtering: ",
    nrow(filtered)
  )
  
  return(filtered)
}


counts_JOC047_filtered <- filter_low_expression(
  counts_JOC047_annotated
)

counts_JOC064_filtered <- filter_low_expression(
  counts_JOC064_annotated
)


counts_combined_filtered <- filter_low_expression(
  counts_combined_annotated
)

################################################################################
## VERIFY MITOCHONDRIAL GENES ARE ABSENT
################################################################################

sum(
  grepl(
    "^mt-",
    counts_JOC047_filtered$gene_name,
    ignore.case = TRUE
  )
)

sum(
  grepl(
    "^mt-",
    counts_JOC064_filtered$gene_name,
    ignore.case = TRUE
  )
)

sum(
  grepl(
    "^mt-",
    counts_combined_filtered$gene_name,
    ignore.case = TRUE
  )
)

#Verify everything retained came from filtered-protein coding genes
all(
  counts_JOC047_filtered$Geneid %in%
    gene_anno_filtered$Geneid
)

all(
  counts_JOC064_filtered$Geneid %in%
    gene_anno_filtered$Geneid
)

all(
  counts_combined_filtered$Geneid %in%
    gene_anno_filtered$Geneid
)

################################################################################
## CREATE DESeq2 INTEGER COUNT MATRICES ####
################################################################################

make_DESeq_matrix <- function(df) {
  
  mat <- df %>%
    dplyr::select(
      -Geneid,
      -gene_name
    ) %>%
    as.matrix()
  
  storage.mode(mat) <- "integer"
  
  rownames(mat) <- df$Geneid
  
  return(mat)
}


count_matrix_JOC047 <- make_DESeq_matrix(
  counts_JOC047_filtered
)

count_matrix_JOC064 <- make_DESeq_matrix(
  counts_JOC064_filtered
)

count_matrix_combined <- make_DESeq_matrix(
  counts_combined_filtered
)

################################################################################
## CREATE METADATA
################################################################################

make_metadata <- function(sample_names) {
  
  metadata <- data.frame(
    
    Sample = sample_names,
    
    experiment = dplyr::case_when(
      grepl("^JOC047_", sample_names) ~ "JOC047",
      grepl("^JOC064_", sample_names) ~ "JOC064",
      TRUE ~ NA_character_
    ),
    
    condition = dplyr::case_when(
      grepl("_WT[0-9]+$", sample_names) ~ "WT",
      grepl("_CyMt[0-9]+$", sample_names) ~ "CyMt",
      TRUE ~ NA_character_
    ),
    
    stringsAsFactors = FALSE
  )
  
  rownames(metadata) <- metadata$Sample
  
  metadata$Sample <- NULL
  
  return(metadata)
}


metadata_JOC047 <- make_metadata(
  colnames(count_matrix_JOC047)
)

metadata_JOC064 <- make_metadata(
  colnames(count_matrix_JOC064)
)

metadata_combined <- make_metadata(
  colnames(count_matrix_combined)
)

metadata_JOC047
metadata_JOC064
metadata_combined

################################################################################
## SET FACTOR LEVELS to normalize data to WT 
################################################################################

metadata_JOC047$condition <- factor(
  metadata_JOC047$condition,
  levels = c("WT", "CyMt")
)

metadata_JOC064$condition <- factor(
  metadata_JOC064$condition,
  levels = c("WT", "CyMt")
)

metadata_combined$condition <- factor(
  metadata_combined$condition,
  levels = c("WT", "CyMt")
)

metadata_combined$experiment <- factor(
  metadata_combined$experiment,
  levels = c("JOC047", "JOC064")
)

#Verify sample alignment before DESeq2
stopifnot(
  identical(
    rownames(metadata_JOC047),
    colnames(count_matrix_JOC047)
  )
)

stopifnot(
  identical(
    rownames(metadata_JOC064),
    colnames(count_matrix_JOC064)
  )
)

stopifnot(
  identical(
    rownames(metadata_combined),
    colnames(count_matrix_combined)
  )
)

#### CREATE the three DESeq2 objects####
################################################################################
## JOC047
################################################################################

dds_JOC047 <- DESeq2::DESeqDataSetFromMatrix(
  countData = count_matrix_JOC047,
  colData = metadata_JOC047,
  design = ~ condition
)


################################################################################
## JOC064
################################################################################

dds_JOC064 <- DESeq2::DESeqDataSetFromMatrix(
  countData = count_matrix_JOC064,
  colData = metadata_JOC064,
  design = ~ condition
)


################################################################################
## COMBINED
##
## Experiment explicitly accounts for JOC047 vs JOC064 batch variation
## with experimental design ~ experiment + condition
################################################################################

dds_combined <- DESeq2::DESeqDataSetFromMatrix(
  countData = count_matrix_combined,
  colData = metadata_combined,
  design = ~ experiment + condition
)

################################################################################
## RUN DESEQ2
################################################################################

dds_JOC047 <- DESeq2::DESeq(
  dds_JOC047
)

dds_JOC064 <- DESeq2::DESeq(
  dds_JOC064
)

dds_combined <- DESeq2::DESeq(
  dds_combined
)

resultsNames(dds_JOC047)

resultsNames(dds_JOC064)

resultsNames(dds_combined)

################################################################################
## EXTRACT CYMT VS WT RESULTS
################################################################################

res_JOC047 <- DESeq2::results(
  dds_JOC047,
  contrast = c("condition", "CyMt", "WT"),
  alpha = 0.05
)

res_JOC064 <- DESeq2::results(
  dds_JOC064,
  contrast = c("condition", "CyMt", "WT"),
  alpha = 0.05
)

res_combined <- DESeq2::results(
  dds_combined,
  contrast = c("condition", "CyMt", "WT"),
  alpha = 0.05
)

################################################################################
## ADD GENE SYMBOLS
################################################################################

annotate_DESeq_results <- function(res_object, annotation_table) {
  
  as.data.frame(res_object) %>%
    tibble::rownames_to_column("Geneid") %>%
    dplyr::left_join(
      annotation_table %>%
        dplyr::select(
          Geneid,
          gene_name
        ),
      by = "Geneid"
    ) %>%
    dplyr::relocate(
      gene_name,
      .after = Geneid
    ) %>%
    dplyr::arrange(padj)
}


res_JOC047_df <- annotate_DESeq_results(
  res_JOC047,
  gene_anno_filtered
)

res_JOC064_df <- annotate_DESeq_results(
  res_JOC064,
  gene_anno_filtered
)

res_combined_df <- annotate_DESeq_results(
  res_combined,
  gene_anno_filtered
)
################################################################################
## SIGNIFICANT DEGs
################################################################################

################################################################################
## SIGNIFICANT DEGs
##
## Criteria:
##   padj < 0.05
##   |log2FC| >= 1.5
################################################################################

make_DEG_table <- function(
    res_df,
    padj_cutoff = 0.05,
    log2FC_cutoff = 0
) {
  
  res_df %>%
    dplyr::filter(
      !is.na(padj),
      padj < padj_cutoff,
      abs(log2FoldChange) >= log2FC_cutoff
    ) %>%
    dplyr::mutate(
      Direction = dplyr::case_when(
        log2FoldChange > 0 ~ "Higher in CyMt",
        log2FoldChange < 0 ~ "Higher in WT"
      )
    )
}
##MAKE DEG TABLES####
DEGs_JOC047 <- make_DEG_table(res_JOC047_df)

DEGs_JOC064 <- make_DEG_table(res_JOC064_df)

DEGs_combined <- make_DEG_table(res_combined_df)

cat("JOC047 DEGs:", nrow(DEGs_JOC047), "\n")
cat("JOC064 DEGs:", nrow(DEGs_JOC064), "\n")
cat("Combined DEGs:", nrow(DEGs_combined), "\n")

table(DEGs_JOC047$Direction)
table(DEGs_JOC064$Direction)
table(DEGs_combined$Direction)

################################################################################
## DEG OVERLAP
################################################################################

genes_JOC047 <- DEGs_JOC047$gene_name
genes_JOC064 <- DEGs_JOC064$gene_name
genes_combined <- DEGs_combined$gene_name

length(intersect(
  genes_JOC047,
  genes_JOC064
))

length(intersect(
  genes_JOC047,
  genes_combined
))

length(intersect(
  genes_JOC064,
  genes_combined
))

length(
  Reduce(
    intersect,
    list(
      genes_JOC047,
      genes_JOC064,
      genes_combined
    )
  )
)

#################################
########Plotting data############

################################################################################
################################################################################
## DOWNSTREAM RNA-SEQ ANALYSIS
##
## JOC047
## JOC064
## Combined: JOC047 + JOC064
##
## Starting objects assumed to exist:
##
## dds_JOC047
## dds_JOC064
## dds_combined
##
## res_JOC047_df
## res_JOC064_df
## res_combined_df
##
## DEGs_JOC047
## DEGs_JOC064
## DEGs_combined
##
## gene_anno_filtered
################################################################################
################################################################################


################################################################################
## 1. PACKAGES
################################################################################

library(DESeq2)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(limma)
library(clusterProfiler)
library(enrichplot)
library(GSVA)


################################################################################
################################################################################
## PART A — VST
################################################################################
################################################################################


################################################################################
## 2. CREATE VST OBJECTS
################################################################################

vsd_JOC047 <- DESeq2::vst(
  dds_JOC047,
  blind = FALSE
)

vsd_JOC064 <- DESeq2::vst(
  dds_JOC064,
  blind = FALSE
)

vsd_combined <- DESeq2::vst(
  dds_combined,
  blind = FALSE
)


################################################################################
## 3. EXTRACT VST MATRICES
################################################################################

vsd_matrix_JOC047 <- assay(
  vsd_JOC047
)

vsd_matrix_JOC064 <- assay(
  vsd_JOC064
)

vsd_matrix_combined <- assay(
  vsd_combined
)

################################################################################
## 4. BATCH-CORRECT COMBINED VST MATRIX
##
## IMPORTANT:
## This is for PCA / heatmaps / GSVA visualization.
##
## DO NOT use this matrix as input to DESeq2.
################################################################################

condition_design <- model.matrix(
  ~ condition,
  data = as.data.frame(
    colData(vsd_combined)
  )
)


vsd_matrix_combined_batchCorrected <-
  limma::removeBatchEffect(
    vsd_matrix_combined,
    batch = vsd_combined$experiment,
    design = condition_design
  )

################################################################################
################################################################################
## PART B — PCA ####
################################################################################
################################################################################


################################################################################
## 5. PCA FUNCTION
################################################################################

make_PCA_plot <- function(
    expression_matrix,
    metadata,
    plot_title
) {
  
  pca <- prcomp(
    t(expression_matrix),
    scale. = FALSE
  )
  
  percentVar <-
    pca$sdev^2 /
    sum(pca$sdev^2) *
    100
  
  
  pca_df <- data.frame(
    Sample = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    metadata[
      rownames(pca$x),
      ,
      drop = FALSE
    ]
  )
  
  
  pca_df$condition <- factor(
    pca_df$condition,
    levels = c(
      "WT",
      "CyMt"
    )
  )
  
  
  p <- ggplot(
    pca_df,
    aes(
      x = PC1,
      y = PC2,
      color = condition,
      shape = condition
    )
  ) +
    
    geom_point(
      size = 4
    ) +
    
    scale_color_manual(
      values = c(
        "WT" = "blue",
        "CyMt" = "red"
      )
    ) +
    
    scale_shape_manual(
      values = c(
        "WT" = 16,
        "CyMt" = 15
      )
    ) +
    
    labs(
      title = plot_title,
      x = paste0(
        "PC1: ",
        round(percentVar[1], 1),
        "% variance"
      ),
      y = paste0(
        "PC2: ",
        round(percentVar[2], 1),
        "% variance"
      ),
      color = "Condition",
      shape = "Condition"
    ) +
    
    theme_bw() +
    
    theme(
      plot.title = element_text(
        hjust = 0.5
      )
    )
  
  
  return(
    list(
      plot = p,
      PCA = pca,
      data = pca_df,
      percentVar = percentVar
    )
  )
}

################################################################################
## 6. PCA — JOC047
################################################################################

PCA_JOC047 <- make_PCA_plot(
  expression_matrix =
    vsd_matrix_JOC047,
  
  metadata =
    metadata_JOC047,
  
  plot_title =
    "JOC047 — PCA"
)


################################################################################
## 7. PCA — JOC064
################################################################################

PCA_JOC064 <- make_PCA_plot(
  expression_matrix =
    vsd_matrix_JOC064,
  
  metadata =
    metadata_JOC064,
  
  plot_title =
    "JOC064 — PCA"
)


################################################################################
## 8. PCA — COMBINED BATCH CORRECTED
################################################################################

PCA_combined_batchCorrected <- make_PCA_plot(
  expression_matrix =
    vsd_matrix_combined_batchCorrected,
  
  metadata =
    metadata_combined,
  
  plot_title =
    "Combined — Batch-Corrected PCA"
)

PCA_JOC047$plot
PCA_JOC064$plot
PCA_combined_batchCorrected$plot

dir.create(
  "PCA_plots",
  showWarnings = FALSE
)

ggsave(
  "PCA_plots/JOC047_PCA.png",
  PCA_JOC047$plot,
  width = 7,
  height = 6,
  dpi = 300
)

ggsave(
  "PCA_plots/JOC064_PCA.png",
  PCA_JOC064$plot,
  width = 7,
  height = 6,
  dpi = 300
)

ggsave(
  "PCA_plots/Combined_PCA_batchCorrected.png",
  PCA_combined_batchCorrected$plot,
  width = 7,
  height = 6,
  dpi = 300
)

################################################################################
################################################################################
## PART C — VOLCANO PLOTS ####
################################################################################
################################################################################


################################################################################
## 9. VOLCANO FUNCTION
################################################################################

make_volcano <- function(
    res_df,
    plot_title,
    padj_cutoff = 0.05,
    label_number = 15
) {
  
  volcano_df <- res_df %>%
    dplyr::mutate(
      significance = dplyr::case_when(
        !is.na(padj) &
          padj < padj_cutoff &
          log2FoldChange > 0.75 ~ "Higher in CyMt",
        
        !is.na(padj) &
          padj < padj_cutoff &
          log2FoldChange < 0.75 ~ "Higher in WT",
        
        TRUE ~ "Not significant"
      )
    )
  
  
  top_genes <- volcano_df %>%
    dplyr::filter(
      !is.na(padj),
      padj < padj_cutoff,
      !is.na(gene_name)
    ) %>%
    dplyr::arrange(padj) %>%
    dplyr::slice_head(
      n = label_number
    )
  
  
  p <- ggplot(
    volcano_df,
    aes(
      x = log2FoldChange,
      y = -log10(padj),
      color = significance
    )
  ) +
    
    geom_point(
      alpha = 0.65,
      size = 1.5
    ) +
    
    # adjusted p-value threshold
    geom_hline(
      yintercept = -log10(padj_cutoff),
      linetype = "dashed"
    ) +
    
    # log2 fold-change reference lines
    geom_vline(
      xintercept = c(-0.75, 0.75),
      linetype = "dashed"
    ) +
    
    scale_color_manual(
      values = c(
        "Higher in WT" = "blue",
        "Not significant" = "grey",
        "Higher in CyMt" = "red"
      )
    ) +
    
    ggrepel::geom_text_repel(
      data = top_genes,
      aes(
        label = gene_name
      ),
      size = 3,
      max.overlaps = Inf
    ) +
    
    labs(
      title = plot_title,
      x = "log2 fold change (CyMt vs WT)",
      y = "-log10 adjusted p-value",
      color = ""
    ) +
    
    theme_bw() +
    
    theme(
      plot.title = element_text(
        hjust = 0.5
      )
    )
  
  
  return(p)
}

Volcano_JOC047 <- make_volcano(
  res_JOC047_df,
  "JOC047 — CyMt vs WT"
)

Volcano_JOC064 <- make_volcano(
  res_JOC064_df,
  "JOC064 — CyMt vs WT"
)

Volcano_combined <- make_volcano(
  res_combined_df,
  "Combined — CyMt vs WT"
)

Volcano_JOC047
Volcano_JOC064
Volcano_combined

dir.create(
  "Volcano_plots",
  showWarnings = FALSE
)

ggsave(
  "Volcano_plots/JOC047_volcano.png",
  Volcano_JOC047,
  width = 8,
  height = 7,
  dpi = 300
)

ggsave(
  "Volcano_plots/JOC064_volcano.png",
  Volcano_JOC064,
  width = 8,
  height = 7,
  dpi = 300
)

ggsave(
  "Volcano_plots/Combined_volcano.png",
  Volcano_combined,
  width = 8,
  height = 7,
  dpi = 300
)

ggsave(
  filename = "Combined_volcano_plot.pdf",
  plot = Volcano_combined,
  device = "pdf",
  width = 8,
  height = 7,
  units = "in"
)

################################################################################
## COMBINED VOLCANO PLOT
## Point color = log10 mean expression (DESeq2 baseMean)
################################################################################

library(ggplot2)
library(ggrepel)
library(dplyr)

#Add gene symbols back into res object

make_expression_volcano <- function(
    res_df,
    plot_title,
    padj_cutoff = 0.05,
    label_number = 15
) {
  
  volcano_df <- res_df %>%
    dplyr::mutate(
      
      # Log-transform mean normalized expression
      log10_baseMean = log10(baseMean + 1),
      
      # Used only for choosing genes to label
      significant = !is.na(padj) & padj < padj_cutoff
      
    )
  
  
  ##############################################################################
  ## Select top significant genes for labeling
  ##############################################################################
  
  top_genes <- volcano_df %>%
    dplyr::filter(
      significant,
      !is.na(gene_name),
      gene_name != ""
    ) %>%
    dplyr::arrange(padj) %>%
    dplyr::slice_head(
      n = label_number
    )
  
  
  ##############################################################################
  ## Plot
  ##############################################################################
  
  p <- ggplot(
    volcano_df,
    aes(
      x = log2FoldChange,
      y = -log10(padj),
      color = log10_baseMean
    )
  ) +
    
    geom_point(
      alpha = 0.75,
      size = 1.7
    ) +
    
    # padj = 0.05 reference line
    geom_hline(
      yintercept = -log10(padj_cutoff),
      linetype = "dashed"
    ) +
    
    # log2FC = -1 and +1 reference lines
    geom_vline(
      xintercept = c(-1, 1),
      linetype = "dashed"
    ) +
    
    # Continuous expression gradient
    scale_color_viridis_c(
      option = "viridis",
      name = "log10 mean\nexpression"
    ) +
    
    # Label top significant genes
    ggrepel::geom_text_repel(
      data = top_genes,
      aes(
        label = gene_name
      ),
      color = "black",
      size = 3,
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    
    labs(
      title = plot_title,
      x = "log2 fold change (CyMt vs WT)",
      y = "-log10 adjusted p-value"
    ) +
    
    theme_bw() +
    
    theme(
      plot.title = element_text(
        hjust = 0.5
      )
    )
  
  return(p)
}

Volcano_combined_expression <- make_expression_volcano(
  
  res_df = res_combined_df,
  
  plot_title = "Combined — CyMt vs WT\nColored by Mean Expression"
  
)

Volcano_combined_expression

### 09102026 Volcano plot labeling genes that are both padj <0.05 and log2FC <1

library(ggplot2)
library(ggrepel)
library(dplyr)

################################################################################
## VOLCANO PLOT — CyMt vs WT
##
## Red  = padj < 0.05 AND log2FC > +1
## Blue = padj < 0.05 AND log2FC < -1
## Gray = all other genes
##
## Labels = only red/blue genes
## Line   = padj 0.05 only
################################################################################

volcano_df <- res_combined_df %>%
  dplyr::filter(
    !is.na(log2FoldChange),
    !is.na(padj)
  ) %>%
  dplyr::mutate(
    neg_log10_padj = -log10(padj),
    
    Significance = dplyr::case_when(
      padj < 0.05 & log2FoldChange >  1 ~ "Higher in CyMt",
      padj < 0.05 & log2FoldChange < -1 ~ "Higher in WT",
      TRUE                              ~ "Not significant"
    ),
    
    Label_gene = padj < 0.05 &
      abs(log2FoldChange) > 1
  )

if (!"gene_name" %in% colnames(volcano_df)) {
  
  volcano_df <- volcano_df %>%
    dplyr::left_join(
      gene_anno_filtered %>%
        dplyr::select(Geneid, gene_name) %>%
        dplyr::distinct(Geneid, .keep_all = TRUE),
      by = "Geneid"
    )
}

volcano_combined <- ggplot(
  volcano_df,
  aes(
    x = log2FoldChange,
    y = neg_log10_padj,
    color = Significance
  )
) +
  
  geom_point(
    alpha = 0.7,
    size = 2
  ) +
  
  # Only padj = 0.05 cutoff line
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dotted",
    linewidth = 0.7,
    color = "black"
  ) +
  
  # Label genes with padj < 0.05 AND |log2FC| > 1
  ggrepel::geom_text_repel(
    data = volcano_df %>%
      dplyr::filter(Label_gene),
    aes(label = gene_name),
    color = "black",
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.3,
    min.segment.length = 0
  ) +
  
  scale_color_manual(
    values = c(
      "Higher in CyMt" = "red",
      "Higher in WT" = "blue",
      "Not significant" = "grey70"
    ),
    breaks = c(
      "Higher in CyMt",
      "Higher in WT",
      "Not significant"
    )
  ) +
  
  labs(
    title = "CyMt vs WT",
    x = expression(log[2]~Fold~Change),
    y = expression(-log[10]~adjusted~italic(P)),
    color = NULL
  ) +
  
  theme_classic(base_size = 14)

volcano_combined

################################################################################
################################################################################
## PART D — DEG HEATMAPS ####
################################################################################
################################################################################


################################################################################
## 10. MAP ENSEMBL ROWS TO GENE SYMBOLS
################################################################################

convert_matrix_to_symbols <- function(
    expression_matrix,
    annotation_table
) {
  
  mapping <- annotation_table %>%
    
    dplyr::select(
      Geneid,
      gene_name
    ) %>%
    
    dplyr::distinct(
      Geneid,
      .keep_all = TRUE
    )
  
  
  symbols <- mapping$gene_name[
    match(
      rownames(expression_matrix),
      mapping$Geneid
    )
  ]
  
  
  keep <-
    !is.na(symbols) &
    symbols != ""
  
  
  expression_matrix <-
    expression_matrix[
      keep,
      ,
      drop = FALSE
    ]
  
  symbols <- symbols[keep]
  
  
  rownames(expression_matrix) <-
    symbols
  
  
  ##########################################################################
  ## Collapse duplicate symbols by mean
  ##########################################################################
  
  expression_matrix <-
    rowsum(
      expression_matrix,
      group = rownames(expression_matrix)
    )
  
  counts_per_symbol <-
    table(symbols)
  
  expression_matrix <-
    expression_matrix /
    as.numeric(
      counts_per_symbol[
        rownames(expression_matrix)
      ]
    )
  
  
  return(
    expression_matrix
  )
}

VST_symbols_JOC047 <-
  convert_matrix_to_symbols(
    vsd_matrix_JOC047,
    gene_anno_filtered
  )

VST_symbols_JOC064 <-
  convert_matrix_to_symbols(
    vsd_matrix_JOC064,
    gene_anno_filtered
  )

VST_symbols_combined_batchCorrected <-
  convert_matrix_to_symbols(
    vsd_matrix_combined_batchCorrected,
    gene_anno_filtered
  )

################################################################################
## 11. DEG HEATMAP FUNCTION
################################################################################

make_DEG_heatmap <- function(
    expression_matrix,
    DEG_table,
    metadata,
    heatmap_title,
    output_file
) {
  
  genes <- unique(
    DEG_table$gene_name
  )
  
  
  genes <- intersect(
    genes,
    rownames(expression_matrix)
  )
  
  
  mat <- expression_matrix[
    genes,
    ,
    drop = FALSE
  ]
  
  
  ##########################################################################
  ## ALL WT FIRST → ALL CYMT SECOND - forced clustering by condition/group
  ##########################################################################
  
  WT_samples <- rownames(metadata)[
    metadata$condition == "WT"
  ]
  
  CyMt_samples <- rownames(metadata)[
    metadata$condition == "CyMt"
  ]
  
  
  sample_order <- c(
    WT_samples,
    CyMt_samples
  )
  
  
  mat <- mat[
    ,
    sample_order,
    drop = FALSE
  ]
  
  
  ##########################################################################
  ## SAMPLE ANNOTATION
  ##########################################################################
  
  annotation_col <- data.frame(
    
    Condition = metadata[
      sample_order,
      "condition"
    ]
  )
  
  rownames(annotation_col) <-
    sample_order
  
  
  ##########################################################################
  ## HEATMAP
  ##########################################################################
  
  ph <- pheatmap::pheatmap(
    
    mat,
    
    scale = "row",
    
    cluster_rows = TRUE,
    
    cluster_cols = FALSE,
    
    annotation_col =
      annotation_col,
    
    show_colnames = TRUE,
    
    show_rownames =
      nrow(mat) <= 100,
    
    border_color = NA,
    
    main = heatmap_title,
    
    fontsize_col = 8,
    
    filename = output_file,
    
    width = 10,
    
    height = 10
    
  )
  
  
  return(ph)
}

dir.create(
  "DEG_heatmaps",
  showWarnings = FALSE
)


Heatmap_JOC047 <- make_DEG_heatmap(
  expression_matrix =
    VST_symbols_JOC047,
  
  DEG_table =
    DEGs_JOC047,
  
  metadata =
    metadata_JOC047,
  
  heatmap_title =
    "JOC047 — Significant DEGs",
  
  output_file =
    "DEG_heatmaps/JOC047_DEGs_heatmap.png"
)


Heatmap_JOC064 <- make_DEG_heatmap(
  expression_matrix =
    VST_symbols_JOC064,
  
  DEG_table =
    DEGs_JOC064,
  
  metadata =
    metadata_JOC064,
  
  heatmap_title =
    "JOC064 — Significant DEGs",
  
  output_file =
    "DEG_heatmaps/JOC064_DEGs_heatmap.png"
)


Heatmap_combined <- make_DEG_heatmap(
  expression_matrix =
    VST_symbols_combined_batchCorrected,
  
  DEG_table =
    DEGs_combined,
  
  metadata =
    metadata_combined,
  
  heatmap_title =
    "Combined — Significant DEGs, Batch Corrected",
  
  output_file =
    "DEG_heatmaps/Combined_DEGs_batchCorrected_heatmap.png"
)

################################################################################
################################################################################
## PART E — GSEA ####
################################################################################
################################################################################


################################################################################
## 12. CREATE GSEA RANKED LIST
################################################################################

make_GSEA_rank <- function(res_df) {
  
  ranked <- res_df %>%
    
    dplyr::filter(
      !is.na(stat),
      !is.na(gene_name),
      gene_name != ""
    ) %>%
    
    dplyr::arrange(
      desc(stat)
    )
  
  
  ##########################################################################
  ## Remove duplicate symbols
  ## Retain gene with largest absolute statistic
  ##########################################################################
  
  ranked <- ranked %>%
    
    dplyr::arrange(
      desc(abs(stat))
    ) %>%
    
    dplyr::distinct(
      gene_name,
      .keep_all = TRUE
    )
  
  
  geneList <- ranked$stat
  
  names(geneList) <-
    ranked$gene_name
  
  
  geneList <- sort(
    geneList,
    decreasing = TRUE
  )
  
  
  return(geneList)
}


GSEA_rank_JOC047 <-
  make_GSEA_rank(
    res_JOC047_df
  )

GSEA_rank_JOC064 <-
  make_GSEA_rank(
    res_JOC064_df
  )

GSEA_rank_combined <-
  make_GSEA_rank(
    res_combined_df
  )

################################################################################
## 13. STANDARDIZE TERM2GENE FORMAT
################################################################################

standardize_TERM2GENE <- function(x) {
  
  if ("term" %in% colnames(x)) {
    
    x <- x %>%
      dplyr::rename(
        Pathway = term
      )
  }
  
  
  if ("gene" %in% colnames(x)) {
    
    x <- x %>%
      dplyr::rename(
        gene_name = gene
      )
  }
  
  
  x %>%
    dplyr::select(
      Pathway,
      gene_name
    ) %>%
    dplyr::distinct()
}

################################################################################
## GSEA GENE-SET SETUP #####
##
## Re-import Hallmark GMT files
## Recreate Hallmark_TERM2GENE
## Recreate Immune_TERM2GENE
################################################################################

library(clusterProfiler)
library(dplyr)
library(readr)

################################################################################
## 1. IMPORT HALLMARK GMT FILES
################################################################################

HALLMARK_IL6_JAK_STAT3 <- clusterProfiler::read.gmt(
  file.choose()
)

HALLMARK_HYPOXIA <- clusterProfiler::read.gmt(
  file.choose()
)

HALLMARK_GLYCOLYSIS <- clusterProfiler::read.gmt(file.choose())

################################################################################
## 2. CREATE HALLMARK TERM2GENE
################################################################################

Hallmark_TERM2GENE <- dplyr::bind_rows(
  HALLMARK_IL6_JAK_STAT3,
  HALLMARK_HYPOXIA,
  HALLMARK_GLYCOLYSIS
) %>%
  dplyr::select(
    term,
    gene
  ) %>%
  dplyr::distinct()

Hallmark_TERM2GENE <- Hallmark_TERM2GENE %>%
  dplyr::rename(
    gene_name = gene
  )

################################################################################
## 3. IMPORT IMMUNE GMT FILES
################################################################################

MYELOID_ACTIVATION <- clusterProfiler::read.gmt(
  file.choose()
)


################################################################################
## 5. STANDARDIZE TO Pathway + gene_name
################################################################################

Hallmark_TERM2GENE <- Hallmark_TERM2GENE %>%
  dplyr::rename(
    Pathway = term
  )

################################################################################
## DEFINE CUSTOM GSEA FUNCTION
################################################################################

run_custom_GSEA <- function(
    geneList,
    TERM2GENE
) {
  
  clusterProfiler::GSEA(
    
    geneList = geneList,
    
    TERM2GENE = TERM2GENE %>%
      dplyr::select(
        Pathway,
        gene_name
      ),
    
    pvalueCutoff = 1,
    
    pAdjustMethod = "BH",
    
    minGSSize = 10,
    
    maxGSSize = 500,
    
    eps = 0,
    
    verbose = FALSE
  )
}

Hallmark_GSEA_JOC047 <- run_custom_GSEA(
  GSEA_rank_JOC047,
  Hallmark_TERM2GENE
)

Hallmark_GSEA_JOC064 <- run_custom_GSEA(
  GSEA_rank_JOC064,
  Hallmark_TERM2GENE
)

Hallmark_GSEA_combined <- run_custom_GSEA(
  GSEA_rank_combined,
  Hallmark_TERM2GENE
)


################################################################################
## LOAD IMMUNE MSigDB METADATA FILES
################################################################################

myeloid_file <- file.choose()
innate_file  <- file.choose()


read_msigdb_metadata_symbols <- function(file) {
  
  lines <- readLines(file)
  
  pathway_line <- grep(
    "^STANDARD_NAME\\t",
    lines,
    value = TRUE
  )
  
  pathway_name <- sub(
    "^STANDARD_NAME\\t",
    "",
    pathway_line
  )
  
  gene_line <- grep(
    "^GENE_SYMBOLS\\t",
    lines,
    value = TRUE
  )
  
  if (length(gene_line) != 1) {
    stop(
      "Could not find exactly one GENE_SYMBOLS row in: ",
      file
    )
  }
  
  gene_string <- sub(
    "^GENE_SYMBOLS\\t",
    "",
    gene_line
  )
  
  genes <- strsplit(
    gene_string,
    ",",
    fixed = TRUE
  )[[1]]
  
  genes <- trimws(genes)
  
  genes <- genes[
    !is.na(genes) &
      genes != ""
  ]
  
  genes <- unique(genes)
  
  data.frame(
    Pathway = pathway_name,
    gene_name = genes,
    stringsAsFactors = FALSE
  )
}

MYELOID_ACTIVATION <- read_msigdb_metadata_symbols(
  myeloid_file
)

INNATE_IMMUNE <- read_msigdb_metadata_symbols(
  innate_file
)

Immune_TERM2GENE <- dplyr::bind_rows(
  MYELOID_ACTIVATION,
  INNATE_IMMUNE
) %>%
  dplyr::distinct()

Immune_GSEA_JOC047 <- run_custom_GSEA(
  GSEA_rank_JOC047,
  Immune_TERM2GENE
)

Immune_GSEA_JOC064 <- run_custom_GSEA(
  GSEA_rank_JOC064,
  Immune_TERM2GENE
)

Immune_GSEA_combined <- run_custom_GSEA(
  GSEA_rank_combined,
  Immune_TERM2GENE
)

################################################################################
################################################################################
## GSEA ENRICHMENT PLOTS FOR EVERY PATHWAY
################################################################################
################################################################################

library(enrichplot)
library(ggplot2)
library(dplyr)


################################################################################
## 1. FUNCTION TO CREATE ONE CLEAN GSEA ENRICHMENT PLOT
################################################################################

make_clean_gsea_plot <- function(
    GSEA_object,
    pathway_id,
    plot_title
) {
  
  ##########################################################################
  ## Get GSEA statistics for this pathway
  ##########################################################################
  
  result_df <- as.data.frame(
    GSEA_object
  )
  
  stats_row <- result_df %>%
    dplyr::filter(
      ID == pathway_id
    )
  
  
  if (nrow(stats_row) != 1) {
    
    warning(
      "Pathway could not be uniquely identified: ",
      pathway_id
    )
    
    return(NULL)
  }
  
  
  ##########################################################################
  ## Format NES
  ##########################################################################
  
  NES_text <- sprintf(
    "%.2f",
    stats_row$NES
  )
  
  
  ##########################################################################
  ## Format adjusted p-value
  ##########################################################################
  
  padj_text <- ifelse(
    
    stats_row$p.adjust < 0.001,
    
    format(
      stats_row$p.adjust,
      scientific = TRUE,
      digits = 2
    ),
    
    sprintf(
      "%.3f",
      stats_row$p.adjust
    )
  )
  
  
  ##########################################################################
  ## Plot title
  ##########################################################################
  
  full_title <- paste0(
    plot_title,
    "\nNES = ",
    NES_text,
    "   |   padj = ",
    padj_text
  )
  
  
  ##########################################################################
  ## Generate standard enrichment plot
  ##########################################################################
  
  p <- enrichplot::gseaplot2(
    
    GSEA_object,
    
    geneSetID =
      pathway_id,
    
    title =
      full_title,
    
    pvalue_table =
      FALSE
    
  )
  
  
  ##########################################################################
  ## Remove gray ranked-list metric curve
  ##
  ## Keep the red-blue gradient.
  ##########################################################################
  
  if (
    length(p) >= 3 &&
    length(p[[3]]$layers) >= 1
  ) {
    
    p[[3]]$layers <-
      p[[3]]$layers[-1]
    
  }
  
  
  ##########################################################################
  ## Remove x/y axes from the gradient panel
  ##########################################################################
  
  if (length(p) >= 3) {
    
    p[[3]] <- p[[3]] +
      
      ggplot2::theme(
        
        axis.title.x =
          ggplot2::element_blank(),
        
        axis.title.y =
          ggplot2::element_blank(),
        
        axis.text.x =
          ggplot2::element_blank(),
        
        axis.text.y =
          ggplot2::element_blank(),
        
        axis.ticks.x =
          ggplot2::element_blank(),
        
        axis.ticks.y =
          ggplot2::element_blank(),
        
        axis.line.x =
          ggplot2::element_blank(),
        
        axis.line.y =
          ggplot2::element_blank(),
        
        plot.margin =
          ggplot2::margin(
            t = 0,
            r = 5,
            b = 0,
            l = 5
          )
      )
  }
  
  
  ##########################################################################
  ## Remove gray grid lines from all panels
  ##########################################################################
  
  for (i in seq_along(p)) {
    
    p[[i]] <- p[[i]] +
      
      ggplot2::theme(
        
        panel.grid.major =
          ggplot2::element_blank(),
        
        panel.grid.minor =
          ggplot2::element_blank()
      )
  }
  
  
  return(p)
}


################################################################################
## 2. FUNCTION TO PLOT ALL PATHWAYS IN ONE GSEA OBJECT
################################################################################

plot_all_GSEA_pathways <- function(
    GSEA_object,
    dataset_name,
    collection_name,
    output_folder
) {
  
  result_df <- as.data.frame(
    GSEA_object
  )
  
  if (nrow(result_df) == 0) {
    
    message(
      "No GSEA pathways available for ",
      dataset_name,
      " — ",
      collection_name
    )
    
    return(
      invisible(NULL)
    )
  }
  
  dir.create(
    output_folder,
    showWarnings = FALSE,
    recursive = TRUE
  )
  
  pathway_ids <- result_df$ID
  
  GSEA_plots <- list()
  
  for (pathway_id in pathway_ids) {
    
    message(
      "Plotting: ",
      dataset_name,
      " — ",
      pathway_id
    )
    
    safe_name <- gsub(
      "[^A-Za-z0-9_]+",
      "_",
      pathway_id
    )
    
    current_plot <- make_clean_gsea_plot(
      
      GSEA_object =
        GSEA_object,
      
      pathway_id =
        pathway_id,
      
      plot_title =
        paste0(
          dataset_name,
          " — ",
          pathway_id
        )
    )
    
    if (is.null(current_plot)) {
      next
    }
    
    GSEA_plots[[pathway_id]] <-
      current_plot
    
    ggplot2::ggsave(
      
      filename =
        file.path(
          output_folder,
          paste0(
            dataset_name,
            "_",
            safe_name,
            "_GSEA.png"
          )
        ),
      
      plot =
        current_plot,
      
      width = 9,
      height = 7,
      dpi = 300
    )
  }
  
  return(
    GSEA_plots
  )
}

################################################################################
## 3. HALLMARK — JOC047
################################################################################

Hallmark_GSEA_plots_JOC047 <-
  plot_all_GSEA_pathways(
    
    GSEA_object =
      Hallmark_GSEA_JOC047,
    
    dataset_name =
      "JOC047",
    
    collection_name =
      "Hallmark",
    
    output_folder =
      "GSEA_enrichment_plots/Hallmark/JOC047"
  )


################################################################################
## 4. HALLMARK — JOC064
################################################################################

Hallmark_GSEA_plots_JOC064 <-
  plot_all_GSEA_pathways(
    
    GSEA_object =
      Hallmark_GSEA_JOC064,
    
    dataset_name =
      "JOC064",
    
    collection_name =
      "Hallmark",
    
    output_folder =
      "GSEA_enrichment_plots/Hallmark/JOC064"
  )


################################################################################
## 5. HALLMARK — COMBINED
################################################################################

Hallmark_GSEA_plots_combined <-
  plot_all_GSEA_pathways(
    
    GSEA_object =
      Hallmark_GSEA_combined,
    
    dataset_name =
      "Combined",
    
    collection_name =
      "Hallmark",
    
    output_folder =
      "GSEA_enrichment_plots/Hallmark/Combined"
  )

################################################################################
## 6. IMMUNE — JOC047
################################################################################

Immune_GSEA_plots_JOC047 <-
  plot_all_GSEA_pathways(
    
    GSEA_object =
      Immune_GSEA_JOC047,
    
    dataset_name =
      "JOC047",
    
    collection_name =
      "Immune",
    
    output_folder =
      "GSEA_enrichment_plots/Immune/JOC047"
  )


################################################################################
## 7. IMMUNE — JOC064
################################################################################

Immune_GSEA_plots_JOC064 <-
  plot_all_GSEA_pathways(
    
    GSEA_object =
      Immune_GSEA_JOC064,
    
    dataset_name =
      "JOC064",
    
    collection_name =
      "Immune",
    
    output_folder =
      "GSEA_enrichment_plots/Immune/JOC064"
  )


################################################################################
## 8. IMMUNE — COMBINED
################################################################################

Immune_GSEA_plots_combined <-
  plot_all_GSEA_pathways(
    
    GSEA_object =
      Immune_GSEA_combined,
    
    dataset_name =
      "Combined",
    
    collection_name =
      "Immune",
    
    output_folder =
      "GSEA_enrichment_plots/Immune/Combined"
  )

################################################################################
## SAVE A LIST OF GSEA PLOTS AS PDF
################################################################################

save_GSEA_plots_pdf <- function(plot_list, outdir) {
  
  dir.create(
    outdir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  for (pathway in names(plot_list)) {
    
    ggsave(
      
      filename = file.path(
        outdir,
        paste0(pathway, ".pdf")
      ),
      
      plot = plot_list[[pathway]],
      
      device = "pdf",
      
      width = 8,
      
      height = 7,
      
      units = "in"
      
    )
    
  }
  
}

save_GSEA_plots_pdf(
  Hallmark_GSEA_plots_combined,
  "GSEA_PDFs/Hallmark/Combined"
)

save_GSEA_plots_pdf(
  Immune_GSEA_plots_combined,
  "GSEA_PDFs/Immune/Combined"
)









################################################################################
################################################################################
## PATHWAY-LEVEL HEATMAPS + GENE TABLES
##
## Works for:
##   Hallmark_TERM2GENE
##   Immune_TERM2GENE
##
## Analyses:
##   JOC047
##   JOC064
##   Combined
##
## Outputs:
##   - one heatmap per pathway per dataset
##   - one TSV per pathway per dataset
##
## Table includes:
##   gene_name
##   baseMean
##   log2FoldChange
##   lfcSE
##   stat
##   pvalue
##   padj
##   Direction
################################################################################
################################################################################


library(dplyr)
library(tidyr)
library(readr)
library(pheatmap)


################################################################################
## 1. HELPER: CREATE PATHWAY GENE TABLE
################################################################################

make_pathway_gene_table <- function(
    pathway_name,
    TERM2GENE,
    res_df
) {
  
  pathway_genes <- TERM2GENE %>%
    dplyr::filter(
      Pathway == pathway_name
    ) %>%
    dplyr::select(
      Pathway,
      gene_name
    ) %>%
    dplyr::distinct()
  
  
  pathway_table <- pathway_genes %>%
    
    dplyr::left_join(
      res_df %>%
        dplyr::select(
          gene_name,
          baseMean,
          log2FoldChange,
          lfcSE,
          stat,
          pvalue,
          padj
        ),
      by = "gene_name"
    ) %>%
    
    dplyr::mutate(
      
      Direction = dplyr::case_when(
        
        is.na(log2FoldChange) ~
          "Not detected",
        
        log2FoldChange > 0 ~
          "Higher in CyMt",
        
        log2FoldChange < 0 ~
          "Higher in WT",
        
        TRUE ~
          "No change"
      ),
      
      Significant = dplyr::case_when(
        
        !is.na(padj) &
          padj < 0.05 ~
          "padj < 0.05",
        
        TRUE ~
          "Not significant"
      )
    ) %>%
    
    dplyr::arrange(
      padj,
      dplyr::desc(
        abs(log2FoldChange)
      )
    )
  
  
  return(
    pathway_table
  )
}

################################################################################
## 2. HELPER: PATHWAY HEATMAP
################################################################################

make_pathway_heatmap <- function(
    pathway_name,
    TERM2GENE,
    expression_matrix,
    metadata,
    dataset_name,
    output_folder
) {
  
  ##########################################################################
  ## Identify genes in pathway
  ##########################################################################
  
  pathway_genes <- TERM2GENE %>%
    
    dplyr::filter(
      Pathway == pathway_name
    ) %>%
    
    dplyr::pull(
      gene_name
    ) %>%
    
    unique()
  
  
  ##########################################################################
  ## Keep pathway genes found in expression matrix
  ##########################################################################
  
  genes_present <- intersect(
    pathway_genes,
    rownames(expression_matrix)
  )
  
  
  if (length(genes_present) == 0) {
    
    message(
      "No pathway genes found for: ",
      pathway_name,
      " — ",
      dataset_name
    )
    
    return(
      invisible(NULL)
    )
  }
  
  
  mat <- expression_matrix[
    genes_present,
    ,
    drop = FALSE
  ]
  
  
  ##########################################################################
  ## Put all WT samples first, then all CyMt samples
  ##########################################################################
  
  WT_samples <- rownames(metadata)[
    metadata$condition == "WT"
  ]
  
  CyMt_samples <- rownames(metadata)[
    metadata$condition == "CyMt"
  ]
  
  
  sample_order <- c(
    WT_samples,
    CyMt_samples
  )
  
  
  mat <- mat[
    ,
    sample_order,
    drop = FALSE
  ]
  
  
  ##########################################################################
  ## Column annotation
  ##########################################################################
  
  annotation_col <- data.frame(
    
    Condition = metadata[
      sample_order,
      "condition"
    ]
    
  )
  
  rownames(annotation_col) <-
    sample_order
  
  
  ##########################################################################
  ## Safe filename
  ##########################################################################
  
  safe_pathway_name <- gsub(
    "[^A-Za-z0-9_]+",
    "_",
    pathway_name
  )
  
  
  dir.create(
    output_folder,
    showWarnings = FALSE,
    recursive = TRUE
  )
  
  
  output_file <- file.path(
    output_folder,
    paste0(
      dataset_name,
      "_",
      safe_pathway_name,
      "_heatmap.png"
    )
  )
  
  
  ##########################################################################
  ## Heatmap
  ##########################################################################
  
  ph <- pheatmap::pheatmap(
    
    mat,
    
    scale = "row",
    
    cluster_rows = TRUE,
    
    cluster_cols = FALSE,
    
    annotation_col =
      annotation_col,
    
    show_colnames = TRUE,
    
    show_rownames = TRUE,
    
    border_color = NA,
    
    main = paste0(
      dataset_name,
      " — ",
      pathway_name
    ),
    
    fontsize_row = 8,
    
    fontsize_col = 8,
    
    filename = output_file,
    
    width = 10,
    
    height = max(
      7,
      0.18 * nrow(mat)
    )
    
  )
  
  
  return(
    ph
  )
}

################################################################################
## 3. RUN ALL PATHWAYS IN A COLLECTION
################################################################################

run_pathway_collection <- function(
    TERM2GENE,
    res_df,
    expression_matrix,
    metadata,
    dataset_name,
    collection_name
) {
  
  pathway_names <- unique(
    TERM2GENE$Pathway
  )
  
  
  heatmap_folder <- file.path(
    "Pathway_heatmaps",
    collection_name,
    dataset_name
  )
  
  
  table_folder <- file.path(
    "Pathway_gene_tables",
    collection_name,
    dataset_name
  )
  
  
  dir.create(
    heatmap_folder,
    showWarnings = FALSE,
    recursive = TRUE
  )
  
  
  dir.create(
    table_folder,
    showWarnings = FALSE,
    recursive = TRUE
  )
  
  
  pathway_tables <- list()
  
  
  for (pathway_name in pathway_names) {
    
    message(
      "Processing ",
      dataset_name,
      " — ",
      pathway_name
    )
    
    
    ########################################################################
    ## Create gene-level statistics table
    ########################################################################
    
    pathway_table <- make_pathway_gene_table(
      
      pathway_name =
        pathway_name,
      
      TERM2GENE =
        TERM2GENE,
      
      res_df =
        res_df
      
    )
    
    
    pathway_tables[[
      pathway_name
    ]] <- pathway_table
    
    
    ########################################################################
    ## Save gene table
    ########################################################################
    
    safe_name <- gsub(
      "[^A-Za-z0-9_]+",
      "_",
      pathway_name
    )
    
    
    readr::write_tsv(
      
      pathway_table,
      
      file.path(
        table_folder,
        paste0(
          dataset_name,
          "_",
          safe_name,
          "_genes.tsv"
        )
      )
    )
    
    
    ########################################################################
    ## Create heatmap
    ########################################################################
    
    make_pathway_heatmap(
      
      pathway_name =
        pathway_name,
      
      TERM2GENE =
        TERM2GENE,
      
      expression_matrix =
        expression_matrix,
      
      metadata =
        metadata,
      
      dataset_name =
        dataset_name,
      
      output_folder =
        heatmap_folder
      
    )
    
  }
  
  
  return(
    pathway_tables
  )
}

################################################################################
## 4. HALLMARK PATHWAYS
################################################################################

Hallmark_tables_JOC047 <- run_pathway_collection(
  
  TERM2GENE =
    Hallmark_TERM2GENE,
  
  res_df =
    res_JOC047_df,
  
  expression_matrix =
    VST_symbols_JOC047,
  
  metadata =
    metadata_JOC047,
  
  dataset_name =
    "JOC047",
  
  collection_name =
    "Hallmark"
)


Hallmark_tables_JOC064 <- run_pathway_collection(
  
  TERM2GENE =
    Hallmark_TERM2GENE,
  
  res_df =
    res_JOC064_df,
  
  expression_matrix =
    VST_symbols_JOC064,
  
  metadata =
    metadata_JOC064,
  
  dataset_name =
    "JOC064",
  
  collection_name =
    "Hallmark"
)


Hallmark_tables_combined <- run_pathway_collection(
  
  TERM2GENE =
    Hallmark_TERM2GENE,
  
  res_df =
    res_combined_df,
  
  expression_matrix =
    VST_symbols_combined_batchCorrected,
  
  metadata =
    metadata_combined,
  
  dataset_name =
    "Combined",
  
  collection_name =
    "Hallmark"
)

################################################################################
## 5. IMMUNE PATHWAYS
################################################################################

Immune_tables_JOC047 <- run_pathway_collection(
  
  TERM2GENE =
    Immune_TERM2GENE,
  
  res_df =
    res_JOC047_df,
  
  expression_matrix =
    VST_symbols_JOC047,
  
  metadata =
    metadata_JOC047,
  
  dataset_name =
    "JOC047",
  
  collection_name =
    "Immune"
)


Immune_tables_JOC064 <- run_pathway_collection(
  
  TERM2GENE =
    Immune_TERM2GENE,
  
  res_df =
    res_JOC064_df,
  
  expression_matrix =
    VST_symbols_JOC064,
  
  metadata =
    metadata_JOC064,
  
  dataset_name =
    "JOC064",
  
  collection_name =
    "Immune"
)


Immune_tables_combined <- run_pathway_collection(
  
  TERM2GENE =
    Immune_TERM2GENE,
  
  res_df =
    res_combined_df,
  
  expression_matrix =
    VST_symbols_combined_batchCorrected,
  
  metadata =
    metadata_combined,
  
  dataset_name =
    "Combined",
  
  collection_name =
    "Immune"
)



################################################################################
## DEG-ONLY PATHWAY HEATMAP — COMBINED ####
##
## Includes only genes that are:
##   1. members of Hallmark or Immune pathways
##   2. significant DEGs (padj < 0.05)
##
## Rows:
##   gene names
##
## Row annotation:
##   pathway / term
##
## Columns:
##   WT samples first, then CyMt
##
## Expression:
##   batch-corrected VST
################################################################################

library(dplyr)
library(pheatmap)


################################################################################
## 1. COMBINE ALL PATHWAY DEFINITIONS
################################################################################

All_PATHWAYS_TERM2GENE <- dplyr::bind_rows(
  
  Hallmark_TERM2GENE %>%
    dplyr::mutate(
      Collection = "Hallmark"
    ),
  
  Immune_TERM2GENE %>%
    dplyr::mutate(
      Collection = "Immune"
    )
  
) %>%
  dplyr::distinct()


################################################################################
## 2. IDENTIFY PATHWAY GENES THAT ARE SIGNIFICANT DEGs
################################################################################

Pathway_DEGs_combined <- All_PATHWAYS_TERM2GENE %>%
  
  dplyr::inner_join(
    
    DEGs_combined %>%
      dplyr::select(
        gene_name,
        log2FoldChange,
        pvalue,
        padj,
        Direction
      ),
    
    by = "gene_name"
    
  ) %>%
  
  dplyr::arrange(
    Pathway,
    padj
  )


################################################################################
## 3. CHECK RESULT
################################################################################

dim(
  Pathway_DEGs_combined
)

head(
  Pathway_DEGs_combined
)

table(
  Pathway_DEGs_combined$Pathway
)

################################################################################
## 4. COLLAPSE MULTIPLE PATHWAY MEMBERSHIPS PER GENE
################################################################################

Pathway_DEG_annotation <- Pathway_DEGs_combined %>%
  
  dplyr::group_by(
    gene_name
  ) %>%
  
  dplyr::summarise(
    
    Term = paste(
      unique(Pathway),
      collapse = " | "
    ),
    
    Collection = paste(
      unique(Collection),
      collapse = " | "
    ),
    
    log2FoldChange = dplyr::first(
      log2FoldChange
    ),
    
    pvalue = dplyr::first(
      pvalue
    ),
    
    padj = dplyr::first(
      padj
    ),
    
    Direction = dplyr::first(
      Direction
    ),
    
    .groups = "drop"
  )

################################################################################
## 5. EXTRACT DEG EXPRESSION MATRIX
################################################################################

pathway_deg_genes <- Pathway_DEG_annotation$gene_name


genes_present <- intersect(
  pathway_deg_genes,
  rownames(
    VST_symbols_combined_batchCorrected
  )
)


Pathway_DEG_matrix_combined <-
  VST_symbols_combined_batchCorrected[
    genes_present,
    ,
    drop = FALSE
  ]

################################################################################
## 6. SAMPLE ORDER
################################################################################

WT_samples <- rownames(
  metadata_combined
)[
  metadata_combined$condition == "WT"
]

CyMt_samples <- rownames(
  metadata_combined
)[
  metadata_combined$condition == "CyMt"
]


sample_order <- c(
  WT_samples,
  CyMt_samples
)


Pathway_DEG_matrix_combined <-
  Pathway_DEG_matrix_combined[
    ,
    sample_order,
    drop = FALSE
  ]

################################################################################
## 7. COLUMN ANNOTATION
################################################################################

annotation_col <- data.frame(
  
  Condition =
    metadata_combined[
      sample_order,
      "condition"
    ],
  
  Experiment =
    metadata_combined[
      sample_order,
      "experiment"
    ]
  
)

rownames(
  annotation_col
) <- sample_order

################################################################################
## 8. ROW ANNOTATION
################################################################################
################################################################################
## REBUILD CLEAN PATHWAY-DEG HEATMAP
##
## IMPORTANT:
## - Gene symbols ONLY as row labels
## - Pathway shown as a short categorical annotation bar
## - WT samples first, CyMt samples second
## - Genes hierarchically clustered
## - Columns NOT clustered, preserving condition grouping
################################################################################


################################################################################
## 1. REBUILD THE MATRIX FRESH
##
## This is necessary because the previous code changed the rownames.
################################################################################

genes_present <- intersect(
  Pathway_DEG_annotation$gene_name,
  rownames(VST_symbols_combined_batchCorrected)
)

Pathway_DEG_matrix_combined <-
  VST_symbols_combined_batchCorrected[
    genes_present,
    ,
    drop = FALSE
  ]


################################################################################
## 2. CREATE SHORT PATHWAY LABELS
################################################################################

Pathway_DEG_annotation_plot <- Pathway_DEG_annotation %>%
  dplyr::filter(
    gene_name %in% genes_present
  ) %>%
  dplyr::mutate(
    
    Term_short = dplyr::case_when(
      
      grepl(
        "IL6_JAK_STAT3",
        Term
      ) ~ "IL6/JAK/STAT3",
      
      grepl(
        "HYPOXIA",
        Term
      ) ~ "Hypoxia",
      
      grepl(
        "GLYCOLYSIS",
        Term
      ) ~ "Glycolysis",
      
      grepl(
        "MYELOID_CELL_ACTIVATION",
        Term
      ) ~ "Myeloid activation",
      
      grepl(
        "INNATE_IMMUNE_RESPONSE",
        Term
      ) ~ "Innate immune",
      
      TRUE ~ "Other"
    )
  )

################################################################################
## 3. MATCH ROW ANNOTATION TO MATRIX ORDER
################################################################################

Pathway_DEG_annotation_plot <-
  Pathway_DEG_annotation_plot[
    match(
      rownames(Pathway_DEG_matrix_combined),
      Pathway_DEG_annotation_plot$gene_name
    ),
    ,
    drop = FALSE
  ]


stopifnot(
  identical(
    rownames(Pathway_DEG_matrix_combined),
    Pathway_DEG_annotation_plot$gene_name
  )
)

################################################################################
## 4. CREATE SIMPLE ROW ANNOTATION
################################################################################

annotation_row <- data.frame(
  
  Pathway =
    Pathway_DEG_annotation_plot$Term_short
  
)

rownames(annotation_row) <-
  Pathway_DEG_annotation_plot$gene_name

################################################################################
## 5. ORDER SAMPLES
## WT first, then CyMt
################################################################################

WT_samples <- rownames(metadata_combined)[
  metadata_combined$condition == "WT"
]

CyMt_samples <- rownames(metadata_combined)[
  metadata_combined$condition == "CyMt"
]

sample_order <- c(
  WT_samples,
  CyMt_samples
)

Pathway_DEG_matrix_combined <-
  Pathway_DEG_matrix_combined[
    ,
    sample_order,
    drop = FALSE
  ]

################################################################################
## 6. COLUMN ANNOTATION
################################################################################

annotation_col <- data.frame(
  
  Condition =
    metadata_combined[
      sample_order,
      "condition"
    ]
  
)

rownames(annotation_col) <-
  sample_order

################################################################################
## 7. DISPLAY CLEAN HEATMAP IN RSTUDIO
################################################################################

Pathway_DEG_heatmap_combined <- pheatmap::pheatmap(
  
  Pathway_DEG_matrix_combined,
  
  scale = "row",
  
  cluster_rows = TRUE,
  
  cluster_cols = FALSE,
  
  annotation_row = annotation_row,
  
  annotation_col = annotation_col,
  
  # Gene names only
  labels_row =
    rownames(Pathway_DEG_matrix_combined),
  
  show_rownames = TRUE,
  
  show_colnames = TRUE,
  
  fontsize_row = 9,
  
  fontsize_col = 9,
  
  angle_col = 45,
  
  cellwidth = 25,
  
  cellheight = 14,
  
  border_color = NA,
  
  treeheight_row = 50,
  
  annotation_names_row = FALSE,
  
  annotation_names_col = TRUE,
  
  main = "Combined — Significant Pathway DEGs"
  
)

################################################################################
## 8. SAVE HIGH-RESOLUTION PDF
################################################################################

pheatmap::pheatmap(
  
  Pathway_DEG_matrix_combined,
  
  scale = "row",
  
  cluster_rows = TRUE,
  
  cluster_cols = FALSE,
  
  annotation_row = annotation_row,
  
  annotation_col = annotation_col,
  
  labels_row =
    rownames(Pathway_DEG_matrix_combined),
  
  show_rownames = TRUE,
  
  show_colnames = TRUE,
  
  fontsize_row = 9,
  
  fontsize_col = 9,
  
  angle_col = 45,
  
  cellwidth = 25,
  
  cellheight = 14,
  
  border_color = NA,
  
  treeheight_row = 50,
  
  annotation_names_row = FALSE,
  
  annotation_names_col = TRUE,
  
  main = "Combined — Significant Pathway DEGs",
  
  filename =
    "Combined_Pathway_DEGs_heatmap_readable.pdf",
  
  width = 14,
  
  height = max(
    10,
    0.25 *
      nrow(Pathway_DEG_matrix_combined)
  )
)


readr::write_tsv(
  
  Pathway_DEG_annotation %>%
    dplyr::arrange(
      Term,
      padj
    ),
  
  "Combined_Pathway_DEGs_gene_term_table.tsv"
)

################################################################################
## PATHWAY SUMMARY STATISTICS
################################################################################

Pathway_summary_combined <-
  Pathway_DEGs_combined %>%
  dplyr::group_by(
    Collection,
    Pathway
  ) %>%
  dplyr::summarise(
    
    n_genes = dplyr::n(),
    
    n_up_CyMt = sum(
      log2FoldChange > 0,
      na.rm = TRUE
    ),
    
    n_up_WT = sum(
      log2FoldChange < 0,
      na.rm = TRUE
    ),
    
    median_log2FC = median(
      log2FoldChange,
      na.rm = TRUE
    ),
    
    mean_log2FC = mean(
      log2FoldChange,
      na.rm = TRUE
    ),
    
    mean_abs_log2FC = mean(
      abs(log2FoldChange),
      na.rm = TRUE
    ),
    
    median_padj = median(
      padj,
      na.rm = TRUE
    ),
    
    min_padj = min(
      padj,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  
  dplyr::mutate(
    
    percent_up_CyMt = round(
      100 * n_up_CyMt / n_genes,
      1
    ),
    
    percent_up_WT = round(
      100 * n_up_WT / n_genes,
      1
    )
  ) %>%
  
  dplyr::arrange(
    Collection,
    median_padj
  )

readr::write_tsv(
  Pathway_summary_combined,
  "Combined_Pathway_summary_statistics.tsv"
)

################################################################################
## 1. ADD GENE-LEVEL Z-SCORES ####
##
## z-score is calculated from log2FoldChange across all pathway DEGs
################################################################################

Pathway_DEGs_combined <- Pathway_DEGs_combined %>%
  dplyr::mutate(
    
    z_log2FC = as.numeric(
      scale(log2FoldChange)
    )
    
  )

Pathway_DEGs_combined %>%
  dplyr::select(
    Collection,
    Pathway,
    gene_name,
    log2FoldChange,
    z_log2FC,
    pvalue,
    padj,
    Direction
  ) %>%
  head()

################################################################################
## 2. GENE-LEVEL PATHWAY STATISTICS TABLE
################################################################################

Gene_table_combined <- Pathway_DEGs_combined %>%
  dplyr::select(
    
    Collection,
    Pathway,
    gene_name,
    
    log2FoldChange,
    z_log2FC,
    
    pvalue,
    padj,
    
    Direction
    
  ) %>%
  dplyr::arrange(
    Collection,
    Pathway,
    padj
  )

readr::write_tsv(
  Gene_table_combined,
  "Combined_Pathway_gene_statistics_with_zscores.tsv"
)

################################################################################
## 3. PATHWAY-LEVEL SUMMARY WITH Z-SCORES
################################################################################

Pathway_summary_combined <- Pathway_DEGs_combined %>%
  
  dplyr::group_by(
    Collection,
    Pathway
  ) %>%
  
  dplyr::summarise(
    
    n_genes =
      dplyr::n(),
    
    n_up_CyMt =
      sum(
        log2FoldChange > 0,
        na.rm = TRUE
      ),
    
    n_up_WT =
      sum(
        log2FoldChange < 0,
        na.rm = TRUE
      ),
    
    percent_up_CyMt =
      round(
        100 *
          sum(log2FoldChange > 0, na.rm = TRUE) /
          dplyr::n(),
        1
      ),
    
    percent_up_WT =
      round(
        100 *
          sum(log2FoldChange < 0, na.rm = TRUE) /
          dplyr::n(),
        1
      ),
    
    mean_log2FC =
      mean(
        log2FoldChange,
        na.rm = TRUE
      ),
    
    median_log2FC =
      median(
        log2FoldChange,
        na.rm = TRUE
      ),
    
    mean_abs_log2FC =
      mean(
        abs(log2FoldChange),
        na.rm = TRUE
      ),
    
    mean_z_log2FC =
      mean(
        z_log2FC,
        na.rm = TRUE
      ),
    
    median_z_log2FC =
      median(
        z_log2FC,
        na.rm = TRUE
      ),
    
    min_padj =
      min(
        padj,
        na.rm = TRUE
      ),
    
    median_padj =
      median(
        padj,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) %>%
  
  dplyr::arrange(
    Collection,
    median_padj
  )

readr::write_tsv(
  Pathway_summary_combined,
  "Combined_Pathway_summary_statistics_with_zscores.tsv"
)

################################################################################
## ADD z-SCORES TO ALL DEG TABLES
##
## z-score is calculated from log2FoldChange across all DEGs
## within each analysis separately
################################################################################

DEGs_JOC047 <- DEGs_JOC047 %>%
  dplyr::mutate(
    z_log2FC = as.numeric(
      scale(log2FoldChange)
    )
  )

DEGs_JOC064 <- DEGs_JOC064 %>%
  dplyr::mutate(
    z_log2FC = as.numeric(
      scale(log2FoldChange)
    )
  )

DEGs_combined <- DEGs_combined %>%
  dplyr::mutate(
    z_log2FC = as.numeric(
      scale(log2FoldChange)
    )
  )


################################################################################
## INSTALL IF NEEDED
################################################################################

if (!requireNamespace("PoiClaClu", quietly = TRUE))
  BiocManager::install("PoiClaClu")

library(PoiClaClu)
library(pheatmap)

################################################################################
## POISSON DISTANCE HEATMAP FUNCTION
## Uses RAW INTEGER COUNTS
################################################################################

plot_poisson_heatmap <- function(
    dds_object,
    metadata,
    title
) {
  
  ##############################################################################
  ## 1. Extract raw integer counts
  ##
  ## These are the original counts stored in the DESeq2 object.
  ## No VST, no normalization, no batch correction.
  ## Hierarchical clustering on both axis
  ##############################################################################
  
  count_matrix <- DESeq2::counts(
    dds_object,
    normalized = FALSE
  )
  
  
  ##############################################################################
  ## 2. Confirm integer count structure
  ##############################################################################
  
  storage.mode(count_matrix) <- "integer"
  
  
  ##############################################################################
  ## 3. Calculate Poisson distance between samples
  ##
  ## PoissonDistance expects:
  ## rows    = samples
  ## columns = genes
  ##
  ## Therefore transpose the count matrix.
  ##############################################################################
  
  pd <- PoiClaClu::PoissonDistance(
    t(count_matrix)
  )
  
  
  ##############################################################################
  ## 4. Convert distance object to matrix for heatmap display
  ##############################################################################
  
  sample_dist_matrix <- as.matrix(
    pd$dd
  )
  
  
  ##############################################################################
  ## 5. Add sample names explicitly
  ##############################################################################
  
  rownames(sample_dist_matrix) <-
    colnames(count_matrix)
  
  colnames(sample_dist_matrix) <-
    colnames(count_matrix)
  
  
  ##############################################################################
  ## 6. Build sample annotation
  ##############################################################################
  
  if ("experiment" %in% colnames(metadata)) {
    
    annotation <- data.frame(
      
      Condition =
        metadata$condition,
      
      Experiment =
        metadata$experiment
      
    )
    
  } else {
    
    annotation <- data.frame(
      
      Condition =
        metadata$condition
      
    )
    
  }
  
  
  rownames(annotation) <-
    rownames(metadata)
  
  
  ##############################################################################
  ## 7. Match annotation order to distance matrix
  ##############################################################################
  
  annotation <- annotation[
    rownames(sample_dist_matrix),
    ,
    drop = FALSE
  ]
  
  
  ##############################################################################
  ## 8. Verify sample names match
  ##############################################################################
  
  stopifnot(
    identical(
      rownames(annotation),
      rownames(sample_dist_matrix)
    )
  )
  
  
  ##############################################################################
  ## 9. Draw Poisson distance heatmap
  ##
  ## pd$dd = Poisson distance object used for hierarchical clustering
  ## sample_dist_matrix = values displayed in the heatmap
  ##############################################################################
  
  p <- pheatmap::pheatmap(
    
    sample_dist_matrix,
    
    clustering_distance_rows =
      pd$dd,
    
    clustering_distance_cols =
      pd$dd,
    
    clustering_method =
      "complete",
    
    annotation_row =
      annotation,
    
    annotation_col =
      annotation,
    
    show_rownames =
      TRUE,
    
    show_colnames =
      TRUE,
    
    fontsize =
      10,
    
    border_color =
      NA,
    
    main =
      title
    
  )
  
  
  ##############################################################################
  ## 10. Return useful objects
  ##############################################################################
  
  return(
    list(
      
      heatmap =
        p,
      
      poisson_distance =
        pd$dd,
      
      distance_matrix =
        sample_dist_matrix,
      
      raw_counts =
        count_matrix
      
    )
  )
}



plot_poisson_heatmap(
  
  dds_JOC047,
  
  as.data.frame(colData(dds_JOC047)),
  
  "JOC047 Poisson Distance"
  
)

plot_poisson_heatmap(
  
  dds_JOC064,
  
  as.data.frame(colData(dds_JOC064)),
  
  "JOC064 Poisson Distance"
  
)

plot_poisson_heatmap(
  
  dds_combined,
  
  as.data.frame(colData(dds_combined)),
  
  "Combined Poisson Distance"
  
)

################################################################################
## SAVE POISSON DISTANCE HEATMAP
################################################################################

pdf(
  "Combined_PoissonDistanceHeatmap.pdf",
  width = 8,
  height = 8
)

plot_poisson_heatmap(
  
  dds_combined,
  
  as.data.frame(colData(dds_combined)),
  
  "Combined Poisson Distance"
  
)

dev.off()

pdf(
  "JOC064_PoissonDistanceHeatmap.pdf",
  width = 8,
  height = 8
)

plot_poisson_heatmap(
  
  dds_JOC064,
  
  as.data.frame(colData(dds_JOC064)),
  
  "JOC064 Poisson Distance"
  
)

dev.off()

pdf(
  "JOC047_PoissonDistanceHeatmap.pdf",
  width = 8,
  height = 8
)

plot_poisson_heatmap(
  
  dds_JOC047,
  
  as.data.frame(colData(dds_JOC047)),
  
  "JOC047 Poisson Distance"
  
)

dev.off()


################################################################################
## CONSTRAINED POISSON DISTANCE HEATMAP — COMBINED DATASET
##
## ROWS / Y AXIS:
##   WT
##      hierarchically clustered within WT
##   CyMt
##      hierarchically clustered within CyMt
##
## COLUMNS / X AXIS:
##   JOC047
##      WT
##      CyMt
##   JOC064
##      WT
##      CyMt
##
## Within each subgroup, samples are hierarchically clustered using
## Poisson distance and complete linkage.
################################################################################

library(DESeq2)
library(PoiClaClu)
library(pheatmap)
library(dplyr)


################################################################################
## 1. GET RAW INTEGER COUNTS
################################################################################

count_matrix <- DESeq2::counts(
  dds_combined,
  normalized = FALSE
)

storage.mode(count_matrix) <- "integer"


################################################################################
## 2. CALCULATE POISSON DISTANCE
################################################################################

pd <- PoiClaClu::PoissonDistance(
  t(count_matrix)
)

poisson_dist_matrix <- as.matrix(
  pd$dd
)


################################################################################
## 3. ADD SAMPLE NAMES
################################################################################

rownames(poisson_dist_matrix) <-
  colnames(count_matrix)

colnames(poisson_dist_matrix) <-
  colnames(count_matrix)


################################################################################
## 4. GET COMBINED METADATA
################################################################################

metadata_poisson <- as.data.frame(
  DESeq2::colData(dds_combined)
)

metadata_poisson$condition <- factor(
  metadata_poisson$condition,
  levels = c(
    "WT",
    "CyMt"
  )
)

metadata_poisson$experiment <- factor(
  metadata_poisson$experiment,
  levels = c(
    "JOC047",
    "JOC064"
  )
)


################################################################################
## 5. HELPER FUNCTION
##
## Hierarchically cluster ONLY the samples supplied to the function.
################################################################################

cluster_within_group <- function(
    sample_names,
    distance_matrix
) {
  
  # If subgroup contains only 0 or 1 sample,
  # there is nothing to cluster.
  if (length(sample_names) <= 1) {
    return(sample_names)
  }
  
  subgroup_matrix <- distance_matrix[
    sample_names,
    sample_names,
    drop = FALSE
  ]
  
  subgroup_dist <- as.dist(
    subgroup_matrix
  )
  
  hc <- hclust(
    subgroup_dist,
    method = "complete"
  )
  
  sample_names[
    hc$order
  ]
}






















################################################################################
## SAVE COMBINED RESULTS for IPA analysis 08312026 ####
################################################################################

res_combined_df <- as.data.frame(res_combined)

res_combined_df$Geneid <- rownames(res_combined_df)

res_combined_df <- res_combined_df %>%
  dplyr::left_join(
    gene_anno_filtered %>%
      dplyr::select(
        Geneid,
        gene_name
      ),
    by = "Geneid"
  ) %>%
  dplyr::relocate(
    Geneid,
    gene_name
  )

readr::write_CSV(
  res_combined_df,
  "Combined_DESeq2_results.tsv"
)

write.csv(
  res_combined_df,
  file = "Combined_DESeq2_results.csv",
  row.names = FALSE
)





################################################################################
## ADD z-SCORE TO COMBINED DESeq2 RESULTS
##
## z-score = standardized log2FoldChange
################################################################################

res_combined_df <- res_combined_df %>%
  dplyr::mutate(
    z_score = as.numeric(
      scale(log2FoldChange)
    )
  ) %>%
  dplyr::relocate(
    z_score,
    .after = log2FoldChange
  )

################################################################################
## SAVE AS CSV
################################################################################

write.csv(
  res_combined_df,
  file = "Combined_DESeq2_results_with_zscores.csv",
  row.names = FALSE
)

################################################################################
# 09092026 code for manuscript figures ####
################################################################################
## PC1 / PC2 LOADINGS → GO / KEGG → GSVA
##
## INPUT EXPRESSION MATRIX:
##     VST_symbols_combined_batchCorrected
##
## Because this matrix is already batch-corrected, downstream GSVA statistics
## use:
##
##     ~ condition
##
## NOT:
##
##     ~ experiment + condition
################################################################################
################################################################################

################################################################################
## 1. PACKAGES
################################################################################

library(dplyr)
library(tidyr)
library(tibble)

library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)

library(GSVA)
library(limma)

library(ggplot2)
library(pheatmap)

################################################################################
## 2. VERIFY REQUIRED OBJECTS
################################################################################

required_objects <- c(
  "PCA_combined_batchCorrected",
  "VST_symbols_combined_batchCorrected",
  "metadata_combined",
  "gene_anno_filtered"
)

sapply(required_objects, exists)

# Confirm matrix/sample matching
dim(VST_symbols_combined_batchCorrected)

head(rownames(VST_symbols_combined_batchCorrected))

colnames(VST_symbols_combined_batchCorrected)

head(metadata_combined)

identical(
  colnames(VST_symbols_combined_batchCorrected),
  rownames(metadata_combined)
)

################################################################################
## 3. EXTRACT PCA LOADINGS
################################################################################

pca_loadings <- as.data.frame(
  PCA_combined_batchCorrected$PCA$rotation
) %>%
  tibble::rownames_to_column(
    "Geneid"
  )
# Join Ensembl IDs to gene_anno_filtered to add gene symbols
pca_loadings <- pca_loadings %>%
  dplyr::left_join(
    gene_anno_filtered %>%
      dplyr::select(
        Geneid,
        gene_name
      ) %>%
      dplyr::distinct(
        Geneid,
        .keep_all = TRUE
      ),
    by = "Geneid"
  )

# Check if any genes failed to map/annotate 
sum(is.na(pca_loadings$gene_name))

sum(!is.na(pca_loadings$gene_name))

head(
  pca_loadings[
    ,
    c("Geneid", "gene_name", "PC1", "PC2")
  ]
)

################################################################################
## SELECT TOP PC1 / PC2 LOADING GENES
################################################################################

top_n_loadings <- 100

PC1_positive <- pca_loadings %>%
  dplyr::arrange(dplyr::desc(PC1)) %>%
  dplyr::slice_head(n = top_n_loadings)

PC1_negative <- pca_loadings %>%
  dplyr::arrange(PC1) %>%
  dplyr::slice_head(n = top_n_loadings)

PC2_positive <- pca_loadings %>%
  dplyr::arrange(dplyr::desc(PC2)) %>%
  dplyr::slice_head(n = top_n_loadings)

PC2_negative <- pca_loadings %>%
  dplyr::arrange(PC2) %>%
  dplyr::slice_head(n = top_n_loadings)

head(
  PC1_positive %>%
    dplyr::select(Geneid, gene_name, PC1)
)

head(
  PC1_negative %>%
    dplyr::select(Geneid, gene_name, PC1)
)

head(
  PC2_positive %>%
    dplyr::select(Geneid, gene_name, PC2)
)

head(
  PC2_negative %>%
    dplyr::select(Geneid, gene_name, PC2)
)

################################################################################
## 6. CONVERT TOP PC LOADING GENES: SYMBOL → ENTREZ ID
################################################################################

library(AnnotationDbi)
library(org.Mm.eg.db)
library(dplyr)

symbols_to_entrez <- function(symbols) {
  
  AnnotationDbi::select(
    org.Mm.eg.db,
    keys = unique(symbols),
    keytype = "SYMBOL",
    columns = c("SYMBOL", "ENTREZID")
  ) %>%
    dplyr::filter(!is.na(ENTREZID)) %>%
    dplyr::distinct(SYMBOL, ENTREZID)
}


PC1_positive_map <- symbols_to_entrez(
  PC1_positive$gene_name
)

PC1_negative_map <- symbols_to_entrez(
  PC1_negative$gene_name
)

PC2_positive_map <- symbols_to_entrez(
  PC2_positive$gene_name
)

PC2_negative_map <- symbols_to_entrez(
  PC2_negative$gene_name
)


################################################################################
## DEFINE GO ENRICHMENT FUNCTION
################################################################################

library(clusterProfiler)
library(org.Mm.eg.db)

run_GO_loading_enrichment <- function(entrez_ids) {
  
  clusterProfiler::enrichGO(
    gene = unique(entrez_ids),
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20,
    readable = TRUE
  )
}

GO_PC1_positive <- run_GO_loading_enrichment(
  PC1_positive_map$ENTREZID
)

GO_PC1_negative <- run_GO_loading_enrichment(
  PC1_negative_map$ENTREZID
)

GO_PC2_positive <- run_GO_loading_enrichment(
  PC2_positive_map$ENTREZID
)

GO_PC2_negative <- run_GO_loading_enrichment(
  PC2_negative_map$ENTREZID
)

################################################################################
## RUN GO BP ENRICHMENT
################################################################################

GO_PC1_positive <- run_GO_loading_enrichment(
  PC1_positive_map$ENTREZID
)

GO_PC1_negative <- run_GO_loading_enrichment(
  PC1_negative_map$ENTREZID
)

GO_PC2_positive <- run_GO_loading_enrichment(
  PC2_positive_map$ENTREZID
)

GO_PC2_negative <- run_GO_loading_enrichment(
  PC2_negative_map$ENTREZID
)

#check results
c(
  GO_PC1_positive = nrow(as.data.frame(GO_PC1_positive)),
  GO_PC1_negative = nrow(as.data.frame(GO_PC1_negative)),
  GO_PC2_positive = nrow(as.data.frame(GO_PC2_positive)),
  GO_PC2_negative = nrow(as.data.frame(GO_PC2_negative))
)

################################################################################
## INSPECT TOP GO BP TERMS
################################################################################

get_top_GO <- function(go_object, n = 20) {
  
  as.data.frame(go_object) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::select(
      ID,
      Description,
      GeneRatio,
      BgRatio,
      pvalue,
      p.adjust,
      qvalue,
      Count,
      geneID
    ) %>%
    dplyr::slice_head(n = n)
}

################################################################################
## INSPECT TOP GO BP TERMS
################################################################################

get_top_GO <- function(go_object, n = 20) {
  
  as.data.frame(go_object) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::select(
      ID,
      Description,
      GeneRatio,
      BgRatio,
      pvalue,
      p.adjust,
      qvalue,
      Count,
      geneID
    ) %>%
    dplyr::slice_head(n = n)
}


GO_PC1_positive_top20 <- get_top_GO(GO_PC1_positive, 20)
GO_PC1_negative_top20 <- get_top_GO(GO_PC1_negative, 20)
GO_PC2_positive_top20 <- get_top_GO(GO_PC2_positive, 20)
GO_PC2_negative_top20 <- get_top_GO(GO_PC2_negative, 20)

View(GO_PC1_positive_top20)
View(GO_PC1_negative_top20)
View(GO_PC2_positive_top20)
View(GO_PC2_negative_top20)

# Make GO plots
library(enrichplot)
library(ggplot2)

dotplot(GO_PC1_positive, showCategory = 20) +
  ggtitle("PC1 Positive Loadings — GO BP")

dotplot(GO_PC1_negative, showCategory = 20) +
  ggtitle("PC1 Negative Loadings — GO BP")

dotplot(GO_PC2_positive, showCategory = 20) +
  ggtitle("PC2 Positive Loadings — GO BP")

dotplot(GO_PC2_negative, showCategory = 20) +
  ggtitle("PC2 Negative Loadings — GO BP")



#Run sample PCA scores to actually see if the actual PC1/PC2 leans more strongly to WT or CyMt
pca_scores <- as.data.frame(
  PCA_combined_batchCorrected$PCA$x
)

pca_scores$condition <- metadata_combined[
  rownames(pca_scores),
  "condition"
]

pca_scores %>%
  dplyr::group_by(condition) %>%
  dplyr::summarise(
    mean_PC1 = mean(PC1),
    mean_PC2 = mean(PC2),
    median_PC1 = median(PC1),
    median_PC2 = median(PC2)
  )

# Adding PCA association to show if gene expression pattern is associated to WT or CyMt more closely. 
# DOES NOT MEANT UPREGULATE OR DOWNREGULATED only that there is an association to the condition/group

PC_direction_key <- data.frame(
  PC_group = c(
    "PC1_positive",
    "PC1_negative",
    "PC2_positive",
    "PC2_negative"
  ),
  PCA_association = c(
    "CyMt-associated",
    "WT-associated",
    "WT-associated",
    "CyMt-associated"
  )
)

PC_direction_key

#####END of GO analysiss

#####START of KEGG analysis

################################################################################
## KEGG ENRICHMENT OF PC1 / PC2 LOADING GENES
################################################################################

library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(dplyr)
library(ggplot2)
library(enrichplot)

run_KEGG_loading_enrichment <- function(entrez_ids) {
  
  clusterProfiler::enrichKEGG(
    gene = unique(entrez_ids),
    organism = "mmu",
    keyType = "ncbi-geneid",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )
}

exists("run_KEGG_loading_enrichment")
################################################################################
## RUN KEGG ENRICHMENT
################################################################################

# PC1 positive = CyMt-associated
KEGG_PC1_positive <- run_KEGG_loading_enrichment(
  PC1_positive_map$ENTREZID
)

# PC1 negative = WT-associated
KEGG_PC1_negative <- run_KEGG_loading_enrichment(
  PC1_negative_map$ENTREZID
)

# PC2 positive = WT-associated
KEGG_PC2_positive <- run_KEGG_loading_enrichment(
  PC2_positive_map$ENTREZID
)

# PC2 negative = CyMt-associated
KEGG_PC2_negative <- run_KEGG_loading_enrichment(
  PC2_negative_map$ENTREZID
)

################################################################################
## NUMBER OF KEGG PATHWAYS
################################################################################

c(
  KEGG_PC1_positive = nrow(as.data.frame(KEGG_PC1_positive)),
  KEGG_PC1_negative = nrow(as.data.frame(KEGG_PC1_negative)),
  KEGG_PC2_positive = nrow(as.data.frame(KEGG_PC2_positive)),
  KEGG_PC2_negative = nrow(as.data.frame(KEGG_PC2_negative))
)

################################################################################
## EXTRACT TOP KEGG PATHWAYS
################################################################################

get_top_KEGG <- function(kegg_object, n = 20) {
  
  as.data.frame(kegg_object) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::select(
      ID,
      Description,
      GeneRatio,
      BgRatio,
      pvalue,
      p.adjust,
      qvalue,
      Count,
      geneID
    ) %>%
    dplyr::slice_head(n = n)
}

KEGG_PC1_positive_top20 <- get_top_KEGG(
  KEGG_PC1_positive,
  20
)

KEGG_PC1_negative_top20 <- get_top_KEGG(
  KEGG_PC1_negative,
  20
)

KEGG_PC2_positive_top20 <- get_top_KEGG(
  KEGG_PC2_positive,
  20
)

KEGG_PC2_negative_top20 <- get_top_KEGG(
  KEGG_PC2_negative,
  20
)

KEGG_PC1_positive_top20 %>%
  dplyr::select(
    ID,
    Description,
    GeneRatio,
    p.adjust,
    Count
  )

KEGG_PC1_negative_top20 %>%
  dplyr::select(
    ID,
    Description,
    GeneRatio,
    p.adjust,
    Count
  )

KEGG_PC2_positive_top20 %>%
  dplyr::select(
    ID,
    Description,
    GeneRatio,
    p.adjust,
    Count
  )

KEGG_PC2_negative_top20 %>%
  dplyr::select(
    ID,
    Description,
    GeneRatio,
    p.adjust,
    Count
  )

KEGG_PC1_positive_top20 <- KEGG_PC1_positive_top20 %>%
  mutate(
    PC = "PC1",
    Loading_direction = "Positive",
    PCA_association = "CyMt-associated"
  )

KEGG_PC1_negative_top20 <- KEGG_PC1_negative_top20 %>%
  mutate(
    PC = "PC1",
    Loading_direction = "Negative",
    PCA_association = "WT-associated"
  )

KEGG_PC2_positive_top20 <- KEGG_PC2_positive_top20 %>%
  mutate(
    PC = "PC2",
    Loading_direction = "Positive",
    PCA_association = "WT-associated"
  )

KEGG_PC2_negative_top20 <- KEGG_PC2_negative_top20 %>%
  mutate(
    PC = "PC2",
    Loading_direction = "Negative",
    PCA_association = "CyMt-associated"
  )

KEGG_PC_loading_summary <- bind_rows(
  KEGG_PC1_positive_top20,
  KEGG_PC1_negative_top20,
  KEGG_PC2_positive_top20,
  KEGG_PC2_negative_top20
)

KEGG_PC_loading_summary %>%
  dplyr::select(
    PC,
    Loading_direction,
    PCA_association,
    ID,
    Description,
    GeneRatio,
    Count,
    p.adjust
  )

write.csv(
  KEGG_PC_loading_summary,
  "Combined_PCA_PC1_PC2_KEGG_top_pathways.csv",
  row.names = FALSE
)

dotplot(
  KEGG_PC1_positive,
  showCategory = 20
) +
  ggtitle(
    "PC1 Positive Loadings — KEGG\nCyMt-associated"
  )

dotplot(
  KEGG_PC1_negative,
  showCategory = 20
) +
  ggtitle(
    "PC1 Negative Loadings — KEGG\nWT-associated"
  )

dotplot(
  KEGG_PC2_positive,
  showCategory = 20
) +
  ggtitle(
    "PC2 Positive Loadings — KEGG\nWT-associated"
  )

dotplot(
  KEGG_PC2_negative,
  showCategory = 20
) +
  ggtitle(
    "PC2 Negative Loadings — KEGG\nCyMt-associated"
  )

################################################################################
################################################################################
## Combine GO and KEGG to run GSVA
##
## INPUT:
##   VST_symbols_combined_batchCorrected
##
## DOWNSTREAM MODEL:
##   ~ condition
################################################################################
################################################################################

library(dplyr)
library(tidyr)
library(tibble)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(GSVA)
library(limma)
library(pheatmap)

################################################################################
## 1. SELECT TOP ENRICHED TERMS
################################################################################

get_top_enrichment_terms <- function(
    enrichment_object,
    n_terms = 10,
    padj_cutoff = 0.05
) {
  
  as.data.frame(enrichment_object) %>%
    dplyr::filter(
      !is.na(p.adjust),
      p.adjust < padj_cutoff
    ) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = n_terms)
}

################################################################################
## 2. TOP GO TERMS
################################################################################

GO_PC1_positive_top <- get_top_enrichment_terms(
  GO_PC1_positive, 10
)

GO_PC1_negative_top <- get_top_enrichment_terms(
  GO_PC1_negative, 10
)

GO_PC2_positive_top <- get_top_enrichment_terms(
  GO_PC2_positive, 10
)

GO_PC2_negative_top <- get_top_enrichment_terms(
  GO_PC2_negative, 10
)

################################################################################
## 3. TOP KEGG TERMS
################################################################################

KEGG_PC1_positive_top <- get_top_enrichment_terms(
  KEGG_PC1_positive, 10
)

KEGG_PC1_negative_top <- get_top_enrichment_terms(
  KEGG_PC1_negative, 10
)

KEGG_PC2_positive_top <- get_top_enrichment_terms(
  KEGG_PC2_positive, 10
)

KEGG_PC2_negative_top <- get_top_enrichment_terms(
  KEGG_PC2_negative, 10
)

################################################################################
## 4. CONVERT GO RESULTS TO GENE SET LISTS
################################################################################

GO_result_to_gene_sets <- function(
    result_df,
    prefix
) {
  
  if (nrow(result_df) == 0) {
    return(list())
  }
  
  sets <- lapply(
    seq_len(nrow(result_df)),
    function(i) {
      
      unique(
        strsplit(
          result_df$geneID[i],
          "/",
          fixed = TRUE
        )[[1]]
      )
    }
  )
  
  names(sets) <- paste0(
    prefix,
    "__",
    result_df$ID,
    "__",
    result_df$Description
  )
  
  sets
}

GO_sets_PC1_positive <- GO_result_to_gene_sets(
  GO_PC1_positive_top,
  "GO_PC1_POS"
)

GO_sets_PC1_negative <- GO_result_to_gene_sets(
  GO_PC1_negative_top,
  "GO_PC1_NEG"
)

GO_sets_PC2_positive <- GO_result_to_gene_sets(
  GO_PC2_positive_top,
  "GO_PC2_POS"
)

GO_sets_PC2_negative <- GO_result_to_gene_sets(
  GO_PC2_negative_top,
  "GO_PC2_NEG"
)
################################################################################
## 5. CONVERT KEGG RESULTS TO GENE SYMBOL SETS
################################################################################

KEGG_result_to_gene_sets <- function(
    result_df,
    prefix
) {
  
  if (nrow(result_df) == 0) {
    return(list())
  }
  
  sets <- lapply(
    seq_len(nrow(result_df)),
    function(i) {
      
      entrez_ids <- strsplit(
        result_df$geneID[i],
        "/",
        fixed = TRUE
      )[[1]]
      
      mapping <- AnnotationDbi::select(
        org.Mm.eg.db,
        keys = unique(entrez_ids),
        keytype = "ENTREZID",
        columns = "SYMBOL"
      )
      
      unique(
        na.omit(
          mapping$SYMBOL
        )
      )
    }
  )
  
  names(sets) <- paste0(
    prefix,
    "__",
    result_df$ID,
    "__",
    result_df$Description
  )
  
  sets
}

KEGG_sets_PC1_positive <- KEGG_result_to_gene_sets(
  KEGG_PC1_positive_top,
  "KEGG_PC1_POS"
)

KEGG_sets_PC1_negative <- KEGG_result_to_gene_sets(
  KEGG_PC1_negative_top,
  "KEGG_PC1_NEG"
)

KEGG_sets_PC2_positive <- KEGG_result_to_gene_sets(
  KEGG_PC2_positive_top,
  "KEGG_PC2_POS"
)

KEGG_sets_PC2_negative <- KEGG_result_to_gene_sets(
  KEGG_PC2_negative_top,
  "KEGG_PC2_NEG"
)

################################################################################
## 6. COMBINE GO + KEGG GENE SETS
################################################################################

PC_GSVA_gene_sets <- c(
  GO_sets_PC1_positive,
  GO_sets_PC1_negative,
  GO_sets_PC2_positive,
  GO_sets_PC2_negative,
  
  KEGG_sets_PC1_positive,
  KEGG_sets_PC1_negative,
  KEGG_sets_PC2_positive,
  KEGG_sets_PC2_negative
)

PC_GSVA_gene_sets <- PC_GSVA_gene_sets[
  !duplicated(names(PC_GSVA_gene_sets))
]

################################################################################
## 7. MATCH GENE SETS TO VST MATRIX
################################################################################

PC_GSVA_gene_sets_filtered <- lapply(
  PC_GSVA_gene_sets,
  function(x) {
    
    intersect(
      unique(x),
      rownames(VST_symbols_combined_batchCorrected)
    )
  }
)

#Keep pathways with atleast 5 genes present
PC_GSVA_gene_sets_filtered <-
  PC_GSVA_gene_sets_filtered[
    lengths(PC_GSVA_gene_sets_filtered) >= 5
  ]

#inspect pathway coverage
GSVA_gene_set_coverage <- data.frame(
  Pathway = names(PC_GSVA_gene_sets),
  n_original = lengths(PC_GSVA_gene_sets),
  n_in_VST = sapply(
    PC_GSVA_gene_sets,
    function(x) {
      length(
        intersect(
          x,
          rownames(VST_symbols_combined_batchCorrected)
        )
      )
    }
  )
)

GSVA_gene_set_coverage %>%
  dplyr::arrange(n_in_VST)

################################################################################
## 8. RUN GSVA
################################################################################

gsva_param_PC <- GSVA::gsvaParam(
  VST_symbols_combined_batchCorrected,
  PC_GSVA_gene_sets_filtered
)

GSVA_PC_scores <- GSVA::gsva(
  gsva_param_PC,
  verbose = FALSE
)

dim(GSVA_PC_scores)

GSVA_PC_scores[
  1:min(5, nrow(GSVA_PC_scores)),
  1:min(5, ncol(GSVA_PC_scores))
]
################################################################################
## 9. ALIGN METADATA
################################################################################

metadata_combined <- metadata_combined[
  colnames(GSVA_PC_scores),
  ,
  drop = FALSE
]

identical(
  colnames(GSVA_PC_scores),
  rownames(metadata_combined)
)

#Set WT as reference
metadata_combined$condition <- factor(
  metadata_combined$condition,
  levels = c("WT", "CyMt")
)

################################################################################
## 10. LIMMA: CyMt vs WT
################################################################################

GSVA_design <- model.matrix(
  ~ condition,
  data = metadata_combined
)

colnames(GSVA_design)

GSVA_fit <- limma::lmFit(
  GSVA_PC_scores,
  GSVA_design
)

GSVA_fit <- limma::eBayes(
  GSVA_fit
)

GSVA_PC_results <- limma::topTable(
  GSVA_fit,
  coef = "conditionCyMt",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
) %>%
  tibble::rownames_to_column("Pathway")

################################################################################
## 11. ADD PATHWAY SOURCE / PC / DIRECTION
################################################################################

GSVA_PC_results <- GSVA_PC_results %>%
  dplyr::mutate(
    
    Source = dplyr::case_when(
      grepl("^GO_", Pathway) ~ "GO",
      grepl("^KEGG_", Pathway) ~ "KEGG",
      TRUE ~ NA_character_
    ),
    
    PC = dplyr::case_when(
      grepl("PC1", Pathway) ~ "PC1",
      grepl("PC2", Pathway) ~ "PC2",
      TRUE ~ NA_character_
    ),
    
    Loading_direction = dplyr::case_when(
      grepl("_POS", Pathway) ~ "Positive loading",
      grepl("_NEG", Pathway) ~ "Negative loading",
      TRUE ~ NA_character_
    ),
    
    PCA_association = dplyr::case_when(
      grepl("PC1_POS", Pathway) ~ "CyMt-associated",
      grepl("PC1_NEG", Pathway) ~ "WT-associated",
      grepl("PC2_POS", Pathway) ~ "WT-associated",
      grepl("PC2_NEG", Pathway) ~ "CyMt-associated",
      TRUE ~ NA_character_
    ),
    
    GSVA_direction = dplyr::case_when(
      logFC > 0 ~ "Higher GSVA score in CyMt",
      logFC < 0 ~ "Higher GSVA score in WT",
      TRUE ~ "No difference"
    )
  )

GSVA_PC_results %>%
  dplyr::select(
    Pathway,
    Source,
    PC,
    Loading_direction,
    PCA_association,
    logFC,
    AveExpr,
    t,
    P.Value,
    adj.P.Val,
    GSVA_direction
  )

#significant pathways
GSVA_PC_significant <- GSVA_PC_results %>%
  dplyr::filter(
    adj.P.Val < 0.05
  )

write.csv(
  GSVA_PC_results,
  "Combined_PC1_PC2_GO_KEGG_GSVA_results.csv",
  row.names = FALSE
)

write.csv(
  GSVA_PC_significant,
  "Combined_PC1_PC2_GO_KEGG_GSVA_significant.csv",
  row.names = FALSE
)


######Analysis of imported xcel file from Qiagen Ingenuity Pathway Analysis ######

################################################################################
## IMPORT QIAGEN IPA CANONICAL PATHWAY RESULTS
################################################################################

library(readxl)
library(dplyr)
library(stringr)

# Select your IPA Excel file
IPA_raw <- readxl::read_excel(
  file.choose()
)

# Inspect column names
colnames(IPA_raw)

# Look at first few rows
head(IPA_raw)

################################################################################
## COUNT MOLECULES IN EACH CANONICAL PATHWAY
################################################################################

IPA_filtered <- IPA_raw %>%
  dplyr::mutate(
    
    Molecule_Count = dplyr::case_when(
      
      is.na(Molecules) ~ 0L,
      
      trimws(Molecules) == "" ~ 0L,
      
      TRUE ~ stringr::str_count(Molecules, ",") + 1L
    )
    
  )

#Inspect counts
IPA <- IPA_raw %>%
  dplyr::rename(
    Pathway = `Ingenuity Canonical Pathways`,
    neg_log10_pvalue = `-log(p-value)`,
    Ratio = Ratio,
    z_score = `z-score`,
    Molecules = Molecules
  )

IPA <- IPA %>%
  dplyr::mutate(
    Molecule_Count = dplyr::case_when(
      is.na(Molecules) ~ 0L,
      trimws(Molecules) == "" ~ 0L,
      TRUE ~ stringr::str_count(Molecules, ",") + 1L
    )
  )

IPA %>%
  dplyr::select(
    Pathway,
    neg_log10_pvalue,
    Ratio,
    z_score,
    Molecule_Count
  ) %>%
  dplyr::arrange(
    dplyr::desc(Molecule_Count)
  )

IPA_selected <- IPA %>%
  dplyr::filter(
    neg_log10_pvalue >= 1.3,   # approximately p <= 0.05
    !is.na(z_score),
    abs(z_score) >= 2
  )


################################################################################
## COMPLETE HALLMARK + IMMUNE GSEA GENE TABLES
## Includes ALL leading-edge genes with DESeq2 statistics
## Combined analysis: CyMt vs WT
################################################################################

library(dplyr)
library(tidyr)
library(tibble)

if (!"gene_name" %in% colnames(res_combined_df)) {
  
  res_combined_df <- res_combined_df %>%
    dplyr::left_join(
      gene_anno_filtered %>%
        dplyr::select(Geneid, gene_name) %>%
        dplyr::distinct(Geneid, .keep_all = TRUE),
      by = "Geneid"
    )
}

make_complete_GSEA_table <- function(
    gsea_object,
    deseq_results,
    collection_name
) {
  
  gsea_df <- as.data.frame(gsea_object)
  
  output <- gsea_df %>%
    dplyr::select(
      ID,
      Description,
      setSize,
      enrichmentScore,
      NES,
      pvalue,
      p.adjust,
      qvalue,
      rank,
      core_enrichment
    ) %>%
    
    tidyr::separate_rows(
      core_enrichment,
      sep = "/"
    ) %>%
    
    dplyr::rename(
      gene_name = core_enrichment,
      GSEA_pvalue = pvalue,
      GSEA_padj = p.adjust,
      GSEA_qvalue = qvalue
    ) %>%
    
    dplyr::mutate(
      Collection = collection_name,
      
      GSEA_direction = dplyr::case_when(
        NES > 0 ~ "CyMt enriched",
        NES < 0 ~ "WT enriched",
        TRUE ~ NA_character_
      )
    ) %>%
    
    dplyr::left_join(
      deseq_results %>%
        dplyr::select(
          Geneid,
          gene_name,
          baseMean,
          log2FoldChange,
          lfcSE,
          stat,
          pvalue,
          padj
        ),
      by = "gene_name"
    ) %>%
    
    dplyr::rename(
      DESeq2_pvalue = pvalue,
      DESeq2_padj = padj
    ) %>%
    
    dplyr::mutate(
      DE_direction = dplyr::case_when(
        log2FoldChange > 0 ~ "Higher in CyMt",
        log2FoldChange < 0 ~ "Higher in WT",
        TRUE ~ NA_character_
      )
    ) %>%
    
    dplyr::select(
      Collection,
      ID,
      Description,
      setSize,
      enrichmentScore,
      NES,
      GSEA_direction,
      GSEA_pvalue,
      GSEA_padj,
      GSEA_qvalue,
      rank,
      Geneid,
      gene_name,
      baseMean,
      log2FoldChange,
      lfcSE,
      stat,
      DESeq2_pvalue,
      DESeq2_padj,
      DE_direction
    )
  
  return(output)
}

Hallmark_GSEA_complete_combined <- make_complete_GSEA_table(
  gsea_object = Hallmark_GSEA_combined,
  deseq_results = res_combined_df,
  collection_name = "HALLMARK"
)

Immune_GSEA_complete_combined <- make_complete_GSEA_table(
  gsea_object = Immune_GSEA_combined,
  deseq_results = res_combined_df,
  collection_name = "IMMUNE"
)

GSEA_complete_combined <- dplyr::bind_rows(
  Hallmark_GSEA_complete_combined,
  Immune_GSEA_complete_combined
)

View(Hallmark_GSEA_complete_combined)
View(Immune_GSEA_complete_combined)
View(GSEA_complete_combined)

write.csv(
  Hallmark_GSEA_complete_combined,
  "Combined_HALLMARK_GSEA_complete_leading_edge_genes.csv",
  row.names = FALSE
)

write.csv(
  Immune_GSEA_complete_combined,
  "Combined_IMMUNE_GSEA_complete_leading_edge_genes.csv",
  row.names = FALSE
)

write.csv(
  GSEA_complete_combined,
  "Combined_HALLMARK_IMMUNE_GSEA_complete_leading_edge_genes.csv",
  row.names = FALSE
)


# Save table as compressed TSV
write.table(
  counts_JOC047_raw,
  file = gzfile("counts1.tsv.gz", "w"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)


# Save table as compressed TSV
write.table(
  counts_JOC064_raw,
  file = gzfile("counts2.tsv.gz", "w"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# Save table as compressed TSV
write.table(
  counts_combined_annotated,
  file = gzfile("counts.tsv.gz", "w"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)


########09292026 Volcano plot labeled with only select DEGs########


library(ggplot2)
library(ggrepel)
library(dplyr)

# 1. Enter ONLY the genes you want labeled
genes_to_label <- c(
  "Trem1",
  "Spp1",
  "Lrg1", "Cxcl13", "Ccr1", "Csf1", "Cd38", "Myd88", "Itgam", "Cd84"
)

# 2. Prepare volcano data
volcano_df <- res_combined_df %>%
  dplyr::filter(
    !is.na(padj),
    !is.na(log2FoldChange)
  ) %>%
  dplyr::mutate(
    neg_log10_padj = -log10(pmax(padj, .Machine$double.xmin)),
    
    Significance = dplyr::case_when(
      padj < 0.05 & log2FoldChange > 0.75 ~ "Higher in CyMt",
      padj < 0.05 & log2FoldChange < -0.75 ~ "Higher in WT",
      TRUE ~ "Other genes"
    )
  )

# Add gene symbols if needed
if (!"gene_name" %in% colnames(volcano_df)) {
  volcano_df <- volcano_df %>%
    dplyr::left_join(
      gene_anno_filtered %>%
        dplyr::select(Geneid, gene_name) %>%
        dplyr::distinct(Geneid, .keep_all = TRUE),
      by = "Geneid"
    )
}

# 3. Select only requested genes for labeling
label_df <- volcano_df %>%
  dplyr::filter(gene_name %in% genes_to_label)

# 4. Create volcano plot
volcano_combined <- ggplot(
  volcano_df,
  aes(
    x = log2FoldChange,
    y = neg_log10_padj,
    color = Significance
  )
) +
  geom_point(
    alpha = 0.7,
    size = 2
  ) +
  
  # Horizontal adjusted p-value cutoff
  geom_hline(
    yintercept = -log10(0.05),
    color = "black",
    linetype = "dotted",
    linewidth = 0.7
  ) +
  
  # Vertical log2FC cutoff
  geom_vline(
    xintercept = c(-1, 1),
    color = "black",
    linetype = "dotted",
    linewidth = 0.7
  ) +
  
  # Label only selected genes
  ggrepel::geom_text_repel(
    data = label_df,
    aes(label = gene_name),
    color = "black",
    size = 3.5,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3
  ) +
  
  scale_color_manual(
    values = c(
      "Higher in CyMt" = "red",
      "Higher in WT" = "blue",
      "Other genes" = "grey70"
    )
  ) +
  
  labs(
    title = "CyMt vs WT",
    x = expression(log[2]~Fold~Change),
    y = expression(-log[10]~adjusted~italic(P)),
    color = NULL
  ) +
  
  theme_classic(base_size = 14)

volcano_combined

############ 09302026 #############
library(ggplot2)
library(ggrepel)
library(dplyr)

################################################################################
## VOLCANO PLOT — CyMt vs WT
##
## Significant:
##   padj < 0.05
##   AND |log2FoldChange| > 0.75
##
## Red  = Higher in CyMt
## Blue = Higher in WT
################################################################################

volcano_df <- res_combined_df %>%
  dplyr::filter(
    !is.na(log2FoldChange),
    !is.na(padj)
  ) %>%
  dplyr::mutate(
    
    # Calculate -log10 adjusted p-value
    neg_log10_padj = -log10(
      pmax(padj, .Machine$double.xmin)
    ),
    
    # Classify genes
    Significance = dplyr::case_when(
      
      padj < 0.05 &
        log2FoldChange > 0.75 ~ "Higher in CyMt",
      
      padj < 0.05 &
        log2FoldChange < -0.75 ~ "Higher in WT",
      
      TRUE ~ "Not significant"
    )
  )


################################################################################
## ADD GENE SYMBOLS IF NOT ALREADY PRESENT
################################################################################

if (!"gene_name" %in% colnames(volcano_df)) {
  
  volcano_df <- volcano_df %>%
    dplyr::left_join(
      gene_anno_filtered %>%
        dplyr::select(
          Geneid,
          gene_name
        ) %>%
        dplyr::distinct(
          Geneid,
          .keep_all = TRUE
        ),
      by = "Geneid"
    )
}


################################################################################
## GENES TO LABEL
##
## Labels ALL genes meeting:
## padj < 0.05 AND |log2FC| > 0.75
################################################################################

label_df <- volcano_df %>%
  dplyr::filter(
    padj < 0.05,
    abs(log2FoldChange) > 0.75,
    !is.na(gene_name)
  )


################################################################################
## CHECK NUMBER OF LABELED GENES
################################################################################

nrow(label_df)

label_df %>%
  dplyr::select(
    gene_name,
    log2FoldChange,
    padj,
    Significance
  ) %>%
  dplyr::arrange(padj)


################################################################################
## VOLCANO PLOT
################################################################################

volcano_combined <- ggplot(
  volcano_df,
  aes(
    x = log2FoldChange,
    y = neg_log10_padj,
    color = Significance
  )
) +
  
  # All genes
  geom_point(
    alpha = 0.7,
    size = 2
  ) +
  
  # padj = 0.05
  geom_hline(
    yintercept = -log10(0.05),
    color = "black",
    linetype = "dotted",
    linewidth = 0.7
  ) +
  
  # log2FC = -0.75 and +0.75
  geom_vline(
    xintercept = c(-0.75, 0.75),
    color = "black",
    linetype = "dotted",
    linewidth = 0.7
  ) +
  
  # Label ALL significant genes
  ggrepel::geom_text_repel(
    data = label_df,
    aes(label = gene_name),
    
    color = "black",
    size = 3,
    
    # Do not automatically discard overlapping labels
    max.overlaps = Inf,
    
    # Spacing
    box.padding = 0.4,
    point.padding = 0.25,
    
    # Lines connecting labels to dots
    min.segment.length = 0,
    segment.color = "black",
    segment.size = 0.3,
    
    # Improve label placement
    force = 2,
    max.iter = 20000,
    seed = 123
  ) +
  
  # Dot colors
  scale_color_manual(
    values = c(
      "Higher in CyMt" = "red",
      "Higher in WT" = "blue",
      "Not significant" = "grey70"
    ),
    breaks = c(
      "Higher in CyMt",
      "Higher in WT",
      "Not significant"
    )
  ) +
  
  # Axis labels
  labs(
    title = "CyMt vs WT",
    x = expression(log[2]~Fold~Change),
    y = expression(-log[10]~adjusted~italic(P)),
    color = NULL
  ) +
  
  # Clean formatting
  theme_classic(base_size = 14) +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  )


################################################################################
## DISPLAY
################################################################################

volcano_combined








