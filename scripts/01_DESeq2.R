# 01_DESeq2.R
# Differential expression analysis of rice RNA-seq data
# Drought vs Water

library(DESeq2)
library(ggplot2)
library(pheatmap)

count_file <- "D:/rnaseq_data/counts/gene_counts.txt"
output_dir <- "D:/rnaseq_data/aligned"

# Read the featureCounts output
counts <- read.delim(
  count_file,
  comment.char = "#",
  check.names = FALSE
)

# Keep the six sample count columns
count_matrix <- counts[, 7:12]
rownames(count_matrix) <- counts$Geneid

# Define the experimental groups
condition <- factor(
  c(
    "Water", "Water", "Water",
    "Drought", "Drought", "Drought"
  ),
  levels = c("Water", "Drought")
)

coldata <- data.frame(
  row.names = colnames(count_matrix),
  condition = condition
)

# Create the DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = coldata,
  design = ~ condition
)

# Run the differential expression analysis
dds <- DESeq(dds)

# Compare drought samples with water controls
res <- results(
  dds,
  contrast = c("condition", "Drought", "Water")
)

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

# Recover the original gene IDs from featureCounts
if (all(grepl("^[0-9]+$", res_df$gene_id))) {
  res_df$gene_id <- counts$Geneid[
    as.integer(res_df$gene_id)
  ]
}

# Select significant DEGs using padj < 0.05 and |log2FC| >= 1
deg <- res_df[
  !is.na(res_df$padj) &
  res_df$padj < 0.05 &
  abs(res_df$log2FoldChange) >= 1,
]

up <- deg[
  deg$log2FoldChange >= 1,
]

down <- deg[
  deg$log2FoldChange <= -1,
]

# Save the differential expression results
write.csv(
  res_df,
  file.path(
    output_dir,
    "DESeq2_all_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  deg,
  file.path(
    output_dir,
    "DEGs_all_with_geneIDs.csv"
  ),
  row.names = FALSE
)

write.csv(
  up,
  file.path(
    output_dir,
    "DEGs_upregulated_with_geneIDs.csv"
  ),
  row.names = FALSE
)

write.csv(
  down,
  file.path(
    output_dir,
    "DEGs_downregulated_with_geneIDs.csv"
  ),
  row.names = FALSE
)

# Apply variance stabilizing transformation for visualization
vsd <- vst(
  dds,
  blind = FALSE
)

vst_matrix <- assay(vsd)

write.csv(
  vst_matrix,
  file.path(
    output_dir,
    "VST_expression_matrix.csv"
  )
)

# Generate the PCA plot
pca_data <- plotPCA(
  vsd,
  intgroup = "condition",
  returnData = TRUE
)

percent_var <- round(
  100 * attr(
    pca_data,
    "percentVar"
  )
)

p_pca <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    color = condition
  )
) +
  geom_point(size = 4) +
  labs(
    title = "PCA: Drought vs Water",
    x = paste0("PC1: ", percent_var[1], "%"),
    y = paste0("PC2: ", percent_var[2], "%")
  ) +
  theme_bw()

ggsave(
  file.path(
    output_dir,
    "PCA.png"
  ),
  p_pca,
  width = 8,
  height = 6,
  dpi = 300
)

# Classify genes for the volcano plot
res_df$significance <- "Not significant"

res_df$significance[
  !is.na(res_df$padj) &
  res_df$padj < 0.05 &
  res_df$log2FoldChange >= 1
] <- "Upregulated"

res_df$significance[
  !is.na(res_df$padj) &
  res_df$padj < 0.05 &
  res_df$log2FoldChange <= -1
] <- "Downregulated"

p_volcano <- ggplot(
  res_df,
  aes(
    x = log2FoldChange,
    y = -log10(padj),
    color = significance
  )
) +
  geom_point(alpha = 0.6) +
  theme_bw() +
  labs(
    title = "Volcano Plot",
    x = "Log2 Fold Change",
    y = "-Log10 adjusted p-value"
  )

ggsave(
  file.path(
    output_dir,
    "Volcano_Plot.png"
  ),
  p_volcano,
  width = 8,
  height = 6,
  dpi = 300
)
