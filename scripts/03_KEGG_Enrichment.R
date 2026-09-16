# 03_KEGG_Enrichment.R
# Visualization of significant KEGG enrichment results
# Rice drought vs water RNA-seq analysis

library(ggplot2)

kegg_dir <- "D:/rnaseq_data/aligned/GO_final_corrected/KEGG_final"
output_dir <- "D:/rnaseq_data/aligned"

# Read significant KEGG pathways for upregulated genes
kegg_up <- read.csv(
  file.path(
    kegg_dir,
    "KEGG_Upregulated_FDR.csv"
  ),
  check.names = FALSE
)

kegg_up <- kegg_up[
  !is.na(kegg_up$padj) &
  kegg_up$padj < 0.05,
]

kegg_up <- kegg_up[
  order(
    kegg_up$Enrichment,
    decreasing = TRUE
  ),
]

# Plot enriched pathways
p_up <- ggplot(
  kegg_up,
  aes(
    x = Enrichment,
    y = reorder(
      Description,
      Enrichment
    )
  )
) +
  geom_col() +
  labs(
    title = "KEGG Pathways - Upregulated Genes",
    x = "Enrichment",
    y = "KEGG Pathway"
  ) +
  theme_bw()

ggsave(
  file.path(
    output_dir,
    "KEGG_Upregulated.png"
  ),
  p_up,
  width = 10,
  height = 6,
  dpi = 300
)

# Read significant KEGG pathways for downregulated genes
kegg_down <- read.csv(
  file.path(
    kegg_dir,
    "KEGG_Downregulated_FDR.csv"
  ),
  check.names = FALSE
)

kegg_down <- kegg_down[
  !is.na(kegg_down$padj) &
  kegg_down$padj < 0.05,
]

kegg_down <- kegg_down[
  order(
    kegg_down$Enrichment,
    decreasing = TRUE
  ),
]

# Plot enriched pathways
p_down <- ggplot(
  kegg_down,
  aes(
    x = Enrichment,
    y = reorder(
      Description,
      Enrichment
    )
  )
) +
  geom_col() +
  labs(
    title = "KEGG Pathways - Downregulated Genes",
    x = "Enrichment",
    y = "KEGG Pathway"
  ) +
  theme_bw()

ggsave(
  file.path(
    output_dir,
    "KEGG_Downregulated.png"
  ),
  p_down,
  width = 11,
  height = 10,
  dpi = 300
)
