# 02_GO_Enrichment.R
# Visualization of significant GO enrichment results
# Rice drought vs water RNA-seq analysis

library(ggplot2)

go_dir <- "D:/rnaseq_data/aligned/GO_final_corrected"
output_dir <- "D:/rnaseq_data/aligned"

# Read significant Biological Process results for downregulated genes
bp_down <- read.csv(
  file.path(
    go_dir,
    "GO_BP_Downregulated_FDR.csv"
  ),
  check.names = FALSE
)

bp_down <- bp_down[
  !is.na(bp_down$padj) &
  bp_down$padj < 0.05,
]

bp_down <- bp_down[
  order(
    bp_down$Enrichment,
    decreasing = TRUE
  ),
]

# Plot enriched Biological Process terms
p_bp <- ggplot(
  bp_down,
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
    title = "GO Biological Process - Downregulated Genes",
    x = "Enrichment",
    y = "Biological Process"
  ) +
  theme_bw()

ggsave(
  file.path(
    output_dir,
    "GO_BP_Downregulated.png"
  ),
  p_bp,
  width = 10,
  height = 7,
  dpi = 300
)

# Read significant Molecular Function results for downregulated genes
mf_down <- read.csv(
  file.path(
    go_dir,
    "GO_MF_Downregulated_FDR.csv"
  ),
  check.names = FALSE
)

mf_down <- mf_down[
  !is.na(mf_down$padj) &
  mf_down$padj < 0.05,
]

mf_down <- mf_down[
  order(
    mf_down$Enrichment,
    decreasing = TRUE
  ),
]

# Keep the most enriched 20 terms so the plot remains readable
mf_down_plot <- head(
  mf_down,
  20
)

p_mf_down <- ggplot(
  mf_down_plot,
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
    title = "GO Molecular Function - Downregulated Genes",
    x = "Enrichment",
    y = "Molecular Function"
  ) +
  theme_bw()

ggsave(
  file.path(
    output_dir,
    "GO_MF_Downregulated.png"
  ),
  p_mf_down,
  width = 11,
  height = 8,
  dpi = 300
)

# Read significant Molecular Function results for upregulated genes
mf_up <- read.csv(
  file.path(
    go_dir,
    "GO_MF_Upregulated_FDR.csv"
  ),
  check.names = FALSE
)

mf_up <- mf_up[
  !is.na(mf_up$padj) &
  mf_up$padj < 0.05,
]

mf_up <- mf_up[
  order(
    mf_up$Enrichment,
    decreasing = TRUE
  ),
]

# Plot enriched Molecular Function terms
p_mf_up <- ggplot(
  mf_up,
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
    title = "GO Molecular Function - Upregulated Genes",
    x = "Enrichment",
    y = "Molecular Function"
  ) +
  theme_bw()

ggsave(
  file.path(
    output_dir,
    "GO_MF_Upregulated.png"
  ),
  p_mf_up,
  width = 10,
  height = 5,
  dpi = 300
)

# Read Cellular Component results for both directions
cc_up <- read.csv(
  file.path(
    go_dir,
    "GO_CC_Upregulated_FDR.csv"
  ),
  check.names = FALSE
)

cc_down <- read.csv(
  file.path(
    go_dir,
    "GO_CC_Downregulated_FDR.csv"
  ),
  check.names = FALSE
)

cc_up$Direction <- "Upregulated"
cc_down$Direction <- "Downregulated"

cc <- rbind(
  cc_up,
  cc_down
)

cc <- cc[
  !is.na(cc$padj) &
  cc$padj < 0.05,
]

# Combine significant Cellular Component terms into one plot
p_cc <- ggplot(
  cc,
  aes(
    x = Enrichment,
    y = reorder(
      Description,
      Enrichment
    ),
    fill = Direction
  )
) +
  geom_col() +
  labs(
    title = "GO Cellular Component",
    x = "Enrichment",
    y = "Cellular Component",
    fill = "Direction"
  ) +
  theme_bw()

ggsave(
  file.path(
    output_dir,
    "GO_CC.png"
  ),
  p_cc,
  width = 9,
  height = 5,
  dpi = 300
)
