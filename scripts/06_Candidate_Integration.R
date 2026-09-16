# 06_Candidate_Integration.R
# Integration of differential expression and functional evidence
# for the final candidate gene list

library(ggplot2)

input_file <- "D:/rnaseq_data/aligned/Final_Candidate_Genes.csv"
output_dir <- "D:/rnaseq_data/aligned"

# Read the final candidate table
candidate_data <- read.csv(
  input_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Keep the main evidence used for candidate prioritization
heatmap_data <- candidate_data[, c(
  "gene_id",
  "log2FoldChange",
  "Specific_GO_count",
  "Relevant_KEGG_count"
)]

colnames(heatmap_data) <- c(
  "gene_id",
  "log2FC",
  "GO_evidence",
  "KEGG_evidence"
)

# Scale the three evidence measures so they can be compared
heatmap_data$log2FC <- as.numeric(
  scale(heatmap_data$log2FC)
)

heatmap_data$GO_evidence <- as.numeric(
  scale(heatmap_data$GO_evidence)
)

heatmap_data$KEGG_evidence <- as.numeric(
  scale(heatmap_data$KEGG_evidence)
)

# Convert the data to long format for ggplot
plot_data <- rbind(
  data.frame(
    gene_id = heatmap_data$gene_id,
    Evidence = "log2FC",
    Value = heatmap_data$log2FC
  ),
  data.frame(
    gene_id = heatmap_data$gene_id,
    Evidence = "GO evidence",
    Value = heatmap_data$GO_evidence
  ),
  data.frame(
    gene_id = heatmap_data$gene_id,
    Evidence = "KEGG evidence",
    Value = heatmap_data$KEGG_evidence
  )
)

# Create the final candidate evidence heatmap
p <- ggplot(
  plot_data,
  aes(
    x = Evidence,
    y = gene_id,
    fill = Value
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.4
  ) +
  geom_text(
    aes(
      label = round(Value, 1)
    ),
    size = 3
  ) +
  scale_fill_gradient2(
    midpoint = 0
  ) +
  labs(
    title = "Final Candidate Evidence Heatmap",
    x = NULL,
    y = "Candidate Gene",
    fill = "Scaled value"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 10),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    panel.grid = element_blank()
  )

# Save the final figure
ggsave(
  file.path(
    output_dir,
    "Final_Candidate_Integration.png"
  ),
  plot = p,
  width = 9,
  height = 10,
  dpi = 300,
  bg = "white"
)
