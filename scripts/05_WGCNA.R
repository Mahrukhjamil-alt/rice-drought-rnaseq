# 05_WGCNA.R
# Weighted gene co-expression network analysis
# Rice drought vs water RNA-seq data

library(WGCNA)

options(stringsAsFactors = FALSE)

input_file <- "D:/rnaseq_data/aligned/VST_expression_matrix.csv"
output_dir <- "D:/rnaseq_data/aligned/WGCNA_results"

dir.create(
  output_dir,
  showWarnings = FALSE
)

# Read the VST expression matrix
vst_matrix <- read.csv(
  input_file,
  row.names = 1,
  check.names = FALSE
)

# Select the 2,000 most variable genes
gene_variance <- apply(
  vst_matrix,
  1,
  var
)

top_genes <- names(
  sort(
    gene_variance,
    decreasing = TRUE
  )
)[1:2000]

# WGCNA expects samples as rows and genes as columns
datExpr_wgcna <- as.data.frame(
  t(
    vst_matrix[
      top_genes,
      drop = FALSE
    ]
  )
)

# Check for problematic samples or genes
gsg <- goodSamplesGenes(
  datExpr_wgcna,
  verbose = 3
)

datExpr_wgcna <- datExpr_wgcna[
  gsg$goodSamples,
  gsg$goodGenes
]

# Define the drought trait
trait <- data.frame(
  Drought = c(
    0, 0, 0,
    1, 1, 1
  )
)

rownames(trait) <- rownames(
  datExpr_wgcna
)

# Build the co-expression network
net_wgcna <- blockwiseModules(
  datExpr_wgcna,
  power = 12,
  networkType = "unsigned",
  TOMType = "unsigned",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 3
)

# Use the module colors returned directly by WGCNA
module_colors <- net_wgcna$colors

# Calculate module-trait relationships
MEs <- orderMEs(
  net_wgcna$MEs
)

moduleTraitCor <- cor(
  MEs,
  trait,
  use = "p"
)

moduleTraitPvalue <- corPvalueStudent(
  moduleTraitCor,
  nSamples = nrow(datExpr_wgcna)
)

# Save module-trait correlations
write.csv(
  moduleTraitCor,
  file.path(
    output_dir,
    "Module_Trait_Correlation.csv"
  )
)

write.csv(
  moduleTraitPvalue,
  file.path(
    output_dir,
    "Module_Trait_Pvalues.csv"
  )
)

# Summarize the number of genes in each module
module_summary <- as.data.frame(
  table(module_colors)
)

colnames(module_summary) <- c(
  "Module",
  "Gene_Count"
)

write.csv(
  module_summary,
  file.path(
    output_dir,
    "WGCNA_Module_Summary.csv"
  ),
  row.names = FALSE
)

# Extract genes from the drought-associated turquoise module
drought_module_color <- "turquoise"

drought_genes <- names(
  module_colors
)[
  module_colors ==
    drought_module_color
]

drought_gene_table <- data.frame(
  gene_id = drought_genes
)

write.csv(
  drought_gene_table,
  file.path(
    output_dir,
    "WGCNA_Drought_Associated_Genes.csv"
  ),
  row.names = FALSE
)
