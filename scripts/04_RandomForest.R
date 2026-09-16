# 04_RandomForest.R
# Exploratory Random Forest analysis of VST expression data

library(randomForest)

input_file <- "D:/rnaseq_data/aligned/VST_expression_matrix.csv"
output_dir <- "D:/rnaseq_data/aligned"

# Read the variance-stabilized expression matrix
vst <- read.csv(
  input_file,
  row.names = 1,
  check.names = FALSE
)

# Transpose the matrix so samples are rows
vst <- as.data.frame(t(vst))

# Define the sample groups
condition <- factor(
  c(
    "Water", "Water", "Water",
    "Drought", "Drought", "Drought"
  )
)

# Train the Random Forest model
set.seed(123)

rf_model <- randomForest(
  x = vst,
  y = condition,
  ntree = 500,
  importance = TRUE
)

# Extract feature importance scores
importance_data <- as.data.frame(
  importance(rf_model)
)

importance_data$gene_id <- rownames(
  importance_data
)

importance_data <- importance_data[
  order(
    importance_data$MeanDecreaseAccuracy,
    decreasing = TRUE
  ),
]

# Save feature importance results
write.csv(
  importance_data,
  file.path(
    output_dir,
    "RF_feature_importance.csv"
  ),
  row.names = FALSE
)

# Save the top 20 ranked features
top20 <- head(
  importance_data,
  20
)

write.csv(
  top20,
  file.path(
    output_dir,
    "RF_top20_biomarkers.csv"
  ),
  row.names = FALSE
)
