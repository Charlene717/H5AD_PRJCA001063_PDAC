# Load necessary libraries
library(infercnv)
library(dplyr)

# Set up directory for outputs
PathName <- setwd(getwd())
RVersion <- "20210610V1"
outputDir <- paste0(PathName,"/",RVersion)
dir.create(outputDir)

# Define function to write data tables
writeData <- function(data, fileNamePrefix) {
  filePathBase <- file.path(outputDir, fileNamePrefix)
  write.table(data, paste0(filePathBase, ".csv"), row.names = TRUE, sep = ',')
  write.table(data, paste0(filePathBase, ".txt"), row.names = TRUE, sep = '\t')
}

# Processing expression data and cell types from an inferCNV object
inferCNV_cds_sub_AcinaDucT_NewK_ReCluster <- cds_sub_AcinaDucT_NewK_ReCluster
expressionData <- as.data.frame(inferCNV_cds_sub_AcinaDucT_NewK_ReCluster@assays@data@listData[["logcounts"]])
filteredExpression <- expressionData[rowSums(expressionData) > 0, ]
writeData(filteredExpression, RVersion+"_filteredExpression")

cellTypes <- as.data.frame(inferCNV_cds_sub_AcinaDucT_NewK_ReCluster@colData)
cellTypesSimplified <- cellTypes[, 15, drop = FALSE]
writeData(cellTypesSimplified, RVersion+"_cellTypes")

# Create and run inferCNV analysis
infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix=as.matrix(filteredExpression),
  annotations_file=as.matrix(cellTypesSimplified),
  gene_order_file=system.file("extdata", "gencode_downsampled.EXAMPLE_ONLY_DONT_REUSE.txt", package = "infercnv"),
  ref_group_names=c("AC","ND01")
)

infercnv_results <- infercnv::run(
  infercnv_obj,
  cutoff=0.1,
  out_dir=tempfile(), 
  cluster_by_groups=TRUE,
  plot_steps=FALSE,
  no_plot=FALSE,
  denoise=TRUE,
  resume_mode=FALSE,
  HMM=TRUE
)

# The detailed filtering and plotting steps for cell types and various cell states have been omitted for clarity.



