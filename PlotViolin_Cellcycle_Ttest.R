# Clean up environment except for 'cds_sub_AcinaDucT_NewK_ReCluster'
rm(list = setdiff(ls(), "cds_sub_AcinaDucT_NewK_ReCluster"))

library(ggpubr) # Load ggpubr for advanced plotting

# Function to extract, plot, and transform data for a given gene and cell cycle phase
processGeneData <- function(cds, gene, phase = NULL, reCluster = NULL) {
  # Filter by gene
  geneData <- cds[rowData(cds)$gene_short_name %in% gene,]
  
  # Optionally filter by cell cycle phase
  if (!is.null(phase)) {
    geneData <- geneData[,geneData@colData@listData[["cell_cycle"]] %in% phase]
  }
  
  # Optionally filter by ReCluster
  if (!is.null(reCluster)) {
    geneData <- geneData[,geneData@colData@listData[["ReCluster"]] %in% reCluster]
  }
  
  # Plot gene expression violin plot
  plot_genes_violin(geneData, group_cells_by = "cell_cycle", ncol = 2) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Transform and label data for plotting
  transformedData <- t(as.data.frame(geneData@assays@data@listData[["counts"]]))
  transformedData <- cbind(transformedData, Type = ifelse(is.null(reCluster), "AcinaDucT", reCluster))
  colnames(transformedData)[2] <- "Type"
  
  transformedData$Value <- as.numeric(rownames(transformedData)) # Ensure values are numeric for plotting
  return(transformedData)
}

# Define genes and cell cycle phases of interest
genesOfInterest <- c("TOP2A")
cellCycles <- c("G1", "G2M", "S")
reClusters <- c(NULL, "CoreCD00")

# Process and plot data for each combination of gene, cell cycle, and reCluster
allData <- list()
for (gene in genesOfInterest) {
  for (phase in cellCycles) {
    for (reCluster in reClusters) {
      data <- processGeneData(cds_sub_AcinaDucT_NewK_ReCluster, gene, phase, reCluster)
      allData[[paste(gene, phase, ifelse(is.null(reCluster), "AcinaDucT", reCluster), sep = "_")]] <- data
    }
  }
}

# Function to plot combined data
plotCombinedData <- function(dataList) {
  for (dataName in names(dataList)) {
    data <- dataList[[dataName]]
    ggplot(data, aes(x = Type, y = Value, fill = Type)) + 
      geom_violin(trim = FALSE) +
      stat_summary(fun = mean, geom = "point", shape = 23, size = 2, color = "blue") +
      scale_fill_manual(values = c("AcinaDucT" = "#edd15f", "CoreCD00" = "#db993b")) +
      ggtitle(dataName) +
      theme_minimal() +
      theme(plot.title = element_text(size = 18, face = "bold"))
  }
}

# Plot all combined data
plotCombinedData(allData)



