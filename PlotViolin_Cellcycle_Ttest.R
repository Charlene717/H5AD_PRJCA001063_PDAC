# Clean up environment except for 'cds_sub_AcinaDucT_NewK_ReCluster'
rm(list = setdiff(ls(), "cds_sub_AcinaDucT_NewK_ReCluster"))

# Set up the environment
memory.limit(150000)
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

################################################################################
## Plot violin of module score
# Colors for Cell-Cycle
colors_cc <- c("#FF9912B3", "#2e6087", "#417034")

# Extract Markers and Cell Cycle information
markers <- c("PDAC", "EMT", "NPC", "NE", "ATR", "Migration", "Metastasis", "ACST")
cellCycleInfo <- as.character(cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["cell_cycle"]])
reClusterInfo <- as.character(cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["ReCluster"]])

# Define genes of interest
genesOfInterest <- c("PTK2", "TOP2A", "CGAS", "TP53")

# Function to process and plot gene data
processAndPlotGeneData <- function(gene) {
  geneData <- cds_sub_AcinaDucT_NewK_ReCluster[rowData(cds_sub_AcinaDucT_NewK_ReCluster)$gene_short_name %in% gene,]
  plot_genes_violin(geneData, group_cells_by = "cell_cycle", ncol = 2) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Convert counts to data frame
  countsDataFrame <- as.data.frame(geneData@assays@data@listData[["counts"]])
  return(countsDataFrame)
}

# Process and plot genes
geneDataFrames <- lapply(genesOfInterest, processAndPlotGeneData)

# Function to plot marker data
plotMarkerData <- function(markerData, markerName, reCluster, cellCycle, colors) {
  markerSum <- as.data.frame(cbind(markerData, reCluster))
  markerSum[,1] <- as.numeric(markerSum[,1])
  
  ggplot(markerSum, aes_string(x = "reCluster", y = names(markerSum)[1], fill = "reCluster")) +
    geom_violin() +
    geom_boxplot(alpha = 0.5, show.legend = FALSE) +
    geom_violin(trim = FALSE) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 2, color = "blue") +
    theme(axis.text.x = element_text(face = "bold", size = 10, angle = 45),
          axis.text.y = element_text(face = "bold"),
          axis.title.x = element_text(size = 14, face = "bold"),
          axis.title.y = element_text(size = 14, face = "bold")) +
    geom_hline(yintercept = mean(markerSum[,1]), linetype = 2) +
    scale_fill_manual(values = colors) +
    ggtitle(paste("Marker:", markerName))
  
  # Plot with cell cycle
  markerSumCC <- as.data.frame(cbind(markerData, reCluster, cellCycle))
  markerSumCC[,1] <- as.numeric(markerSumCC[,1])
  
  ggplot(markerSumCC, aes_string(x = "reCluster", y = names(markerSumCC)[1], fill = "cellCycle")) +
    geom_violin() +
    scale_fill_manual(values = colors) +
    theme(axis.text.x = element_text(face = "bold", size = 10, angle = 45),
          axis.text.y = element_text(face = "bold"),
          axis.title.x = element_text(size = 14, face = "bold"),
          axis.title.y = element_text(size = 14, face = "bold")) +
    ggtitle(paste("Marker:", markerName, "by Cell Cycle"))
}

# Example usage for a marker
plotMarkerData(cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["PDAC"]], "PDAC", reClusterInfo, cellCycleInfo, colors_cc)





