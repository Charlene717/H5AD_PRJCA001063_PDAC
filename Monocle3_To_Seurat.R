# Load necessary libraries
library(Seurat)
library(SummarizedExperiment) 
library(AnnotationDbi)
library(org.Mm.eg.db)
library('org.Hs.eg.db')
library(Hmisc)
library(monocle3)

# Convert Monocle3 Object to Seurat Object function
Monocle3_To_Seurat <- function(cds, FigureName, outputPath) {
  # Extract and prepare count data
  CCdata <- cds@assays@data@listData$counts
  colnames(CCdata) <- colnames(CCdata)
  rownames(CCdata) <- rownames(CCdata)
  DataCellcycle <- as(as.matrix(CCdata), "dgCMatrix")
  
  # Create Seurat object and preprocess data
  marrow <- CreateSeuratObject(counts = DataCellcycle)
  marrow <- NormalizeData(marrow)
  marrow <- FindVariableFeatures(marrow, selection.method = "vst")
  marrow <- ScaleData(marrow, features = rownames(marrow))
  
  # Ensure consistency of row names
  rowDataNames <- rownames(CCdata)
  marrow@assays[["RNA"]]@counts@Dimnames[[1]] <- rowDataNames
  marrow@assays[["RNA"]]@data@Dimnames[[1]] <- rowDataNames
  marrow@commands[["ScaleData.RNA"]]@params[["features"]] <- rowDataNames
  
  # Run PCA and plot Dimensional Heatmaps
  marrow <- RunPCA(marrow, features = VariableFeatures(marrow), ndims.print = 1:10, nfeatures.print = 10)
  
  # Plotting Dimensional Heatmaps for each specified dimension pair and save
  for (dims in list(c(1, 2), c(3, 4), c(5, 6), c(7, 8), c(9, 10))) {
    dimPairName <- paste(dims, collapse = "&")
    plotPath <- paste0(outputPath, "/", RVersion, "_", FigureName, "_CC_DimHP", dimPairName, ".png")
    png(plotPath)
    DimHeatmap(marrow, dims = dims)
    dev.off()
  }
  
  return(marrow)
}

# Example usage:
# marrow <- Monocle3_To_Seurat(cds = myCDS, FigureName = "MyFigure", outputPath = "path/to/output")



