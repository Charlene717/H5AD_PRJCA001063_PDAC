##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(300000)

##### Loading the data ######
  load("D:/Dropbox/##_GitHub/##_PHH_Lab/#_H5AD_PRJCA001063_PDAC/2022-07-14_Com_PDAC/scRNA.SeuObj_CDS_PRJCA001063_Combine_Anno_ReDR.RData")
  
##### Load Packages #####
  if(!require("Seurat")) install.packages("Seurat")
  if(!require("tidyverse")) install.packages("tidyverse")
  library(tidyverse)
  library(Seurat)


##### Plot #####
  FeaturePlot(scRNA.SeuObj, features = c("TOP2A")) %>% BeautifyggPlot(.,LegPos = c(0.1, 0.2))
  FeaturePlot(scRNA.SeuObj, features = c("MARCO")) %>% BeautifyggPlot(.,LegPos = c(0.1, 0.2))

  DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "Cell_type")  %>% BeautifyggPlot(.,LegPos = c(0.05, 0.15))
  DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "DataSetID")  %>% BeautifyggPlot(.,LegPos = c(0.05, 0.15))
  DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "CONDITION")  %>% BeautifyggPlot(.,LegPos = c(0.05, 0.15))
  DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "cell_cycle")  %>% BeautifyggPlot(.,LegPos = c(0.05, 0.15))
  DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "seurat_clusters")  %>% BeautifyggPlot(.,LegPos = c(1.05, 0.5))
  
##### Recluster #####
  # Seurat re-clustering a cell subset but cell identity numbers are not completely showing up
  # https://www.biostars.org/p/9485834/#9486162
  DefaultAssay(scRNA.SeuObj) <- "integrated"
  scRNA.SeuObj <- FindClusters(scRNA.SeuObj, resolution = 2)

  ## Plot
  DefaultAssay(scRNA.SeuObj) <- "RNA"
  DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "seurat_clusters")  %>% BeautifyggPlot(.,LegPos = c(1.05, 0.5))
  FeaturePlot(scRNA.SeuObj, features = c("MARCO")) %>% BeautifyggPlot(.,LegPos = c(0.1, 0.2))
  

##### Extract macrophage 2 ReDR & Recluster #####
  scRNA_Mac.SeuObj <- scRNA.SeuObj[,scRNA.SeuObj@meta.data[["Cell_type"]] %in% "Macrophage cell"]
  # # subset Seurat to only Macrophage cell
  # scRNA_Mac.SeuObj <- subset(x = ,scRNA.SeuObj, idents = "Macrophage cell")
  
  DefaultAssay(scRNA_Mac.SeuObj) <- "integrated"
  scRNA_Mac.SeuObj <- FindVariableFeatures(scRNA_Mac.SeuObj)
  scRNA_Mac.SeuObj <- ScaleData(scRNA_Mac.SeuObj, verbose = FALSE)
  scRNA_Mac.SeuObj <- RunPCA(scRNA_Mac.SeuObj, npcs = 300, verbose = FALSE)
  # scRNA_Mac.SeuObj <- RunUMAP(scRNA_Mac.SeuObj, reduction = "pca", dims = 1:160,n.neighbors = 20,min.dist = 0.3)
  scRNA_Mac.SeuObj <- RunUMAP(scRNA_Mac.SeuObj, reduction = "pca", dims = 1:300,n.neighbors = 30,min.dist = 0.3)
  
  scRNA_Mac.SeuObj <- FindNeighbors(scRNA_Mac.SeuObj, reduction = "pca", dims = 1:300)
  scRNA_Mac.SeuObj <- FindClusters(scRNA_Mac.SeuObj, resolution = 0.2)
  DimPlot(scRNA_Mac.SeuObj, reduction = "umap",group.by = "seurat_clusters")  %>% BeautifyggPlot(.,LegPos = c(1.05, 0.5))
  FeaturePlot(scRNA_Mac.SeuObj, features = c("MARCO")) %>% BeautifyggPlot(.,LegPos = c(0.05, 0.2))
  FeaturePlot(scRNA_Mac.SeuObj, features = c("IL1B")) %>% BeautifyggPlot(.,LegPos = c(0.05, 0.2))
  
  scRNA_Mac.SeuObj <- RenameIdents(scRNA_Mac.SeuObj, `0` = "MarcoP", `1` = "MarcoN", `2` = "MarcoN",
                               `3` = "MarcoP", `4` = "MarcoP", `5` = "MarcoN", `6` = "MarcoN")
  DimPlot(scRNA_Mac.SeuObj, reduction = "umap")  %>% BeautifyggPlot(.,LegPos = c(0.05, 0.9))
  
##### Volcano plot #####  
  scRNA_Mac.SeuObj$MarcoType <- Idents(scRNA_Mac.SeuObj)
  #### Define group by different phenotype ####
  source("FUN_Find_Markers.R")
  dir.create(paste0("FindMarkers"))
  PDACMarker.lt <- Find_Markers(scRNA_Mac.SeuObj, ident1="MarcoP", ident2="MarcoN", CellType ="Macro",Path = getwd(),
               log2FC=0.5, Pval=0.05,ResultFolder = "FindMarkers",ProjectTitle="Tar")
  
  source("FUN_VolcanoPlot.R")
  VolcanoPlot(PDACMarker.lt[["TarMarker.S"]],
              PDACMarker.lt[["TarMarker.S_Pos_List"]],
              PDACMarker.lt[["TarMarker.S_Neg_List"]], 
              log2FC = 0.5,PValue = 0.05,
              ShowGeneNum = 10
              ) + 
              ggtitle(paste0("MarcoP vs. MarcoN")) -> Plot.Volcano
  Plot.Volcano
  
  tiff(file = paste0("FindMarkers/","Marco.tif"), width = 17, height = 17, units = "cm", res = 200)
    print(Plot.Volcano)
  graphics.off()


  
