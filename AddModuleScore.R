# https://rdrr.io/github/satijalab/seurat/man/AddModuleScore.html
# https://academic.oup.com/nar/article/47/21/e133/5531181#185740462

cd_features <- list(c(
  'CD79B',
  'CD79A',
  'CD19',
  'CD180',
  'CD200',
  'CD3D',
  'CD2',
  'CD3E',
  'CD7',
  'CD8A',
  'CD14',
  'CD1C',
  'CD68',
  'CD9',
  'CD247'
))
pbmc_small <- AddModuleScore(
  object = marrow,
  features = cd_features,
  ctrl = 5,
  name = 'CD_Features'
)
# head(x = pbmc_small[])
# 
# p1 <- DimPlot(pbmc_small, reduction = "pca", group.by = "CD_Features1")
# 
# pbmc_small <- RunUMAP(pbmc_small, reduction = "pca", dims = 1:30)
# p2 <- DimPlot(pbmc_small, reduction = "umap", group.by = "CD_Features1")
# p2
# p1+p2


#將Seurat跑出的Cell cycle結果寫入Monocle3的cds檔

#cell_cycle <- marrow@active.ident
cds@colData@listData$CD_Features1 <- pbmc_small@meta.data[["CD_Features1"]]
plot_cells(cds, color_cells_by="CD_Features1", label_cell_groups=FALSE, show_trajectory_graph = FALSE)

##
PathName = setwd(getwd())
RVersion = "20210501V1"
dir.create(paste0(PathName,"/",RVersion))

#PDAC_Marker
Marker_file_Name <- c("GRUETZMANN_PANCREATIC_CANCER_UP")
PDAC_Marker_file <- paste0(PathName,"/",Marker_file_Name,".txt")

PDAC_Marker <- read.delim(PDAC_Marker_file, header = TRUE, stringsAsFactors = FALSE)
PDAC_Marker <- as.data.frame(PDAC_Marker[-1,])

marrow_PDAC_Marker <- AddModuleScore(
                      object = marrow,
                      features = PDAC_Marker,
                      ctrl = 5,
                      name = 'PDAC_Marker')
cds@colData@listData$PDAC_Marker <- marrow_PDAC_Marker@meta.data[["PDAC_Marker1"]]
plot_cells(cds, color_cells_by="PDAC_Marker", label_cell_groups=FALSE, show_trajectory_graph = FALSE)

### Metastasis
Meta_Marker_file_Name <- c("NAKAMURA_METASTASIS_MODEL_UP")
Meta_Marker_file <- paste0(PathName,"/",Meta_Marker_file_Name,".txt")

Meta_Marker <- read.delim(Meta_Marker_file, header = TRUE, stringsAsFactors = FALSE)
Meta_Marker <- as.data.frame(Meta_Marker[-1,])

marrow_Meta_Marker <- AddModuleScore(
                      object = marrow,
                      features = Meta_Marker,
                      ctrl = 5,
                      name = 'Meta_Marker')
cds@colData@listData$Meta_Marker <- marrow_Meta_Marker@meta.data[["Meta_Marker1"]]
plot_cells(cds, color_cells_by="Meta_Marker", label_cell_groups=FALSE, show_trajectory_graph = FALSE)

### Migration
Mig_Marker_file_Name <- c("REACTOME_SEMA4D_MEDIATED_INHIBITION_OF_CELL_ATTACHMENT_AND_MIGRATION")
Mig_Marker_file <- paste0(PathName,"/",Mig_Marker_file_Name,".txt")

Mig_Marker <- read.delim(Mig_Marker_file, header = TRUE, stringsAsFactors = FALSE)
Mig_Marker <- as.data.frame(Mig_Marker[-1,])

marrow_Mig_Marker <- AddModuleScore(
  object = marrow,
  features = Mig_Marker,
  ctrl = 5,
  name = 'Mig_Marker')
cds@colData@listData$Mig_Marker <- marrow_Mig_Marker@meta.data[["Mig_Marker1"]]
plot_cells(cds, color_cells_by="Mig_Marker", label_cell_groups=FALSE, show_trajectory_graph = FALSE)

### EMT
EMT_Marker_file_Name <- c("ALONSO_METASTASIS_EMT_UP")
EMT_Marker_file <- paste0(PathName,"/",EMT_Marker_file_Name,".txt")

EMT_Marker <- read.delim(EMT_Marker_file, header = TRUE, stringsAsFactors = FALSE)
EMT_Marker <- as.data.frame(EMT_Marker[-1,])

marrow_EMT_Marker <- AddModuleScore(
  object = marrow,
  features = EMT_Marker,
  ctrl = 5,
  name = 'EMT_Marker')
cds@colData@listData$EMT_Marker <- marrow_EMT_Marker@meta.data[["EMT_Marker1"]]
plot_cells(cds, color_cells_by="EMT_Marker", label_cell_groups=FALSE, show_trajectory_graph = FALSE)
######################################  cds_subset_PDAC ########################################
########################  PDAC ##########################
cds_sub_PDAC <- choose_cells(cds)
#cds_subset <- reduce_dimension(cds_subset)
plot_cells(cds_sub_PDAC, label_cell_groups=FALSE, show_trajectory_graph = FALSE)
plot_cells(cds_sub_PDAC, color_cells_by="PDAC_Marker", label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1.5)
plot_cells(cds_sub_PDAC, color_cells_by="Meta_Marker", label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1.5)
plot_cells(cds_sub_PDAC, color_cells_by="Mig_Marker", label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1.5)

plot_cells(cds_sub_PDAC, genes ="TOP2A", label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1.5)
plot_cells(cds_sub_PDAC, genes ="TOP2B", label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1.5)


############ cds_subset_PDAC to  Seurat
library(Seurat)
# https://github.com/rstudio/rstudio/issues/4741
library(SummarizedExperiment) 
library(AnnotationDbi)
library(org.Mm.eg.db)
library('org.Hs.eg.db')
library(Hmisc)

CCdata_cds_sub_PDAC <- cds_sub_PDAC@assays@data@listData$counts
colnames(CCdata_cds_sub_PDAC) = cds_sub_PDAC@assays@data@listData[["counts"]]@Dimnames[[2]]
rownames(CCdata_cds_sub_PDAC) = cds_sub_PDAC@assays@data@listData[["counts"]]@Dimnames[[1]]


DataCellcycle_cds_sub_PDAC <- CCdata_cds_sub_PDAC
#rownames(DataCellcycle) <- make.names(DataCellcycle[,1], unique = TRUE)
#DataCellcycle <- DataCellcycle[1:length(data[,1]), 2:length(data[1,])]
DataCellcycle_cds_sub_PDAC <- as(as.matrix(DataCellcycle_cds_sub_PDAC), "dgCMatrix")

marrow_cds_sub_PDAC <- CreateSeuratObject(counts = DataCellcycle_cds_sub_PDAC)
marrow_cds_sub_PDAC <- NormalizeData(marrow_cds_sub_PDAC)
marrow_cds_sub_PDAC <- FindVariableFeatures(marrow_cds_sub_PDAC, selection.method = "vst")
marrow_cds_sub_PDAC <- ScaleData(marrow_cds_sub_PDAC, features = rownames(marrow_cds_sub_PDAC))

#PDAC_Marker
Marker_file_Name <- c("GRUETZMANN_PANCREATIC_CANCER_UP")
PDAC_Marker_file <- paste0(PathName,"/",Marker_file_Name,".txt")

PDAC_Marker <- read.delim(PDAC_Marker_file, header = TRUE, stringsAsFactors = FALSE)
PDAC_Marker <- as.data.frame(PDAC_Marker[-1,])

marrow_PDAC_Marker <- AddModuleScore(
  object = marrow_cds_sub_PDAC,
  features = PDAC_Marker,
  ctrl = 5,
  name = 'PDAC_Marker')
cds_sub_PDAC@colData@listData$PDAC_Marker <- marrow_PDAC_Marker@meta.data[["PDAC_Marker1"]]
plot_cells(cds_sub_PDAC, color_cells_by="PDAC_Marker", label_cell_groups=FALSE,
           show_trajectory_graph = FALSE,cell_size = 1)

### Metastasis
Meta_Marker_file_Name <- c("NAKAMURA_METASTASIS_MODEL_UP")
Meta_Marker_file <- paste0(PathName,"/",Meta_Marker_file_Name,".txt")

Meta_Marker <- read.delim(Meta_Marker_file, header = TRUE, stringsAsFactors = FALSE)
Meta_Marker <- as.data.frame(Meta_Marker[-1,])

marrow_Meta_Marker <- AddModuleScore(
  object = marrow_cds_sub_PDAC,
  features = Meta_Marker,
  ctrl = 5,
  name = 'Meta_Marker')
cds_sub_PDAC@colData@listData$Meta_Marker <- marrow_Meta_Marker@meta.data[["Meta_Marker1"]]
plot_cells(cds_sub_PDAC, color_cells_by="Meta_Marker", label_cell_groups=FALSE, 
           show_trajectory_graph = FALSE,cell_size = 1)


### Migration
Mig_Marker_file_Name <- c("REACTOME_SEMA4D_MEDIATED_INHIBITION_OF_CELL_ATTACHMENT_AND_MIGRATION")
Mig_Marker_file <- paste0(PathName,"/",Mig_Marker_file_Name,".txt")

Mig_Marker <- read.delim(Mig_Marker_file, header = TRUE, stringsAsFactors = FALSE)
Mig_Marker <- as.data.frame(Mig_Marker[-1,])

marrow_Mig_Marker <- AddModuleScore(
  object = marrow_cds_sub_PDAC,
  features = Mig_Marker,
  ctrl = 5,
  name = 'Mig_Marker')
cds_sub_PDAC@colData@listData$Mig_Marker <- marrow_Mig_Marker@meta.data[["Mig_Marker1"]]
plot_cells(cds_sub_PDAC, color_cells_by="Mig_Marker", label_cell_groups=FALSE,  
           show_trajectory_graph = FALSE,cell_size = 1)

### EMT
EMT_Marker_file_Name <- c("ALONSO_METASTASIS_EMT_UP")
EMT_Marker_file <- paste0(PathName,"/",EMT_Marker_file_Name,".txt")

EMT_Marker <- read.delim(EMT_Marker_file, header = TRUE, stringsAsFactors = FALSE)
EMT_Marker <- as.data.frame(EMT_Marker[-1,])

marrow_EMT_Marker <- AddModuleScore(
  object = marrow_cds_sub_PDAC,
  features = EMT_Marker,
  ctrl = 5,
  name = 'EMT_Marker')
cds_sub_PDAC@colData@listData$EMT_Marker <- marrow_EMT_Marker@meta.data[["EMT_Marker1"]]
plot_cells(cds_sub_PDAC, color_cells_by="EMT_Marker", label_cell_groups=FALSE, 
           show_trajectory_graph = FALSE,cell_size = 1)
