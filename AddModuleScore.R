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

Marker_file_Name <- c("GRUETZMANN_PANCREATIC_CANCER_UP")
PDAC_Marker_file <- paste0(PathName,"/",Marker_file_Name,".txt")

PDAC_Marker <- read.delim(PDAC_Marker_file, header = TRUE, stringsAsFactors = FALSE)
PDAC_Marker <- as.data.frame(PDAC_Marker[-1,])

marrow_PDAC_Marker <- AddModuleScore(
                      object = marrow,
                      features = PDAC_Marker,
                      ctrl = 5,
                      name = 'CD_Features')
cds@colData@listData$PDAC_Marker <- marrow_PDAC_Marker@meta.data[["CD_Features1"]]
plot_cells(cds, color_cells_by="PDAC_Marker", label_cell_groups=FALSE, show_trajectory_graph = FALSE)

