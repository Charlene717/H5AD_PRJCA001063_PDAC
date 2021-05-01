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
head(x = pbmc_small[])

p1 <- DimPlot(pbmc_small, reduction = "pca", group.by = "CD_Features1")

pbmc_small <- RunUMAP(pbmc_small, reduction = "pca", dims = 1:30)
p2 <- DimPlot(pbmc_small, reduction = "umap", group.by = "CD_Features1")
p2
p1+p2


#將Seurat跑出的Cell cycle結果寫入Monocle3的cds檔

#cell_cycle <- marrow@active.ident
cds@colData@listData$CD_Features1 <- pbmc_small@meta.data[["CD_Features1"]]
plot_cells(cds, color_cells_by="CD_Features1", label_cell_groups=FALSE, show_trajectory_graph = FALSE)
