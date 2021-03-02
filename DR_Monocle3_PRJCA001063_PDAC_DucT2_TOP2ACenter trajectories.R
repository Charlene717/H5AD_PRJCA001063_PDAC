########################  DucT2_TOP2ACenter trajectories_T2 ##########################
cds_sub_DucT2_TOP2ACenter_T2 <- choose_graph_segments(cds_sub_DucT2 ,clear_cds = FALSE)

plot_cells(cds_sub_DucT2_TOP2ACenter_T2, color_cells_by="cell_cycle",cell_size=2, label_cell_groups=FALSE) + scale_color_manual(values = colorsT)


CCdata_sub_DucT2_TOP2ACenter_T2 <- cds_sub_DucT2_TOP2ACenter_T2@assays@data@listData$counts
colnames(CCdata_sub_DucT2_TOP2ACenter_T2) = cds_sub_DucT2_TOP2ACenter_T2@assays@data@listData[["counts"]]@Dimnames[[2]]
rownames(CCdata_sub_DucT2_TOP2ACenter_T2) = cds_sub_DucT2_TOP2ACenter_T2@assays@data@listData[["counts"]]@Dimnames[[1]]


DataCellcycle_sub_DucT2_TOP2ACenter_T2 <- CCdata_sub_DucT2_TOP2ACenter_T2
DataCellcycle_sub_DucT2_TOP2ACenter_T2 <- as(as.matrix(DataCellcycle_sub_DucT2_TOP2ACenter_T2), "dgCMatrix")

marrow_sub_DucT2_TOP2ACenter_T2 <- CreateSeuratObject(counts = DataCellcycle_sub_DucT2_TOP2ACenter_T2)
marrow_sub_DucT2_TOP2ACenter_T2 <- NormalizeData(marrow_sub_DucT2_TOP2ACenter_T2)
marrow_sub_DucT2_TOP2ACenter_T2 <- FindVariableFeatures(marrow_sub_DucT2_TOP2ACenter_T2, selection.method = "vst")
marrow_sub_DucT2_TOP2ACenter_T2 <- ScaleData(marrow_sub_DucT2_TOP2ACenter_T2, features = rownames(marrow_sub_DucT2_TOP2ACenter_T2))
#錯誤: 無法配置大小為 3.9 Gb 的向量
#https://d.cosx.org/d/413001-413001
#memory.limit(15000)


marrow_sub_DucT2_TOP2ACenter_T2 <- RunPCA(marrow_sub_DucT2_TOP2ACenter_T2, features = VariableFeatures(marrow_sub_DucT2_TOP2ACenter_T2), ndims.print = 1:10, nfeatures.print = 25)

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T2_DimHeatmap1.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T2, dims = c(1,2),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T2_DimHeatmap2.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T2, dims = c(3,4),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T2_DimHeatmap3.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T2, dims = c(5,6),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T2_DimHeatmap4.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T2, dims = c(7,8),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T2_DimHeatmap5.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T2, dims = c(9,10),nfeatures = 30)
dev.off() # 關閉輸出圖檔

DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T2, dims = c(1:2),nfeatures = 50)

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T2_DimHeatmapMa.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T2, dims = c(1:18))
dev.off() # 關閉輸出圖檔



########################  DucT2_TOP2ACenter trajectories_T3 ##########################
cds_sub_DucT2_TOP2ACenter_T3 <- choose_graph_segments(cds_sub_DucT2 ,clear_cds = FALSE)

plot_cells(cds_sub_DucT2_TOP2ACenter_T3, color_cells_by="cell_cycle",cell_size=2, label_cell_groups=FALSE) + scale_color_manual(values = colorsT)

CCdata_sub_DucT2_TOP2ACenter_T3 <- cds_sub_DucT2_TOP2ACenter_T3@assays@data@listData$counts
colnames(CCdata_sub_DucT2_TOP2ACenter_T3) = cds_sub_DucT2_TOP2ACenter_T3@assays@data@listData[["counts"]]@Dimnames[[2]]
rownames(CCdata_sub_DucT2_TOP2ACenter_T3) = cds_sub_DucT2_TOP2ACenter_T3@assays@data@listData[["counts"]]@Dimnames[[1]]


DataCellcycle_sub_DucT2_TOP2ACenter_T3 <- CCdata_sub_DucT2_TOP2ACenter_T3
DataCellcycle_sub_DucT2_TOP2ACenter_T3 <- as(as.matrix(DataCellcycle_sub_DucT2_TOP2ACenter_T3), "dgCMatrix")

marrow_sub_DucT2_TOP2ACenter_T3 <- CreateSeuratObject(counts = DataCellcycle_sub_DucT2_TOP2ACenter_T3)
marrow_sub_DucT2_TOP2ACenter_T3 <- NormalizeData(marrow_sub_DucT2_TOP2ACenter_T3)
marrow_sub_DucT2_TOP2ACenter_T3 <- FindVariableFeatures(marrow_sub_DucT2_TOP2ACenter_T3, selection.method = "vst")
marrow_sub_DucT2_TOP2ACenter_T3 <- ScaleData(marrow_sub_DucT2_TOP2ACenter_T3, features = rownames(marrow_sub_DucT2_TOP2ACenter_T3))
#錯誤: 無法配置大小為 3.9 Gb 的向量
#https://d.cosx.org/d/413001-413001
#memory.limit(15000)


marrow_sub_DucT2_TOP2ACenter_T3 <- RunPCA(marrow_sub_DucT2_TOP2ACenter_T3, features = VariableFeatures(marrow_sub_DucT2_TOP2ACenter_T3), ndims.print = 1:10, nfeatures.print = 25)

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T3_DimHeatmap1.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T3, dims = c(1,2),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T3_DimHeatmap2.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T3, dims = c(3,4),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T3_DimHeatmap3.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T3, dims = c(5,6),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T3_DimHeatmap4.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T3, dims = c(7,8),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T3_DimHeatmap5.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T3, dims = c(9,10),nfeatures = 30)
dev.off() # 關閉輸出圖檔

DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T3, dims = c(1:2),nfeatures = 50)

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T3_DimHeatmapMa.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T3, dims = c(1:18))
dev.off() # 關閉輸出圖檔


########################  DucT2_TOP2ACenter trajectories_T4 ##########################
cds_sub_DucT2_TOP2ACenter_T4 <- choose_graph_segments(cds_sub_DucT2 ,clear_cds = FALSE)

plot_cells(cds_sub_DucT2_TOP2ACenter_T4, color_cells_by="cell_cycle",cell_size=2, label_cell_groups=FALSE) + scale_color_manual(values = colorsT)

CCdata_sub_DucT2_TOP2ACenter_T4 <- cds_sub_DucT2_TOP2ACenter_T4@assays@data@listData$counts
colnames(CCdata_sub_DucT2_TOP2ACenter_T4) = cds_sub_DucT2_TOP2ACenter_T4@assays@data@listData[["counts"]]@Dimnames[[2]]
rownames(CCdata_sub_DucT2_TOP2ACenter_T4) = cds_sub_DucT2_TOP2ACenter_T4@assays@data@listData[["counts"]]@Dimnames[[1]]


DataCellcycle_sub_DucT2_TOP2ACenter_T4 <- CCdata_sub_DucT2_TOP2ACenter_T4
DataCellcycle_sub_DucT2_TOP2ACenter_T4 <- as(as.matrix(DataCellcycle_sub_DucT2_TOP2ACenter_T4), "dgCMatrix")

marrow_sub_DucT2_TOP2ACenter_T4 <- CreateSeuratObject(counts = DataCellcycle_sub_DucT2_TOP2ACenter_T4)
marrow_sub_DucT2_TOP2ACenter_T4 <- NormalizeData(marrow_sub_DucT2_TOP2ACenter_T4)
marrow_sub_DucT2_TOP2ACenter_T4 <- FindVariableFeatures(marrow_sub_DucT2_TOP2ACenter_T4, selection.method = "vst")
marrow_sub_DucT2_TOP2ACenter_T4 <- ScaleData(marrow_sub_DucT2_TOP2ACenter_T4, features = rownames(marrow_sub_DucT2_TOP2ACenter_T4))
#錯誤: 無法配置大小為 3.9 Gb 的向量
#https://d.cosx.org/d/413001-413001
#memory.limit(15000)


marrow_sub_DucT2_TOP2ACenter_T4 <- RunPCA(marrow_sub_DucT2_TOP2ACenter_T4, features = VariableFeatures(marrow_sub_DucT2_TOP2ACenter_T4), ndims.print = 1:10, nfeatures.print = 25)

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T4_DimHeatmap1.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T4, dims = c(1,2),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T4_DimHeatmap2.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T4, dims = c(3,4),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T4_DimHeatmap3.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T4, dims = c(5,6),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T4_DimHeatmap4.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T4, dims = c(7,8),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T4_DimHeatmap5.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T4, dims = c(9,10),nfeatures = 30)
dev.off() # 關閉輸出圖檔

DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T4, dims = c(1:2),nfeatures = 50)

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T4_DimHeatmapMa.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T4, dims = c(1:18))
dev.off() # 關閉輸出圖檔


########################  DucT2_TOP2ACenter trajectories_T5 ##########################
cds_sub_DucT2_TOP2ACenter_T5 <- choose_graph_segments(cds_sub_DucT2 ,clear_cds = FALSE)

plot_cells(cds_sub_DucT2_TOP2ACenter_T5, color_cells_by="cell_cycle",cell_size=2, label_cell_groups=FALSE) + scale_color_manual(values = colorsT)

CCdata_sub_DucT2_TOP2ACenter_T5 <- cds_sub_DucT2_TOP2ACenter_T5@assays@data@listData$counts
colnames(CCdata_sub_DucT2_TOP2ACenter_T5) = cds_sub_DucT2_TOP2ACenter_T5@assays@data@listData[["counts"]]@Dimnames[[2]]
rownames(CCdata_sub_DucT2_TOP2ACenter_T5) = cds_sub_DucT2_TOP2ACenter_T5@assays@data@listData[["counts"]]@Dimnames[[1]]


DataCellcycle_sub_DucT2_TOP2ACenter_T5 <- CCdata_sub_DucT2_TOP2ACenter_T5
DataCellcycle_sub_DucT2_TOP2ACenter_T5 <- as(as.matrix(DataCellcycle_sub_DucT2_TOP2ACenter_T5), "dgCMatrix")

marrow_sub_DucT2_TOP2ACenter_T5 <- CreateSeuratObject(counts = DataCellcycle_sub_DucT2_TOP2ACenter_T5)
marrow_sub_DucT2_TOP2ACenter_T5 <- NormalizeData(marrow_sub_DucT2_TOP2ACenter_T5)
marrow_sub_DucT2_TOP2ACenter_T5 <- FindVariableFeatures(marrow_sub_DucT2_TOP2ACenter_T5, selection.method = "vst")
marrow_sub_DucT2_TOP2ACenter_T5 <- ScaleData(marrow_sub_DucT2_TOP2ACenter_T5, features = rownames(marrow_sub_DucT2_TOP2ACenter_T5))
#錯誤: 無法配置大小為 3.9 Gb 的向量
#https://d.cosx.org/d/413001-413001
#memory.limit(15000)


marrow_sub_DucT2_TOP2ACenter_T5 <- RunPCA(marrow_sub_DucT2_TOP2ACenter_T5, features = VariableFeatures(marrow_sub_DucT2_TOP2ACenter_T5), ndims.print = 1:10, nfeatures.print = 25)

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T5_DimHeatmap1.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T5, dims = c(1,2),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T5_DimHeatmap2.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T5, dims = c(3,4),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T5_DimHeatmap3.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T5, dims = c(5,6),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T5_DimHeatmap4.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T5, dims = c(7,8),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T5_DimHeatmap5.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T5, dims = c(9,10),nfeatures = 30)
dev.off() # 關閉輸出圖檔

DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T5, dims = c(1:2),nfeatures = 50)

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T5_DimHeatmapMa.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T5, dims = c(1:18))
dev.off() # 關閉輸出圖檔


########################  DucT2_TOP2ACenter trajectories_T6 ##########################
cds_sub_DucT2_TOP2ACenter_T6 <- choose_graph_segments(cds_sub_DucT2 ,clear_cds = FALSE)

plot_cells(cds_sub_DucT2_TOP2ACenter_T6, color_cells_by="cell_cycle",cell_size=2, label_cell_groups=FALSE) + scale_color_manual(values = colorsT)

CCdata_sub_DucT2_TOP2ACenter_T6 <- cds_sub_DucT2_TOP2ACenter_T6@assays@data@listData$counts
colnames(CCdata_sub_DucT2_TOP2ACenter_T6) = cds_sub_DucT2_TOP2ACenter_T6@assays@data@listData[["counts"]]@Dimnames[[2]]
rownames(CCdata_sub_DucT2_TOP2ACenter_T6) = cds_sub_DucT2_TOP2ACenter_T6@assays@data@listData[["counts"]]@Dimnames[[1]]


DataCellcycle_sub_DucT2_TOP2ACenter_T6 <- CCdata_sub_DucT2_TOP2ACenter_T6
DataCellcycle_sub_DucT2_TOP2ACenter_T6 <- as(as.matrix(DataCellcycle_sub_DucT2_TOP2ACenter_T6), "dgCMatrix")

marrow_sub_DucT2_TOP2ACenter_T6 <- CreateSeuratObject(counts = DataCellcycle_sub_DucT2_TOP2ACenter_T6)
marrow_sub_DucT2_TOP2ACenter_T6 <- NormalizeData(marrow_sub_DucT2_TOP2ACenter_T6)
marrow_sub_DucT2_TOP2ACenter_T6 <- FindVariableFeatures(marrow_sub_DucT2_TOP2ACenter_T6, selection.method = "vst")
marrow_sub_DucT2_TOP2ACenter_T6 <- ScaleData(marrow_sub_DucT2_TOP2ACenter_T6, features = rownames(marrow_sub_DucT2_TOP2ACenter_T6))
#錯誤: 無法配置大小為 3.9 Gb 的向量
#https://d.cosx.org/d/413001-413001
#memory.limit(15000)


marrow_sub_DucT2_TOP2ACenter_T6 <- RunPCA(marrow_sub_DucT2_TOP2ACenter_T6, features = VariableFeatures(marrow_sub_DucT2_TOP2ACenter_T6), ndims.print = 1:10, nfeatures.print = 25)

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T6_DimHeatmap1.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T6, dims = c(1,2),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T6_DimHeatmap2.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T6, dims = c(3,4),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T6_DimHeatmap3.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T6, dims = c(5,6),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T6_DimHeatmap4.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T6, dims = c(7,8),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T6_DimHeatmap5.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T6, dims = c(9,10),nfeatures = 30)
dev.off() # 關閉輸出圖檔

DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T6, dims = c(1:2),nfeatures = 50)

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T6_DimHeatmapMa.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T6, dims = c(1:18))
dev.off() # 關閉輸出圖檔


########################  DucT2_TOP2ACenter trajectories_T7 ##########################
cds_sub_DucT2_TOP2ACenter_T7 <- choose_graph_segments(cds_sub_DucT2 ,clear_cds = FALSE)

plot_cells(cds_sub_DucT2_TOP2ACenter_T7, color_cells_by="cell_cycle",cell_size=2, label_cell_groups=FALSE) + scale_color_manual(values = colorsT)

CCdata_sub_DucT2_TOP2ACenter_T7 <- cds_sub_DucT2_TOP2ACenter_T7@assays@data@listData$counts
colnames(CCdata_sub_DucT2_TOP2ACenter_T7) = cds_sub_DucT2_TOP2ACenter_T7@assays@data@listData[["counts"]]@Dimnames[[2]]
rownames(CCdata_sub_DucT2_TOP2ACenter_T7) = cds_sub_DucT2_TOP2ACenter_T7@assays@data@listData[["counts"]]@Dimnames[[1]]


DataCellcycle_sub_DucT2_TOP2ACenter_T7 <- CCdata_sub_DucT2_TOP2ACenter_T7
DataCellcycle_sub_DucT2_TOP2ACenter_T7 <- as(as.matrix(DataCellcycle_sub_DucT2_TOP2ACenter_T7), "dgCMatrix")

marrow_sub_DucT2_TOP2ACenter_T7 <- CreateSeuratObject(counts = DataCellcycle_sub_DucT2_TOP2ACenter_T7)
marrow_sub_DucT2_TOP2ACenter_T7 <- NormalizeData(marrow_sub_DucT2_TOP2ACenter_T7)
marrow_sub_DucT2_TOP2ACenter_T7 <- FindVariableFeatures(marrow_sub_DucT2_TOP2ACenter_T7, selection.method = "vst")
marrow_sub_DucT2_TOP2ACenter_T7 <- ScaleData(marrow_sub_DucT2_TOP2ACenter_T7, features = rownames(marrow_sub_DucT2_TOP2ACenter_T7))
#錯誤: 無法配置大小為 3.9 Gb 的向量
#https://d.cosx.org/d/413001-413001
#memory.limit(15000)


marrow_sub_DucT2_TOP2ACenter_T7 <- RunPCA(marrow_sub_DucT2_TOP2ACenter_T7, features = VariableFeatures(marrow_sub_DucT2_TOP2ACenter_T7), ndims.print = 1:10, nfeatures.print = 25)

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T7_DimHeatmap1.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T7, dims = c(1,2),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T7_DimHeatmap2.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T7, dims = c(3,4),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T7_DimHeatmap3.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T7, dims = c(5,6),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T7_DimHeatmap4.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T7, dims = c(7,8),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T7_DimHeatmap5.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T7, dims = c(9,10),nfeatures = 30)
dev.off() # 關閉輸出圖檔

DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T7, dims = c(1:2),nfeatures = 50)

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T7_DimHeatmapMa.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T7, dims = c(1:18))
dev.off() # 關閉輸出圖檔

########################  DucT2_TOP2ACenter trajectories_T8 ##########################
cds_sub_DucT2_TOP2ACenter_T8 <- choose_graph_segments(cds_sub_DucT2 ,clear_cds = FALSE)

plot_cells(cds_sub_DucT2_TOP2ACenter_T8, color_cells_by="cell_cycle",cell_size=2, label_cell_groups=FALSE) + scale_color_manual(values = colorsT)

CCdata_sub_DucT2_TOP2ACenter_T8 <- cds_sub_DucT2_TOP2ACenter_T8@assays@data@listData$counts
colnames(CCdata_sub_DucT2_TOP2ACenter_T8) = cds_sub_DucT2_TOP2ACenter_T8@assays@data@listData[["counts"]]@Dimnames[[2]]
rownames(CCdata_sub_DucT2_TOP2ACenter_T8) = cds_sub_DucT2_TOP2ACenter_T8@assays@data@listData[["counts"]]@Dimnames[[1]]


DataCellcycle_sub_DucT2_TOP2ACenter_T8 <- CCdata_sub_DucT2_TOP2ACenter_T8
DataCellcycle_sub_DucT2_TOP2ACenter_T8 <- as(as.matrix(DataCellcycle_sub_DucT2_TOP2ACenter_T8), "dgCMatrix")

marrow_sub_DucT2_TOP2ACenter_T8 <- CreateSeuratObject(counts = DataCellcycle_sub_DucT2_TOP2ACenter_T8)
marrow_sub_DucT2_TOP2ACenter_T8 <- NormalizeData(marrow_sub_DucT2_TOP2ACenter_T8)
marrow_sub_DucT2_TOP2ACenter_T8 <- FindVariableFeatures(marrow_sub_DucT2_TOP2ACenter_T8, selection.method = "vst")
marrow_sub_DucT2_TOP2ACenter_T8 <- ScaleData(marrow_sub_DucT2_TOP2ACenter_T8, features = rownames(marrow_sub_DucT2_TOP2ACenter_T8))
#錯誤: 無法配置大小為 3.9 Gb 的向量
#https://d.cosx.org/d/413001-413001
#memory.limit(15000)


marrow_sub_DucT2_TOP2ACenter_T8 <- RunPCA(marrow_sub_DucT2_TOP2ACenter_T8, features = VariableFeatures(marrow_sub_DucT2_TOP2ACenter_T8), ndims.print = 1:10, nfeatures.print = 25)

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T8_DimHeatmap1.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T8, dims = c(1,2),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T8_DimHeatmap2.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T8, dims = c(3,4),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T8_DimHeatmap3.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T8, dims = c(5,6),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T8_DimHeatmap4.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T8, dims = c(7,8),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T8_DimHeatmap5.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T8, dims = c(9,10),nfeatures = 30)
dev.off() # 關閉輸出圖檔

DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T8, dims = c(1:2),nfeatures = 50)

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T8_DimHeatmapMa.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T8, dims = c(1:18))
dev.off() # 關閉輸出圖檔


