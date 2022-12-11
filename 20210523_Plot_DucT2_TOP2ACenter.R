########################  DucT2_TOP2ACenter ##########################
cds_sub_DucT2_TOP2ACenter <- choose_cells(cds_subset_NewK)
plot_cells(cds_sub_DucT2_TOP2ACenter, label_cell_groups=FALSE)
plot_cells(cds_sub_DucT2_TOP2ACenter, label_cell_groups=FALSE, show_trajectory_graph = FALSE, cell_size = 2)
plot_cells(cds_sub_DucT2_TOP2ACenter, label_cell_groups=FALSE, show_trajectory_graph = FALSE, cell_size = 2, color_cells_by="cell_cycle")+ scale_color_manual(values = colors_cc)
plot_cells(cds_sub_DucT2_TOP2ACenter, genes=c(Main), label_cell_groups=FALSE, show_trajectory_graph = FALSE, cell_size = 2)
#plot_cells(cds_sub_DucT2_TOP2ACenter, genes=c(Main),cell_size=1,label_cell_groups = FALSE, show_trajectory_graph = FALSE)

###### Convert Monocle3 Object to Seurat Object ######
# # getFilePath("Monocle3_To_Seurat.R")
# marrow_sub_DucT2_TOP2ACenter <- Monocle3_To_Seurat(cds_sub_DucT2_TOP2ACenter,"sub_DT2TOP2ACTR") #sub_DT2TOP2ACTR:sub_DucT2_TOP2ACenter
# 
# ###### Assign Cell-Cycle Scores ######
# #getFilePath("Cell-Cycle Scoring and Regression.R")
# marrow_sub_DucT2_TOP2ACenter <- CCScorReg(GeneNAFMT,marrow_sub_DucT2_TOP2ACenter) #這個function存在於Cell-Cycle Scoring and Regression.R裡面
# 
# RidgePlot(marrow_sub_DucT2_TOP2ACenter,cols = colors_cc, features = c(Main), ncol = 1)
# RidgePlot(marrow_sub_DucT2_TOP2ACenter,cols = colors_cc, features = c(Main_Group), ncol = 2,log=TRUE) 
# RidgePlot(marrow_sub_DucT2_TOP2ACenter,cols = colors_cc, features = c(Main_Group), ncol = 2,y.max = 100) 

CCdata_sub_DucT2_TOP2ACenter <- cds_sub_DucT2_TOP2ACenter@assays@data@listData$counts
colnames(CCdata_sub_DucT2_TOP2ACenter) = cds_sub_DucT2_TOP2ACenter@assays@data@listData[["counts"]]@Dimnames[[2]]
rownames(CCdata_sub_DucT2_TOP2ACenter) = cds_sub_DucT2_TOP2ACenter@assays@data@listData[["counts"]]@Dimnames[[1]]


DataCellcycle_sub_DucT2_TOP2ACenter <- CCdata_sub_DucT2_TOP2ACenter
DataCellcycle_sub_DucT2_TOP2ACenter <- as(as.matrix(DataCellcycle_sub_DucT2_TOP2ACenter), "dgCMatrix")

marrow_sub_DucT2_TOP2ACenter <- CreateSeuratObject(counts = DataCellcycle_sub_DucT2_TOP2ACenter)
marrow_sub_DucT2_TOP2ACenter <- NormalizeData(marrow_sub_DucT2_TOP2ACenter)
marrow_sub_DucT2_TOP2ACenter <- FindVariableFeatures(marrow_sub_DucT2_TOP2ACenter, selection.method = "vst")
marrow_sub_DucT2_TOP2ACenter <- ScaleData(marrow_sub_DucT2_TOP2ACenter, features = rownames(marrow_sub_DucT2_TOP2ACenter))


marrow_sub_DucT2_TOP2ACenter <- RunPCA(marrow_sub_DucT2_TOP2ACenter, features = VariableFeatures(marrow_sub_DucT2_TOP2ACenter), ndims.print = 6:10, nfeatures.print = 10)
marrow_sub_DucT2_TOP2ACenter@assays[["RNA"]]@counts@Dimnames[[1]] <- cds_sub_DucT2_TOP2ACenter@assays@data@listData[["counts"]]@Dimnames[[1]]
marrow_sub_DucT2_TOP2ACenter@assays[["RNA"]]@data@Dimnames[[1]] <- cds_sub_DucT2_TOP2ACenter@assays@data@listData[["counts"]]@Dimnames[[1]]

#png(paste0(PathName,"/",RVersion,"/",RVersion,"_","DucT2_TOP2ACenter_CcDimHeatmap.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter, dims = c(3, 4))
#dev.off() # 關閉輸出圖檔

#png(paste0(PathName,"/",RVersion,"/",RVersion,"_","DucT2_TOP2ACenter_CcDimHeatmap2.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter, dims = c(1, 2))
#dev.off() # 關閉輸出圖檔

## Assign Cell-Cycle Scores
marrow_sub_DucT2_TOP2ACenter <- CellCycleScoring(marrow_sub_DucT2_TOP2ACenter, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(marrow_sub_DucT2_TOP2ACenter[[]])


## Plot the RidgePlot
RidgePlot(marrow_sub_DucT2_TOP2ACenter,cols = colorsT, features = c("TOP2A"), ncol = 1)
RidgePlot(marrow_sub_DucT2_TOP2ACenter,cols = colorsT, features = c(Main), ncol = 2,log=TRUE)
RidgePlot(marrow_sub_DucT2_TOP2ACenter,cols = colorsT, features = c(Main), ncol = 2,y.max = 100)

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_RidgePlot_Sub_DucT2_TOP2ACenter_V2.png")) # 設定輸出圖檔
RidgePlot(marrow_sub_DucT2_TOP2ACenter,cols = colorsT, features = c(Main), ncol = 1)
dev.off() # 關閉輸出圖檔


## 將Seurat跑出的Cell cycle結果寫入Monocle3的cds檔
cds_sub_DucT2_TOP2ACenter@colData@listData$cell_cycle <- marrow_sub_DucT2_TOP2ACenter@active.ident

## Plot the Violin Plot 
Maingroup_ciliated_genes <- c("TOP2A")
cds_sub_DucT2_TOP2ACenter <- cds_sub_DucT2_TOP2ACenter[rowData(cds_sub_DucT2_TOP2ACenter)$gene_short_name %in% Maingroup_ciliated_genes,]
plot_genes_violin(cds_sub_DucT2_TOP2ACenter, group_cells_by="cell_cycle", ncol=2, log_scale = FALSE) +
  theme(axis.text.x=element_text(angle=45, hjust=1))


png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_Violin_Main_Sub_DucT2_TOP2ACenter_V2.png")) # 設定輸出圖檔
plot_genes_violin(cds_sub_DucT2_TOP2ACenter, group_cells_by="cell_cycle", ncol=2, log_scale = FALSE) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
dev.off() # 關閉輸出圖檔

##
png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_CellCycle_Sub_DucT2_TOP2ACenter_V2.png")) # 設定輸出圖檔
plot_cells(cds_sub_DucT2_TOP2ACenter, color_cells_by="cell_cycle",cell_size=3, label_cell_groups=FALSE, show_trajectory_graph = FALSE) + scale_color_manual(values = colorsT)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_Sub_DucT2_TOP2ACenter_V2.png")) # 設定輸出圖檔
plot_cells(cds_sub_DucT2_TOP2ACenter, genes=c("TOP2A"),cell_size=3, label_cell_groups=FALSE, show_trajectory_graph = FALSE)
dev.off() # 關閉輸出圖檔


plot_cells(cds_sub_DucT2_TOP2ACenter, color_cells_by="cell_cycle",cell_size=2, label_cell_groups=FALSE) + scale_color_manual(values = colorsT)



plot_cells(cds_sub_DucT2_TOP2ACenter, genes=c("TOP2A"),cell_size=2, label_cell_groups=FALSE, show_trajectory_graph = FALSE)

AFD_lineage_cds_sub_DucT2_TOP2ACenter <- cds_sub_DucT2_TOP2ACenter[rowData(cds_sub_DucT2_TOP2ACenter)$gene_short_name %in% Main]
plot_genes_in_pseudotime(AFD_lineage_cds_sub_DucT2_TOP2ACenter,
                         color_cells_by="cell_cycle",cell_size=2,
                         min_expr=0.5)+ scale_color_manual(values = colorsT)




########################  DucT2_TOP2ACenter trajectories ##########################
cds_sub_DucT2_TOP2ACenter_T1 <- choose_graph_segments(cds_sub_DucT2 ,clear_cds = FALSE)

plot_cells(cds_sub_DucT2_TOP2ACenter_T1, color_cells_by="cell_cycle",cell_size=2, label_cell_groups=FALSE) + scale_color_manual(values = colors_cc)


CCdata_sub_DucT2_TOP2ACenter_T1 <- cds_sub_DucT2_TOP2ACenter_T1@assays@data@listData$counts
colnames(CCdata_sub_DucT2_TOP2ACenter_T1) = cds_sub_DucT2_TOP2ACenter_T1@assays@data@listData[["counts"]]@Dimnames[[2]]
rownames(CCdata_sub_DucT2_TOP2ACenter_T1) = cds_sub_DucT2_TOP2ACenter_T1@assays@data@listData[["counts"]]@Dimnames[[1]]


DataCellcycle_sub_DucT2_TOP2ACenter_T1 <- CCdata_sub_DucT2_TOP2ACenter_T1
DataCellcycle_sub_DucT2_TOP2ACenter_T1 <- as(as.matrix(DataCellcycle_sub_DucT2_TOP2ACenter_T1), "dgCMatrix")

marrow_sub_DucT2_TOP2ACenter_T1 <- CreateSeuratObject(counts = DataCellcycle_sub_DucT2_TOP2ACenter_T1)
marrow_sub_DucT2_TOP2ACenter_T1 <- NormalizeData(marrow_sub_DucT2_TOP2ACenter_T1)
marrow_sub_DucT2_TOP2ACenter_T1 <- FindVariableFeatures(marrow_sub_DucT2_TOP2ACenter_T1, selection.method = "vst")
marrow_sub_DucT2_TOP2ACenter_T1 <- ScaleData(marrow_sub_DucT2_TOP2ACenter_T1, features = rownames(marrow_sub_DucT2_TOP2ACenter_T1))
#錯誤: 無法配置大小為 3.9 Gb 的向量
#https://d.cosx.org/d/413001-413001
#memory.limit(15000)


marrow_sub_DucT2_TOP2ACenter_T1 <- RunPCA(marrow_sub_DucT2_TOP2ACenter_T1, features = VariableFeatures(marrow_sub_DucT2_TOP2ACenter_T1), ndims.print = 1:10, nfeatures.print = 25)
# marrow_sub_DucT2_TOP2ACenter_T1@assays[["RNA"]]@counts@Dimnames[[1]] <- cds_sub_DucT2_TOP2ACenter_T1@assays@data@listData[["counts"]]@Dimnames[[1]]
# marrow_sub_DucT2_TOP2ACenter_T1@assays[["RNA"]]@data@Dimnames[[1]] <- cds_sub_DucT2_TOP2ACenter_T1@assays@data@listData[["counts"]]@Dimnames[[1]]
png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T1_DimHeatmap1.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T1, dims = c(1,2),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T1_DimHeatmap2.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T1, dims = c(3,4),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T1_DimHeatmap3.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T1, dims = c(5,6),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T1_DimHeatmap4.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T1, dims = c(7,8),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T1_DimHeatmap5.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T1, dims = c(9,10),nfeatures = 30)
dev.off() # 關閉輸出圖檔

DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T1, dims = c(1:2),nfeatures = 50)

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_T1_DimHeatmapMa.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter_T1, dims = c(1:18))
dev.off() # 關閉輸出圖檔

