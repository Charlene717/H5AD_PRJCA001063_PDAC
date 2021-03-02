######################## cds_sub_Normal_TOP2A_Traj  ##########################
#cds_sub_DucT2_TOP2A_Cluster7 <- choose_graph_segments(cds_sub_DucT2 ,clear_cds = FALSE)
cds_sub_Normal_TOP2A_Traj <- choose_cells(cds ,clear_cds = FALSE)

plot_cells(cds_sub_Normal_TOP2A_Traj, color_cells_by="cell_cycle",cell_size=2, label_cell_groups=FALSE, show_trajectory_graph = FALSE) + scale_color_manual(values = colorsT)


CCdata_sub_Normal_TOP2A_Traj <- cds_sub_Normal_TOP2A_Traj@assays@data@listData$counts
colnames(CCdata_sub_Normal_TOP2A_Traj) = cds_sub_Normal_TOP2A_Traj@assays@data@listData[["counts"]]@Dimnames[[2]]
rownames(CCdata_sub_Normal_TOP2A_Traj) = cds_sub_Normal_TOP2A_Traj@assays@data@listData[["counts"]]@Dimnames[[1]]


DataCellcycle_sub_Normal_TOP2A_Traj <- CCdata_sub_Normal_TOP2A_Traj
DataCellcycle_sub_Normal_TOP2A_Traj <- as(as.matrix(DataCellcycle_sub_Normal_TOP2A_Traj), "dgCMatrix")

marrow_sub_Normal_TOP2A_Traj <- CreateSeuratObject(counts = DataCellcycle_sub_Normal_TOP2A_Traj)
marrow_sub_Normal_TOP2A_Traj <- NormalizeData(marrow_sub_Normal_TOP2A_Traj)
marrow_sub_Normal_TOP2A_Traj <- FindVariableFeatures(marrow_sub_Normal_TOP2A_Traj, selection.method = "vst")
marrow_sub_Normal_TOP2A_Traj <- ScaleData(marrow_sub_Normal_TOP2A_Traj, features = rownames(marrow_sub_Normal_TOP2A_Traj))
#錯誤: 無法配置大小為 3.9 Gb 的向量
#https://d.cosx.org/d/413001-413001
#memory.limit(15000)


marrow_sub_Normal_TOP2A_Traj <- RunPCA(marrow_sub_Normal_TOP2A_Traj, features = VariableFeatures(marrow_sub_Normal_TOP2A_Traj), ndims.print = 1:10, nfeatures.print = 25)

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_Normal_TOP2A_Traj_DimHeatmap1.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_Normal_TOP2A_Traj, dims = c(1,2),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_Normal_TOP2A_Traj_DimHeatmap2.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_Normal_TOP2A_Traj, dims = c(3,4),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_Normal_TOP2A_Traj_DimHeatmap3.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_Normal_TOP2A_Traj, dims = c(5,6),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_Normal_TOP2A_Traj_DimHeatmap4.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_Normal_TOP2A_Traj, dims = c(7,8),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_Normal_TOP2A_Traj_DimHeatmap5.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_Normal_TOP2A_Traj, dims = c(9,10),nfeatures = 30)
dev.off() # 關閉輸出圖檔

DimHeatmap(marrow_sub_Normal_TOP2A_Traj, dims = c(1:2),nfeatures = 50)

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_Normal_TOP2A_Traj_DimHeatmapMa.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_Normal_TOP2A_Traj, dims = c(1:18))
dev.off() # 關閉輸出圖檔


################################################




# Also read in a list of cell cycle markers, from Tirosh et al, 2015
# cc.genes <- readLines(con = "~/Downloads/seurat_resources/regev_lab_cell_cycle_genes.txt")
# cc.genes <- read.csv(paste0(PathName,"/Cell cycle/regev_lab_cell_cycle_genesCh.csv"))
# # We can segregate this list into markers of G2/M phase and markers of S
# # phase
# s.genes <- cc.genes[1:43,]
# s.genes2 <- capitalize(tolower(s.genes))
# # https://www.360kuai.com/pc/9ae634cb511f34aee?cota=4&kuai_so=1&tj_url=so_rec&sign=360_7bc3b157
# library(tidyverse)
# G_listCCs.genes <- AnnotationDbi::select(org.Mm.eg.db, keys=s.genes2, columns='ENSEMBL', keytype='SYMBOL')
# # G_listCCs.genes <- getBM(filters= "hgnc_symbol", attributes= c("ensembl_gene_id","hgnc_symbol"),values=s.genes,mart= martCCG,uniqueRows= TRUE,
# #                          quote = "\"")
# G_listCCs.genes2 = G_listCCs.genes[!duplicated(G_listCCs.genes[2]),]
# G_listCCs.genes3 <- na.omit(G_listCCs.genes2[2])
# G_listCCs.genes4 <- G_listCCs.genes3[,1]
# # (Ori)s.genes <- cc.genes[,1:43]
# g2m.genes <- cc.genes[44:97,]
# g2m.genes2 <- capitalize(tolower(g2m.genes))
# G_listCCg2m.genes <- AnnotationDbi::select(org.Mm.eg.db, keys=g2m.genes2, columns='ENSEMBL', keytype='SYMBOL')
# G_listCCg2m.genes2 = G_listCCg2m.genes[!duplicated(G_listCCg2m.genes[2]),]
# G_listCCg2m.genes3 <- na.omit(G_listCCg2m.genes2[2])
# G_listCCg2m.genes4 <-G_listCCg2m.genes3[,1]
# #marrow <- CellCycleScoring(marrow, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# #marrow <- CellCycleScoring(marrow, s.features = G_listCCs.genes2, g2m.features = G_listCCg2m.genes2, set.ident = TRUE)
# 
# 
# 
# marrow_sub_DucT2_TOP2ACenter <- CellCycleScoring(marrow_sub_DucT2_TOP2ACenter, s.features = G_listCCs.genes4, g2m.features = G_listCCg2m.genes4, set.ident = TRUE)
marrow_sub_Normal_TOP2A_Traj <- CellCycleScoring(marrow_sub_Normal_TOP2A_Traj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(marrow_sub_Normal_TOP2A_Traj[[]])


##############
colorsT <- c("#FF9912B3", "#32CD3299", "#4169E1B3")
RidgePlot(marrow_sub_Normal_TOP2A_Traj,cols = colorsT, features = c("TOP2A"), ncol = 1)
Main=("TOP2A")

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_RidgePlot_Sub_Normal_TOP2A_Traj_V2.png")) # 設定輸出圖檔

RidgePlot(marrow_sub_Normal_TOP2A_Traj,cols = colorsT, features = c(Main), ncol = 1)
dev.off() # 關閉輸出圖檔

RidgePlot(marrow_sub_Normal_TOP2A_Traj,cols = colorsT, features = c("TOP2A"), ncol = 1)
# https://datavizpyr.com/ridgeline-plot-with-ggridges-in-r/
# https://www.r-graph-gallery.com/294-basic-ridgeline-plot.html
RidgePlot(marrow_sub_Normal_TOP2A_Traj,cols = colorsT, features = c(Main), ncol = 2,log=TRUE)
RidgePlot(marrow_sub_Normal_TOP2A_Traj,cols = colorsT, features = c(Main), ncol = 2,y.max = 100)

#將Seurat跑出的Cell cycle結果寫入Monocle3的cds檔
cds_sub_Normal_TOP2A_Traj@colData@listData$cell_cycle <- marrow_sub_Normal_TOP2A_Traj@active.ident


##
Maingroup_ciliated_genes <- c("TOP2A")
cds_sub_Normal_TOP2A_Traj <- cds_sub_Normal_TOP2A_Traj[rowData(cds_sub_Normal_TOP2A_Traj)$gene_short_name %in% Maingroup_ciliated_genes,]

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_Violin_Main_Sub_Normal_TOP2A_Traj_V2.png")) # 設定輸出圖檔
plot_genes_violin(cds_sub_Normal_TOP2A_Traj, group_cells_by="cell_cycle", ncol=2, log_scale = FALSE) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
dev.off() # 關閉輸出圖檔

##
png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_CellCycle_Sub_Normal_TOP2A_Traj_V2.png")) # 設定輸出圖檔
plot_cells(cds_sub_Normal_TOP2A_Traj, color_cells_by="cell_cycle",cell_size=3, label_cell_groups=FALSE, show_trajectory_graph = FALSE) + scale_color_manual(values = colorsT)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_Sub_Normal_TOP2A_Traj_V2.png")) # 設定輸出圖檔
plot_cells(cds_sub_Normal_TOP2A_Traj, genes=c("TOP2A"),cell_size=3, label_cell_groups=FALSE, show_trajectory_graph = FALSE)
dev.off() # 關閉輸出圖檔


plot_cells(cds_sub_Normal_TOP2A_Traj, color_cells_by="cell_cycle",cell_size=2, label_cell_groups=FALSE) + scale_color_manual(values = colorsT)



plot_cells(cds_sub_Normal_TOP2A_Traj, genes=c("TOP2A"),cell_size=2, label_cell_groups=FALSE, show_trajectory_graph = FALSE)

# AFD_lineage_cds_sub_DucT2_TOP2A_Cluster7 <- cds_sub_DucT2_TOP2A_Cluster7[rowData(cds_sub_DucT2_TOP2A_Cluster7)$gene_short_name %in% Main]
# plot_genes_in_pseudotime(AFD_lineage_cds_sub_DucT2_TOP2A_Cluster7,
#                          color_cells_by="cell_cycle",cell_size=2,
#                          min_expr=0.5)+ scale_color_manual(values = colorsT)




########################  cds_sub_Abnormal_TOP2A_Traj  ##########################
#cds_sub_DucT2_TOP2A_Cluster7 <- choose_graph_segments(cds_sub_DucT2 ,clear_cds = FALSE)
cds_sub_Abnormal_TOP2A_Traj <- choose_cells(cds ,clear_cds = FALSE)

plot_cells(cds_sub_Abnormal_TOP2A_Traj, color_cells_by="cell_cycle",cell_size=2, label_cell_groups=FALSE, show_trajectory_graph = FALSE) + scale_color_manual(values = colorsT)


CCdata_sub_Abnormal_TOP2A_Traj <- cds_sub_Abnormal_TOP2A_Traj@assays@data@listData$counts
colnames(CCdata_sub_Abnormal_TOP2A_Traj) = cds_sub_Abnormal_TOP2A_Traj@assays@data@listData[["counts"]]@Dimnames[[2]]
rownames(CCdata_sub_Abnormal_TOP2A_Traj) = cds_sub_Abnormal_TOP2A_Traj@assays@data@listData[["counts"]]@Dimnames[[1]]


DataCellcycle_sub_Abnormal_TOP2A_Traj <- CCdata_sub_Abnormal_TOP2A_Traj
DataCellcycle_sub_Abnormal_TOP2A_Traj <- as(as.matrix(DataCellcycle_sub_Abnormal_TOP2A_Traj), "dgCMatrix")

marrow_sub_Abnormal_TOP2A_Traj <- CreateSeuratObject(counts = DataCellcycle_sub_Abnormal_TOP2A_Traj)
marrow_sub_Abnormal_TOP2A_Traj <- NormalizeData(marrow_sub_Abnormal_TOP2A_Traj)
marrow_sub_Abnormal_TOP2A_Traj <- FindVariableFeatures(marrow_sub_Abnormal_TOP2A_Traj, selection.method = "vst")
marrow_sub_Abnormal_TOP2A_Traj <- ScaleData(marrow_sub_Abnormal_TOP2A_Traj, features = rownames(marrow_sub_Abnormal_TOP2A_Traj))
#錯誤: 無法配置大小為 3.9 Gb 的向量
#https://d.cosx.org/d/413001-413001
#memory.limit(15000)


marrow_sub_Abnormal_TOP2A_Traj <- RunPCA(marrow_sub_Abnormal_TOP2A_Traj, features = VariableFeatures(marrow_sub_Abnormal_TOP2A_Traj), ndims.print = 1:10, nfeatures.print = 25)

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_Abnormal_TOP2A_Traj_DimHeatmap1.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_Abnormal_TOP2A_Traj, dims = c(1,2),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_Abnormal_TOP2A_Traj_DimHeatmap2.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_Abnormal_TOP2A_Traj, dims = c(3,4),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_Abnormal_TOP2A_Traj_DimHeatmap3.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_Abnormal_TOP2A_Traj, dims = c(5,6),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_Abnormal_TOP2A_Traj_DimHeatmap4.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_Abnormal_TOP2A_Traj, dims = c(7,8),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_Abnormal_TOP2A_Traj_DimHeatmap5.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_Abnormal_TOP2A_Traj, dims = c(9,10),nfeatures = 30)
dev.off() # 關閉輸出圖檔

DimHeatmap(marrow_sub_Abnormal_TOP2A_Traj, dims = c(1:2),nfeatures = 50)

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_Abnormal_TOP2A_Traj_DimHeatmapMa.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_Abnormal_TOP2A_Traj, dims = c(1:18))
dev.off() # 關閉輸出圖檔


################################################




# Also read in a list of cell cycle markers, from Tirosh et al, 2015
# cc.genes <- readLines(con = "~/Downloads/seurat_resources/regev_lab_cell_cycle_genes.txt")
# cc.genes <- read.csv(paste0(PathName,"/Cell cycle/regev_lab_cell_cycle_genesCh.csv"))
# # We can segregate this list into markers of G2/M phase and markers of S
# # phase
# s.genes <- cc.genes[1:43,]
# s.genes2 <- capitalize(tolower(s.genes))
# # https://www.360kuai.com/pc/9ae634cb511f34aee?cota=4&kuai_so=1&tj_url=so_rec&sign=360_7bc3b157
# library(tidyverse)
# G_listCCs.genes <- AnnotationDbi::select(org.Mm.eg.db, keys=s.genes2, columns='ENSEMBL', keytype='SYMBOL')
# # G_listCCs.genes <- getBM(filters= "hgnc_symbol", attributes= c("ensembl_gene_id","hgnc_symbol"),values=s.genes,mart= martCCG,uniqueRows= TRUE,
# #                          quote = "\"")
# G_listCCs.genes2 = G_listCCs.genes[!duplicated(G_listCCs.genes[2]),]
# G_listCCs.genes3 <- na.omit(G_listCCs.genes2[2])
# G_listCCs.genes4 <- G_listCCs.genes3[,1]
# # (Ori)s.genes <- cc.genes[,1:43]
# g2m.genes <- cc.genes[44:97,]
# g2m.genes2 <- capitalize(tolower(g2m.genes))
# G_listCCg2m.genes <- AnnotationDbi::select(org.Mm.eg.db, keys=g2m.genes2, columns='ENSEMBL', keytype='SYMBOL')
# G_listCCg2m.genes2 = G_listCCg2m.genes[!duplicated(G_listCCg2m.genes[2]),]
# G_listCCg2m.genes3 <- na.omit(G_listCCg2m.genes2[2])
# G_listCCg2m.genes4 <-G_listCCg2m.genes3[,1]
# #marrow <- CellCycleScoring(marrow, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# #marrow <- CellCycleScoring(marrow, s.features = G_listCCs.genes2, g2m.features = G_listCCg2m.genes2, set.ident = TRUE)
# 
# 
# 
# marrow_sub_DucT2_TOP2ACenter <- CellCycleScoring(marrow_sub_DucT2_TOP2ACenter, s.features = G_listCCs.genes4, g2m.features = G_listCCg2m.genes4, set.ident = TRUE)
marrow_sub_Abnormal_TOP2A_Traj <- CellCycleScoring(marrow_sub_Abnormal_TOP2A_Traj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(marrow_sub_Abnormal_TOP2A_Traj[[]])


##############
colorsT <- c("#FF9912B3", "#32CD3299", "#4169E1B3")
RidgePlot(marrow_sub_Abnormal_TOP2A_Traj,cols = colorsT, features = c("TOP2A"), ncol = 1)
Main=("TOP2A")

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_RidgePlot_Sub_Abnormal_TOP2A_Traj_V2.png")) # 設定輸出圖檔

RidgePlot(marrow_sub_Abnormal_TOP2A_Traj,cols = colorsT, features = c(Main), ncol = 1)
dev.off() # 關閉輸出圖檔

RidgePlot(marrow_sub_Abnormal_TOP2A_Traj,cols = colorsT, features = c("TOP2A"), ncol = 1)
# https://datavizpyr.com/ridgeline-plot-with-ggridges-in-r/
# https://www.r-graph-gallery.com/294-basic-ridgeline-plot.html
RidgePlot(marrow_sub_Abnormal_TOP2A_Traj,cols = colorsT, features = c(Main), ncol = 2,log=TRUE)
RidgePlot(marrow_sub_Abnormal_TOP2A_Traj,cols = colorsT, features = c(Main), ncol = 2,y.max = 100)

#將Seurat跑出的Cell cycle結果寫入Monocle3的cds檔
cds_sub_Abnormal_TOP2A_Traj@colData@listData$cell_cycle <- marrow_sub_Abnormal_TOP2A_Traj@active.ident


##
Maingroup_ciliated_genes <- c("TOP2A")
cds_sub_Abnormal_TOP2A_Traj <- cds_sub_Abnormal_TOP2A_Traj[rowData(cds_sub_Abnormal_TOP2A_Traj)$gene_short_name %in% Maingroup_ciliated_genes,]

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_Violin_Main_Sub_Abnormal_TOP2A_Traj_V2.png")) # 設定輸出圖檔
plot_genes_violin(cds_sub_Abnormal_TOP2A_Traj, group_cells_by="cell_cycle", ncol=2, log_scale = FALSE) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
dev.off() # 關閉輸出圖檔

##
png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_CellCycle_Sub_Abnormal_TOP2A_Traj_V2.png")) # 設定輸出圖檔
plot_cells(cds_sub_Abnormal_TOP2A_Traj, color_cells_by="cell_cycle",cell_size=3, label_cell_groups=FALSE, show_trajectory_graph = FALSE) + scale_color_manual(values = colorsT)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_Sub_Abnormal_TOP2A_Traj_V2.png")) # 設定輸出圖檔
plot_cells(cds_sub_Abnormal_TOP2A_Traj, genes=c("TOP2A"),cell_size=3, label_cell_groups=FALSE, show_trajectory_graph = FALSE)
dev.off() # 關閉輸出圖檔


plot_cells(cds_sub_Abnormal_TOP2A_Traj, color_cells_by="cell_cycle",cell_size=2, label_cell_groups=FALSE) + scale_color_manual(values = colorsT)



plot_cells(cds_sub_Abnormal_TOP2A_Traj, genes=c("TOP2A"),cell_size=2, label_cell_groups=FALSE, show_trajectory_graph = FALSE)

# AFD_lineage_cds_sub_DucT2_TOP2A_Cluster7 <- cds_sub_DucT2_TOP2A_Cluster7[rowData(cds_sub_DucT2_TOP2A_Cluster7)$gene_short_name %in% Main]
# plot_genes_in_pseudotime(AFD_lineage_cds_sub_DucT2_TOP2A_Cluster7,
#                          color_cells_by="cell_cycle",cell_size=2,
#                          min_expr=0.5)+ scale_color_manual(values = colorsT)


########################################  cds_sub_Malignant_End1  ##########################################
#cds_sub_DucT2_TOP2A_Cluster7 <- choose_graph_segments(cds_sub_DucT2 ,clear_cds = FALSE)
cds_sub_Malignant_End1 <- choose_cells(cds ,clear_cds = FALSE)

plot_cells(cds_sub_Malignant_End1, color_cells_by="cell_cycle",cell_size=2, label_cell_groups=FALSE, show_trajectory_graph = FALSE) + scale_color_manual(values = colorsT)


CCdata_sub_Malignant_End1 <- cds_sub_Malignant_End1@assays@data@listData$counts
colnames(CCdata_sub_Malignant_End1) = cds_sub_Malignant_End1@assays@data@listData[["counts"]]@Dimnames[[2]]
rownames(CCdata_sub_Malignant_End1) = cds_sub_Malignant_End1@assays@data@listData[["counts"]]@Dimnames[[1]]


DataCellcycle_sub_Malignant_End1 <- CCdata_sub_Malignant_End1
DataCellcycle_sub_Malignant_End1 <- as(as.matrix(DataCellcycle_sub_Malignant_End1), "dgCMatrix")

marrow_sub_Malignant_End1 <- CreateSeuratObject(counts = DataCellcycle_sub_Malignant_End1)
marrow_sub_Malignant_End1 <- NormalizeData(marrow_sub_Malignant_End1)
marrow_sub_Malignant_End1 <- FindVariableFeatures(marrow_sub_Malignant_End1, selection.method = "vst")
marrow_sub_Malignant_End1 <- ScaleData(marrow_sub_Malignant_End1, features = rownames(marrow_sub_Malignant_End1))
#錯誤: 無法配置大小為 3.9 Gb 的向量
#https://d.cosx.org/d/413001-413001
#memory.limit(15000)


marrow_sub_Malignant_End1 <- RunPCA(marrow_sub_Malignant_End1, features = VariableFeatures(marrow_sub_Malignant_End1), ndims.print = 1:10, nfeatures.print = 25)

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_Malignant_End1_DimHeatmap1.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_Malignant_End1, dims = c(1,2),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_Malignant_End1_DimHeatmap2.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_Malignant_End1, dims = c(3,4),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_Malignant_End1_DimHeatmap3.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_Malignant_End1, dims = c(5,6),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_Malignant_End1_DimHeatmap4.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_Malignant_End1, dims = c(7,8),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_Malignant_End1_DimHeatmap5.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_Malignant_End1, dims = c(9,10),nfeatures = 30)
dev.off() # 關閉輸出圖檔

DimHeatmap(marrow_sub_Malignant_End1, dims = c(1:2),nfeatures = 50)

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_Malignant_End1_DimHeatmapMa.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_Malignant_End1, dims = c(1:18))
dev.off() # 關閉輸出圖檔


################################################




# Also read in a list of cell cycle markers, from Tirosh et al, 2015
# cc.genes <- readLines(con = "~/Downloads/seurat_resources/regev_lab_cell_cycle_genes.txt")
# cc.genes <- read.csv(paste0(PathName,"/Cell cycle/regev_lab_cell_cycle_genesCh.csv"))
# # We can segregate this list into markers of G2/M phase and markers of S
# # phase
# s.genes <- cc.genes[1:43,]
# s.genes2 <- capitalize(tolower(s.genes))
# # https://www.360kuai.com/pc/9ae634cb511f34aee?cota=4&kuai_so=1&tj_url=so_rec&sign=360_7bc3b157
# library(tidyverse)
# G_listCCs.genes <- AnnotationDbi::select(org.Mm.eg.db, keys=s.genes2, columns='ENSEMBL', keytype='SYMBOL')
# # G_listCCs.genes <- getBM(filters= "hgnc_symbol", attributes= c("ensembl_gene_id","hgnc_symbol"),values=s.genes,mart= martCCG,uniqueRows= TRUE,
# #                          quote = "\"")
# G_listCCs.genes2 = G_listCCs.genes[!duplicated(G_listCCs.genes[2]),]
# G_listCCs.genes3 <- na.omit(G_listCCs.genes2[2])
# G_listCCs.genes4 <- G_listCCs.genes3[,1]
# # (Ori)s.genes <- cc.genes[,1:43]
# g2m.genes <- cc.genes[44:97,]
# g2m.genes2 <- capitalize(tolower(g2m.genes))
# G_listCCg2m.genes <- AnnotationDbi::select(org.Mm.eg.db, keys=g2m.genes2, columns='ENSEMBL', keytype='SYMBOL')
# G_listCCg2m.genes2 = G_listCCg2m.genes[!duplicated(G_listCCg2m.genes[2]),]
# G_listCCg2m.genes3 <- na.omit(G_listCCg2m.genes2[2])
# G_listCCg2m.genes4 <-G_listCCg2m.genes3[,1]
# #marrow <- CellCycleScoring(marrow, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# #marrow <- CellCycleScoring(marrow, s.features = G_listCCs.genes2, g2m.features = G_listCCg2m.genes2, set.ident = TRUE)
# 
# 
# 
# marrow_sub_DucT2_TOP2ACenter <- CellCycleScoring(marrow_sub_DucT2_TOP2ACenter, s.features = G_listCCs.genes4, g2m.features = G_listCCg2m.genes4, set.ident = TRUE)
marrow_sub_Malignant_End1 <- CellCycleScoring(marrow_sub_Malignant_End1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(marrow_sub_Malignant_End1[[]])


##############
colorsT <- c("#FF9912B3", "#32CD3299", "#4169E1B3")
RidgePlot(marrow_sub_Malignant_End1,cols = colorsT, features = c("TOP2A"), ncol = 1)
Main=("TOP2A")

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_RidgePlot_Sub_Malignant_End1_V2.png")) # 設定輸出圖檔

RidgePlot(marrow_sub_Malignant_End1,cols = colorsT, features = c(Main), ncol = 1)
dev.off() # 關閉輸出圖檔

RidgePlot(marrow_sub_Malignant_End1,cols = colorsT, features = c("TOP2A"), ncol = 1)
# https://datavizpyr.com/ridgeline-plot-with-ggridges-in-r/
# https://www.r-graph-gallery.com/294-basic-ridgeline-plot.html
RidgePlot(marrow_sub_Malignant_End1,cols = colorsT, features = c(Main), ncol = 2,log=TRUE)
RidgePlot(marrow_sub_Malignant_End1,cols = colorsT, features = c(Main), ncol = 2,y.max = 100)

#將Seurat跑出的Cell cycle結果寫入Monocle3的cds檔
cds_sub_Malignant_End1@colData@listData$cell_cycle <- marrow_sub_Malignant_End1@active.ident


##
Maingroup_ciliated_genes <- c("TOP2A")
cds_sub_Malignant_End1 <- cds_sub_Malignant_End1[rowData(cds_sub_Malignant_End1)$gene_short_name %in% Maingroup_ciliated_genes,]

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_Violin_Main_Sub_Malignant_End1_V2.png")) # 設定輸出圖檔
plot_genes_violin(cds_sub_Malignant_End1, group_cells_by="cell_cycle", ncol=2, log_scale = FALSE) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
dev.off() # 關閉輸出圖檔

##
png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_CellCycle_Sub_Malignant_End1_V2.png")) # 設定輸出圖檔
plot_cells(cds_sub_Malignant_End1, color_cells_by="cell_cycle",cell_size=3, label_cell_groups=FALSE, show_trajectory_graph = FALSE) + scale_color_manual(values = colorsT)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_Sub_Malignant_End1_V2.png")) # 設定輸出圖檔
plot_cells(cds_sub_Malignant_End1, genes=c("TOP2A"),cell_size=3, label_cell_groups=FALSE, show_trajectory_graph = FALSE)
dev.off() # 關閉輸出圖檔


plot_cells(cds_sub_Malignant_End1, color_cells_by="cell_cycle",cell_size=2, label_cell_groups=FALSE) + scale_color_manual(values = colorsT)



plot_cells(cds_sub_Malignant_End1, genes=c("TOP2A"),cell_size=2, label_cell_groups=FALSE, show_trajectory_graph = FALSE)

# AFD_lineage_cds_sub_DucT2_TOP2A_Cluster7 <- cds_sub_DucT2_TOP2A_Cluster7[rowData(cds_sub_DucT2_TOP2A_Cluster7)$gene_short_name %in% Main]
# plot_genes_in_pseudotime(AFD_lineage_cds_sub_DucT2_TOP2A_Cluster7,
#                          color_cells_by="cell_cycle",cell_size=2,
#                          min_expr=0.5)+ scale_color_manual(values = colorsT)



########################################  cds_sub_Malignant_End2  ##########################################
#cds_sub_DucT2_TOP2A_Cluster7 <- choose_graph_segments(cds_sub_DucT2 ,clear_cds = FALSE)
cds_sub_Malignant_End2 <- choose_cells(cds ,clear_cds = FALSE)

plot_cells(cds_sub_Malignant_End2, color_cells_by="cell_cycle",cell_size=2, label_cell_groups=FALSE, show_trajectory_graph = FALSE) + scale_color_manual(values = colorsT)


CCdata_sub_Malignant_End2 <- cds_sub_Malignant_End2@assays@data@listData$counts
colnames(CCdata_sub_Malignant_End2) = cds_sub_Malignant_End2@assays@data@listData[["counts"]]@Dimnames[[2]]
rownames(CCdata_sub_Malignant_End2) = cds_sub_Malignant_End2@assays@data@listData[["counts"]]@Dimnames[[1]]


DataCellcycle_sub_Malignant_End2 <- CCdata_sub_Malignant_End2
DataCellcycle_sub_Malignant_End2 <- as(as.matrix(DataCellcycle_sub_Malignant_End2), "dgCMatrix")

marrow_sub_Malignant_End2 <- CreateSeuratObject(counts = DataCellcycle_sub_Malignant_End2)
marrow_sub_Malignant_End2 <- NormalizeData(marrow_sub_Malignant_End2)
marrow_sub_Malignant_End2 <- FindVariableFeatures(marrow_sub_Malignant_End2, selection.method = "vst")
marrow_sub_Malignant_End2 <- ScaleData(marrow_sub_Malignant_End2, features = rownames(marrow_sub_Malignant_End2))
#錯誤: 無法配置大小為 3.9 Gb 的向量
#https://d.cosx.org/d/413001-413001
#memory.limit(15000)


marrow_sub_Malignant_End2 <- RunPCA(marrow_sub_Malignant_End2, features = VariableFeatures(marrow_sub_Malignant_End2), ndims.print = 1:10, nfeatures.print = 25)

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_Malignant_End2_DimHeatmap1.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_Malignant_End2, dims = c(1,2),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_Malignant_End2_DimHeatmap2.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_Malignant_End2, dims = c(3,4),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_Malignant_End2_DimHeatmap3.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_Malignant_End2, dims = c(5,6),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_Malignant_End2_DimHeatmap4.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_Malignant_End2, dims = c(7,8),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_Malignant_End2_DimHeatmap5.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_Malignant_End2, dims = c(9,10),nfeatures = 30)
dev.off() # 關閉輸出圖檔

DimHeatmap(marrow_sub_Malignant_End2, dims = c(1:2),nfeatures = 50)

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_Malignant_End2_DimHeatmapMa.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_Malignant_End2, dims = c(1:18))
dev.off() # 關閉輸出圖檔


################################################




# Also read in a list of cell cycle markers, from Tirosh et al, 2015
# cc.genes <- readLines(con = "~/Downloads/seurat_resources/regev_lab_cell_cycle_genes.txt")
# cc.genes <- read.csv(paste0(PathName,"/Cell cycle/regev_lab_cell_cycle_genesCh.csv"))
# # We can segregate this list into markers of G2/M phase and markers of S
# # phase
# s.genes <- cc.genes[1:43,]
# s.genes2 <- capitalize(tolower(s.genes))
# # https://www.360kuai.com/pc/9ae634cb511f34aee?cota=4&kuai_so=1&tj_url=so_rec&sign=360_7bc3b157
# library(tidyverse)
# G_listCCs.genes <- AnnotationDbi::select(org.Mm.eg.db, keys=s.genes2, columns='ENSEMBL', keytype='SYMBOL')
# # G_listCCs.genes <- getBM(filters= "hgnc_symbol", attributes= c("ensembl_gene_id","hgnc_symbol"),values=s.genes,mart= martCCG,uniqueRows= TRUE,
# #                          quote = "\"")
# G_listCCs.genes2 = G_listCCs.genes[!duplicated(G_listCCs.genes[2]),]
# G_listCCs.genes3 <- na.omit(G_listCCs.genes2[2])
# G_listCCs.genes4 <- G_listCCs.genes3[,1]
# # (Ori)s.genes <- cc.genes[,1:43]
# g2m.genes <- cc.genes[44:97,]
# g2m.genes2 <- capitalize(tolower(g2m.genes))
# G_listCCg2m.genes <- AnnotationDbi::select(org.Mm.eg.db, keys=g2m.genes2, columns='ENSEMBL', keytype='SYMBOL')
# G_listCCg2m.genes2 = G_listCCg2m.genes[!duplicated(G_listCCg2m.genes[2]),]
# G_listCCg2m.genes3 <- na.omit(G_listCCg2m.genes2[2])
# G_listCCg2m.genes4 <-G_listCCg2m.genes3[,1]
# #marrow <- CellCycleScoring(marrow, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# #marrow <- CellCycleScoring(marrow, s.features = G_listCCs.genes2, g2m.features = G_listCCg2m.genes2, set.ident = TRUE)
# 
# 
# 
# marrow_sub_DucT2_TOP2ACenter <- CellCycleScoring(marrow_sub_DucT2_TOP2ACenter, s.features = G_listCCs.genes4, g2m.features = G_listCCg2m.genes4, set.ident = TRUE)
marrow_sub_Malignant_End2 <- CellCycleScoring(marrow_sub_Malignant_End2, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(marrow_sub_Malignant_End2[[]])


##############
colorsT <- c("#FF9912B3", "#32CD3299", "#4169E1B3")
RidgePlot(marrow_sub_Malignant_End2,cols = colorsT, features = c("TOP2A"), ncol = 1)
Main=("TOP2A")

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_RidgePlot_Sub_Malignant_End2_V2.png")) # 設定輸出圖檔

RidgePlot(marrow_sub_Malignant_End2,cols = colorsT, features = c(Main), ncol = 1)
dev.off() # 關閉輸出圖檔

RidgePlot(marrow_sub_Malignant_End2,cols = colorsT, features = c("TOP2A"), ncol = 1)
# https://datavizpyr.com/ridgeline-plot-with-ggridges-in-r/
# https://www.r-graph-gallery.com/294-basic-ridgeline-plot.html
RidgePlot(marrow_sub_Malignant_End2,cols = colorsT, features = c(Main), ncol = 2,log=TRUE)
RidgePlot(marrow_sub_Malignant_End2,cols = colorsT, features = c(Main), ncol = 2,y.max = 100)

#將Seurat跑出的Cell cycle結果寫入Monocle3的cds檔
cds_sub_Malignant_End2@colData@listData$cell_cycle <- marrow_sub_Malignant_End2@active.ident


##
Maingroup_ciliated_genes <- c("TOP2A")
cds_sub_Malignant_End2 <- cds_sub_Malignant_End2[rowData(cds_sub_Malignant_End2)$gene_short_name %in% Maingroup_ciliated_genes,]

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_Violin_Main_Sub_Malignant_End2_V2.png")) # 設定輸出圖檔
plot_genes_violin(cds_sub_Malignant_End2, group_cells_by="cell_cycle", ncol=2, log_scale = FALSE) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
dev.off() # 關閉輸出圖檔

##
png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_CellCycle_Sub_Malignant_End2_V2.png")) # 設定輸出圖檔
plot_cells(cds_sub_Malignant_End2, color_cells_by="cell_cycle",cell_size=3, label_cell_groups=FALSE, show_trajectory_graph = FALSE) + scale_color_manual(values = colorsT)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_Sub_Malignant_End2_V2.png")) # 設定輸出圖檔
plot_cells(cds_sub_Malignant_End2, genes=c("TOP2A"),cell_size=3, label_cell_groups=FALSE, show_trajectory_graph = FALSE)
dev.off() # 關閉輸出圖檔


plot_cells(cds_sub_Malignant_End2, color_cells_by="cell_cycle",cell_size=2, label_cell_groups=FALSE) + scale_color_manual(values = colorsT)



plot_cells(cds_sub_Malignant_End2, genes=c("TOP2A"),cell_size=2, label_cell_groups=FALSE, show_trajectory_graph = FALSE)

# AFD_lineage_cds_sub_DucT2_TOP2A_Cluster7 <- cds_sub_DucT2_TOP2A_Cluster7[rowData(cds_sub_DucT2_TOP2A_Cluster7)$gene_short_name %in% Main]
# plot_genes_in_pseudotime(AFD_lineage_cds_sub_DucT2_TOP2A_Cluster7,
#                          color_cells_by="cell_cycle",cell_size=2,
#                          min_expr=0.5)+ scale_color_manual(values = colorsT)




