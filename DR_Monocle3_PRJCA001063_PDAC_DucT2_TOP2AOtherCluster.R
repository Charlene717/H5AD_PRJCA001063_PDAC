########################  DucT2_TOP2A Other cluster7 ##########################
#cds_sub_DucT2_TOP2A_Cluster7 <- choose_graph_segments(cds_sub_DucT2 ,clear_cds = FALSE)
cds_sub_DucT2_TOP2A_Cluster7 <- choose_cells(cds_subset_K100 ,clear_cds = FALSE)

plot_cells(cds_sub_DucT2_TOP2A_Cluster7, color_cells_by="cell_cycle",cell_size=2, label_cell_groups=FALSE, show_trajectory_graph = FALSE) + scale_color_manual(values = colorsT)


CCdata_sub_DucT2_TOP2A_Cluster7 <- cds_sub_DucT2_TOP2A_Cluster7@assays@data@listData$counts
colnames(CCdata_sub_DucT2_TOP2A_Cluster7) = cds_sub_DucT2_TOP2A_Cluster7@assays@data@listData[["counts"]]@Dimnames[[2]]
rownames(CCdata_sub_DucT2_TOP2A_Cluster7) = cds_sub_DucT2_TOP2A_Cluster7@assays@data@listData[["counts"]]@Dimnames[[1]]


DataCellcycle_sub_DucT2_TOP2A_Cluster7 <- CCdata_sub_DucT2_TOP2A_Cluster7
DataCellcycle_sub_DucT2_TOP2A_Cluster7 <- as(as.matrix(DataCellcycle_sub_DucT2_TOP2A_Cluster7), "dgCMatrix")

marrow_sub_DucT2_TOP2A_Cluster7 <- CreateSeuratObject(counts = DataCellcycle_sub_DucT2_TOP2A_Cluster7)
marrow_sub_DucT2_TOP2A_Cluster7 <- NormalizeData(marrow_sub_DucT2_TOP2A_Cluster7)
marrow_sub_DucT2_TOP2A_Cluster7 <- FindVariableFeatures(marrow_sub_DucT2_TOP2A_Cluster7, selection.method = "vst")
marrow_sub_DucT2_TOP2A_Cluster7 <- ScaleData(marrow_sub_DucT2_TOP2A_Cluster7, features = rownames(marrow_sub_DucT2_TOP2A_Cluster7))
#錯誤: 無法配置大小為 3.9 Gb 的向量
#https://d.cosx.org/d/413001-413001
#memory.limit(15000)


marrow_sub_DucT2_TOP2A_Cluster7 <- RunPCA(marrow_sub_DucT2_TOP2A_Cluster7, features = VariableFeatures(marrow_sub_DucT2_TOP2A_Cluster7), ndims.print = 1:10, nfeatures.print = 25)

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2A_Cluster7_DimHeatmap1.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2A_Cluster7, dims = c(1,2),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2A_Cluster7_DimHeatmap2.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2A_Cluster7, dims = c(3,4),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2A_Cluster7_DimHeatmap3.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2A_Cluster7, dims = c(5,6),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2A_Cluster7_DimHeatmap4.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2A_Cluster7, dims = c(7,8),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2A_Cluster7_DimHeatmap5.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2A_Cluster7, dims = c(9,10),nfeatures = 30)
dev.off() # 關閉輸出圖檔

DimHeatmap(marrow_sub_DucT2_TOP2A_Cluster7, dims = c(1:2),nfeatures = 50)

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2A_Cluster7_DimHeatmapMa.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2A_Cluster7, dims = c(1:18))
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
marrow_sub_DucT2_TOP2A_Cluster7 <- CellCycleScoring(marrow_sub_DucT2_TOP2A_Cluster7, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(marrow_sub_DucT2_TOP2A_Cluster7[[]])


##############
colorsT <- c("#FF9912B3", "#32CD3299", "#4169E1B3")
RidgePlot(marrow_sub_DucT2_TOP2A_Cluster7,cols = colorsT, features = c("TOP2A"), ncol = 1)
Main=("TOP2A")

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_RidgePlot_Sub_DucT2_TOP2A_Cluster7_V2.png")) # 設定輸出圖檔

RidgePlot(marrow_sub_DucT2_TOP2A_Cluster7,cols = colorsT, features = c(Main), ncol = 1)
dev.off() # 關閉輸出圖檔

RidgePlot(marrow_sub_DucT2_TOP2A_Cluster7,cols = colorsT, features = c("TOP2A"), ncol = 1)
# https://datavizpyr.com/ridgeline-plot-with-ggridges-in-r/
# https://www.r-graph-gallery.com/294-basic-ridgeline-plot.html
RidgePlot(marrow_sub_DucT2_TOP2A_Cluster7,cols = colorsT, features = c(Main), ncol = 2,log=TRUE)
RidgePlot(marrow_sub_DucT2_TOP2A_Cluster7,cols = colorsT, features = c(Main), ncol = 2,y.max = 100)

#將Seurat跑出的Cell cycle結果寫入Monocle3的cds檔
cds_sub_DucT2_TOP2A_Cluster7@colData@listData$cell_cycle <- marrow_sub_DucT2_TOP2A_Cluster7@active.ident


##
Maingroup_ciliated_genes <- c("TOP2A")
cds_sub_DucT2_TOP2A_Cluster7 <- cds_sub_DucT2_TOP2A_Cluster7[rowData(cds_sub_DucT2_TOP2A_Cluster7)$gene_short_name %in% Maingroup_ciliated_genes,]

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_Violin_Main_Sub_DucT2_TOP2A_Cluster7_V2.png")) # 設定輸出圖檔
plot_genes_violin(cds_sub_DucT2_TOP2A_Cluster7, group_cells_by="cell_cycle", ncol=2, log_scale = FALSE) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
dev.off() # 關閉輸出圖檔

##
png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_CellCycle_Sub_DucT2_TOP2A_Cluster7_V2.png")) # 設定輸出圖檔
plot_cells(cds_sub_DucT2_TOP2A_Cluster7, color_cells_by="cell_cycle",cell_size=3, label_cell_groups=FALSE, show_trajectory_graph = FALSE) + scale_color_manual(values = colorsT)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_Sub_DucT2_TOP2A_Cluster7_V2.png")) # 設定輸出圖檔
plot_cells(cds_sub_DucT2_TOP2A_Cluster7, genes=c("TOP2A"),cell_size=3, label_cell_groups=FALSE, show_trajectory_graph = FALSE)
dev.off() # 關閉輸出圖檔


plot_cells(cds_sub_DucT2_TOP2A_Cluster7, color_cells_by="cell_cycle",cell_size=2, label_cell_groups=FALSE) + scale_color_manual(values = colorsT)



plot_cells(cds_sub_DucT2_TOP2A_Cluster7, genes=c("TOP2A"),cell_size=2, label_cell_groups=FALSE, show_trajectory_graph = FALSE)

# AFD_lineage_cds_sub_DucT2_TOP2A_Cluster7 <- cds_sub_DucT2_TOP2A_Cluster7[rowData(cds_sub_DucT2_TOP2A_Cluster7)$gene_short_name %in% Main]
# plot_genes_in_pseudotime(AFD_lineage_cds_sub_DucT2_TOP2A_Cluster7,
#                          color_cells_by="cell_cycle",cell_size=2,
#                          min_expr=0.5)+ scale_color_manual(values = colorsT)


########################  DucT2_TOP2A Other cluster10 ##########################
#cds_sub_DucT2_TOP2A_Cluster10 <- choose_graph_segments(cds_sub_DucT2 ,clear_cds = FALSE)
cds_sub_DucT2_TOP2A_Cluster10 <- choose_cells(cds_subset_K100 ,clear_cds = FALSE)

plot_cells(cds_sub_DucT2_TOP2A_Cluster10, color_cells_by="cell_cycle",cell_size=2, label_cell_groups=FALSE, show_trajectory_graph = FALSE) + scale_color_manual(values = colorsT)


CCdata_sub_DucT2_TOP2A_Cluster10 <- cds_sub_DucT2_TOP2A_Cluster10@assays@data@listData$counts
colnames(CCdata_sub_DucT2_TOP2A_Cluster10) = cds_sub_DucT2_TOP2A_Cluster10@assays@data@listData[["counts"]]@Dimnames[[2]]
rownames(CCdata_sub_DucT2_TOP2A_Cluster10) = cds_sub_DucT2_TOP2A_Cluster10@assays@data@listData[["counts"]]@Dimnames[[1]]


DataCellcycle_sub_DucT2_TOP2A_Cluster10 <- CCdata_sub_DucT2_TOP2A_Cluster10
DataCellcycle_sub_DucT2_TOP2A_Cluster10 <- as(as.matrix(DataCellcycle_sub_DucT2_TOP2A_Cluster10), "dgCMatrix")

marrow_sub_DucT2_TOP2A_Cluster10 <- CreateSeuratObject(counts = DataCellcycle_sub_DucT2_TOP2A_Cluster10)
marrow_sub_DucT2_TOP2A_Cluster10 <- NormalizeData(marrow_sub_DucT2_TOP2A_Cluster10)
marrow_sub_DucT2_TOP2A_Cluster10 <- FindVariableFeatures(marrow_sub_DucT2_TOP2A_Cluster10, selection.method = "vst")
marrow_sub_DucT2_TOP2A_Cluster10 <- ScaleData(marrow_sub_DucT2_TOP2A_Cluster10, features = rownames(marrow_sub_DucT2_TOP2A_Cluster10))
#錯誤: 無法配置大小為 3.9 Gb 的向量
#https://d.cosx.org/d/413001-413001
#memory.limit(15000)


marrow_sub_DucT2_TOP2A_Cluster10 <- RunPCA(marrow_sub_DucT2_TOP2A_Cluster10, features = VariableFeatures(marrow_sub_DucT2_TOP2A_Cluster10), ndims.print = 1:10, nfeatures.print = 25)

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2A_Cluster10_DimHeatmap1.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2A_Cluster10, dims = c(1,2),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2A_Cluster10_DimHeatmap2.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2A_Cluster10, dims = c(3,4),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2A_Cluster10_DimHeatmap3.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2A_Cluster10, dims = c(5,6),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2A_Cluster10_DimHeatmap4.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2A_Cluster10, dims = c(7,8),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2A_Cluster10_DimHeatmap5.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2A_Cluster10, dims = c(9,10),nfeatures = 30)
dev.off() # 關閉輸出圖檔

DimHeatmap(marrow_sub_DucT2_TOP2A_Cluster10, dims = c(1:2),nfeatures = 50)

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2A_Cluster10_DimHeatmapMa.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2A_Cluster10, dims = c(1:18))
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
marrow_sub_DucT2_TOP2A_Cluster10 <- CellCycleScoring(marrow_sub_DucT2_TOP2A_Cluster10, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(marrow_sub_DucT2_TOP2A_Cluster10[[]])


##############
colorsT <- c("#FF9912B3", "#32CD3299", "#4169E1B3")
RidgePlot(marrow_sub_DucT2_TOP2A_Cluster10,cols = colorsT, features = c("TOP2A"), ncol = 1)
Main=("TOP2A")

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_RidgePlot_Sub_DucT2_TOP2A_Cluster10_V2.png")) # 設定輸出圖檔

RidgePlot(marrow_sub_DucT2_TOP2A_Cluster10,cols = colorsT, features = c(Main), ncol = 1)
dev.off() # 關閉輸出圖檔

RidgePlot(marrow_sub_DucT2_TOP2A_Cluster10,cols = colorsT, features = c("TOP2A"), ncol = 1)
# https://datavizpyr.com/ridgeline-plot-with-ggridges-in-r/
# https://www.r-graph-gallery.com/294-basic-ridgeline-plot.html
RidgePlot(marrow_sub_DucT2_TOP2A_Cluster10,cols = colorsT, features = c(Main), ncol = 2,log=TRUE)
RidgePlot(marrow_sub_DucT2_TOP2A_Cluster10,cols = colorsT, features = c(Main), ncol = 2,y.max = 100)

#將Seurat跑出的Cell cycle結果寫入Monocle3的cds檔
cds_sub_DucT2_TOP2A_Cluster10@colData@listData$cell_cycle <- marrow_sub_DucT2_TOP2A_Cluster10@active.ident


##
Maingroup_ciliated_genes <- c("TOP2A")
cds_sub_DucT2_TOP2A_Cluster10 <- cds_sub_DucT2_TOP2A_Cluster10[rowData(cds_sub_DucT2_TOP2A_Cluster10)$gene_short_name %in% Maingroup_ciliated_genes,]

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_Violin_Main_Sub_DucT2_TOP2A_Cluster10_V2.png")) # 設定輸出圖檔
plot_genes_violin(cds_sub_DucT2_TOP2A_Cluster10, group_cells_by="cell_cycle", ncol=2, log_scale = FALSE) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
dev.off() # 關閉輸出圖檔

##
png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_CellCycle_Sub_DucT2_TOP2A_Cluster10_V2.png")) # 設定輸出圖檔
plot_cells(cds_sub_DucT2_TOP2A_Cluster10, color_cells_by="cell_cycle",cell_size=3, label_cell_groups=FALSE, show_trajectory_graph = FALSE) + scale_color_manual(values = colorsT)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_Sub_DucT2_TOP2A_Cluster10_V2.png")) # 設定輸出圖檔
plot_cells(cds_sub_DucT2_TOP2A_Cluster10, genes=c("TOP2A"),cell_size=3, label_cell_groups=FALSE, show_trajectory_graph = FALSE)
dev.off() # 關閉輸出圖檔


plot_cells(cds_sub_DucT2_TOP2A_Cluster10, color_cells_by="cell_cycle",cell_size=2, label_cell_groups=FALSE) + scale_color_manual(values = colorsT)



plot_cells(cds_sub_DucT2_TOP2A_Cluster10, genes=c("TOP2A"),cell_size=2, label_cell_groups=FALSE, show_trajectory_graph = FALSE)

# AFD_lineage_cds_sub_DucT2_TOP2A_Cluster7 <- cds_sub_DucT2_TOP2A_Cluster7[rowData(cds_sub_DucT2_TOP2A_Cluster7)$gene_short_name %in% Main]
# plot_genes_in_pseudotime(AFD_lineage_cds_sub_DucT2_TOP2A_Cluster7,
#                          color_cells_by="cell_cycle",cell_size=2,
#                          min_expr=0.5)+ scale_color_manual(values = colorsT)






########################  DucT2_TOP2A Other cluster2 ##########################
#cds_sub_DucT2_TOP2A_Cluster2 <- choose_graph_segments(cds_sub_DucT2 ,clear_cds = FALSE)
cds_sub_DucT2_TOP2A_Cluster2 <- choose_cells(cds_subset_K100 ,clear_cds = FALSE)

plot_cells(cds_sub_DucT2_TOP2A_Cluster2, color_cells_by="cell_cycle",cell_size=2, label_cell_groups=FALSE, show_trajectory_graph = FALSE) + scale_color_manual(values = colorsT)


CCdata_sub_DucT2_TOP2A_Cluster2 <- cds_sub_DucT2_TOP2A_Cluster2@assays@data@listData$counts
colnames(CCdata_sub_DucT2_TOP2A_Cluster2) = cds_sub_DucT2_TOP2A_Cluster2@assays@data@listData[["counts"]]@Dimnames[[2]]
rownames(CCdata_sub_DucT2_TOP2A_Cluster2) = cds_sub_DucT2_TOP2A_Cluster2@assays@data@listData[["counts"]]@Dimnames[[1]]


DataCellcycle_sub_DucT2_TOP2A_Cluster2 <- CCdata_sub_DucT2_TOP2A_Cluster2
DataCellcycle_sub_DucT2_TOP2A_Cluster2 <- as(as.matrix(DataCellcycle_sub_DucT2_TOP2A_Cluster2), "dgCMatrix")

marrow_sub_DucT2_TOP2A_Cluster2 <- CreateSeuratObject(counts = DataCellcycle_sub_DucT2_TOP2A_Cluster2)
marrow_sub_DucT2_TOP2A_Cluster2 <- NormalizeData(marrow_sub_DucT2_TOP2A_Cluster2)
marrow_sub_DucT2_TOP2A_Cluster2 <- FindVariableFeatures(marrow_sub_DucT2_TOP2A_Cluster2, selection.method = "vst")
marrow_sub_DucT2_TOP2A_Cluster2 <- ScaleData(marrow_sub_DucT2_TOP2A_Cluster2, features = rownames(marrow_sub_DucT2_TOP2A_Cluster2))
#錯誤: 無法配置大小為 3.9 Gb 的向量
#https://d.cosx.org/d/413001-413001
#memory.limit(15000)


marrow_sub_DucT2_TOP2A_Cluster2 <- RunPCA(marrow_sub_DucT2_TOP2A_Cluster2, features = VariableFeatures(marrow_sub_DucT2_TOP2A_Cluster2), ndims.print = 1:10, nfeatures.print = 25)

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2A_Cluster2_DimHeatmap1.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2A_Cluster2, dims = c(1,2),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2A_Cluster2_DimHeatmap2.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2A_Cluster2, dims = c(3,4),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2A_Cluster2_DimHeatmap3.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2A_Cluster2, dims = c(5,6),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2A_Cluster2_DimHeatmap4.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2A_Cluster2, dims = c(7,8),nfeatures = 30)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2A_Cluster2_DimHeatmap5.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2A_Cluster2, dims = c(9,10),nfeatures = 30)
dev.off() # 關閉輸出圖檔

DimHeatmap(marrow_sub_DucT2_TOP2A_Cluster2, dims = c(1:2),nfeatures = 50)

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2A_Cluster2_DimHeatmapMa.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2A_Cluster2, dims = c(1:18))
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
marrow_sub_DucT2_TOP2A_Cluster2 <- CellCycleScoring(marrow_sub_DucT2_TOP2A_Cluster2, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(marrow_sub_DucT2_TOP2A_Cluster2[[]])


##############
colorsT <- c("#FF9912B3", "#32CD3299", "#4169E1B3")
RidgePlot(marrow_sub_DucT2_TOP2A_Cluster2,cols = colorsT, features = c("TOP2A"), ncol = 1)
Main=("TOP2A")

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_RidgePlot_Sub_DucT2_TOP2A_Cluster2_V2.png")) # 設定輸出圖檔

RidgePlot(marrow_sub_DucT2_TOP2A_Cluster2,cols = colorsT, features = c(Main), ncol = 1)
dev.off() # 關閉輸出圖檔

RidgePlot(marrow_sub_DucT2_TOP2A_Cluster2,cols = colorsT, features = c("TOP2A"), ncol = 1)
# https://datavizpyr.com/ridgeline-plot-with-ggridges-in-r/
# https://www.r-graph-gallery.com/294-basic-ridgeline-plot.html
RidgePlot(marrow_sub_DucT2_TOP2A_Cluster2,cols = colorsT, features = c(Main), ncol = 2,log=TRUE)
RidgePlot(marrow_sub_DucT2_TOP2A_Cluster2,cols = colorsT, features = c(Main), ncol = 2,y.max = 100)

#將Seurat跑出的Cell cycle結果寫入Monocle3的cds檔
cds_sub_DucT2_TOP2A_Cluster2@colData@listData$cell_cycle <- marrow_sub_DucT2_TOP2A_Cluster2@active.ident


##
Maingroup_ciliated_genes <- c("TOP2A")
cds_sub_DucT2_TOP2A_Cluster2 <- cds_sub_DucT2_TOP2A_Cluster2[rowData(cds_sub_DucT2_TOP2A_Cluster2)$gene_short_name %in% Maingroup_ciliated_genes,]

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_Violin_Main_Sub_DucT2_TOP2A_Cluster2_V2.png")) # 設定輸出圖檔
plot_genes_violin(cds_sub_DucT2_TOP2A_Cluster2, group_cells_by="cell_cycle", ncol=2, log_scale = FALSE) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
dev.off() # 關閉輸出圖檔

##
png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_CellCycle_Sub_DucT2_TOP2A_Cluster2_V2.png")) # 設定輸出圖檔
plot_cells(cds_sub_DucT2_TOP2A_Cluster2, color_cells_by="cell_cycle",cell_size=3, label_cell_groups=FALSE, show_trajectory_graph = FALSE) + scale_color_manual(values = colorsT)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_Sub_DucT2_TOP2A_Cluster2_V2.png")) # 設定輸出圖檔
plot_cells(cds_sub_DucT2_TOP2A_Cluster2, genes=c("TOP2A"),cell_size=3, label_cell_groups=FALSE, show_trajectory_graph = FALSE)
dev.off() # 關閉輸出圖檔


plot_cells(cds_sub_DucT2_TOP2A_Cluster2, color_cells_by="cell_cycle",cell_size=2, label_cell_groups=FALSE) + scale_color_manual(values = colorsT)



plot_cells(cds_sub_DucT2_TOP2A_Cluster2, genes=c("TOP2A"),cell_size=2, label_cell_groups=FALSE, show_trajectory_graph = FALSE)

# AFD_lineage_cds_sub_DucT2_TOP2A_Cluster7 <- cds_sub_DucT2_TOP2A_Cluster7[rowData(cds_sub_DucT2_TOP2A_Cluster7)$gene_short_name %in% Main]
# plot_genes_in_pseudotime(AFD_lineage_cds_sub_DucT2_TOP2A_Cluster7,
#                          color_cells_by="cell_cycle",cell_size=2,
#                          min_expr=0.5)+ scale_color_manual(values = colorsT)
