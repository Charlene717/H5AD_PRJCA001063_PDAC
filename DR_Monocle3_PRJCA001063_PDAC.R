# https://rdrr.io/github/theislab/zellkonverter/man/readH5AD.html
# https://rdrr.io/cran/Seurat/man/h5ad.html
# https://www.rdocumentation.org/packages/Seurat/versions/3.1.4/topics/ReadH5AD

#############
rm(list = ls()) #Delete variable

############# Library list #############

library(SummarizedExperiment)
library(Seurat)
library(SeuratDisk)
library(stringr)
library(SeuratWrappers)
library(monocle3)
library(AnnotationDbi)
library(org.Mm.eg.db)
library('org.Hs.eg.db')
library(Hmisc)
library(tidyverse)
library(garnett)
library(cicero)

############# Library list #############



############# Read file settings #############

PathName = setwd(getwd())
RVersion = "20210501V1"
dir.create(paste0(PathName,"/",RVersion))

# Marker gene file
Marker_file_Name <- c("NAKAMURA_METASTASIS_MODEL_M18483")
Marker_file <- paste0(PathName,"/marker_file_",Marker_file_Name,".txt")
Marker_List <- read.table(Marker_file,header=F,sep= c(","),stringsAsFactors = FALSE, fill = TRUE)
library(stringr)
Marker_List_1 <- Marker_List[2,1]
Marker_List_2 <- str_replace_all(Marker_List_1,"expressed: ","")
Marker_List <- str_trim(Marker_List[2,-1], side = c("both"))
Marker_List <- c(Marker_List_2,Marker_List)

############# Read file settings #############


######################################## Gene list of interest ########################################
Main = c("TOP2A")
Main_Group = c("TOP2A","TOP2B","TP53","CCNE1")
Main_Group2 = c("KRAS","EXO1","NSUN2","MUC1","AMBP","FXYD2")
EMT_Meta = c("ANLN","APLP2","CD63","CDH2","CLIC4","CTSB","CX3CR1","DSG2","EDNRB")
candidates14 = c("BRIP1","KIF23","TOP2A","FOSL1","FAM25A","ANLN","NCAPH","KRT9","MCM4","CKAP2L","CENPE","RACGAP1","DTL","RAD51AP1")

DREAM_complex= c("RBL2","E2F4","E2F5","TFDP1","TFDP2")
Regulators= c("TP53","YBX1","E2F1")
#######################################################################################################

# https://rdrr.io/github/satijalab/seurat/man/AddModuleScore.html


library(SummarizedExperiment)
library(Seurat)
# https://satijalab.org/seurat/install.html # https://rdrr.io/rforge/matrixStats/man/matrixStats-package.html
# https://blog.csdn.net/zengwanqin/article/details/114895417
# https://www.youtube.com/watch?v=Eucn_BJ8EJI

library(SeuratDisk)
# https://github.com/mojaveazure/seurat-disk
# library(cellexalvrR)
# # https://cellexalvr.med.lu.se/cellexalvrr-vignette

####20210108 https://github.com/satijalab/seurat/issues/3414
library(SeuratDisk)

Convert(paste0(PathName,"/StdWf1_PRJCA001063_CRC_besca2.raw.h5ad"), "PRJCA001063.h5seurat")
# This creates a copy of this .h5ad object reformatted into .h5seurat inside the example_dir directory

# This .d5seurat object can then be read in manually
seuratObject <- LoadH5Seurat(paste0(PathName,"/PRJCA001063.h5seurat"))
# 

# seurat Object to Monocle3
# https://satijalab.org/signac/articles/monocle.html
# as.cell_data_set: Convert objects to Monocle3 'cell_data_set' objects
# https://rdrr.io/github/satijalab/seurat-wrappers/man/as.cell_data_set.html
library(SeuratWrappers)
erythroid.cds <- as.cell_data_set(seuratObject)


# run Monocle3
library(monocle3)
cds <- erythroid.cds
# issues with cds object in monocle3 #54
# https://github.com/satijalab/seurat-wrappers/issues/54
cds <- estimate_size_factors(cds)
cds <- preprocess_cds(cds, num_dim = 100)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds)

plot_cells(cds)

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","Cell_type.png")) # 設定輸出圖檔
plot_cells(cds, color_cells_by="Cell_type", label_cell_groups=FALSE, show_trajectory_graph = FALSE)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","Type.png")) # 設定輸出圖檔
plot_cells(cds, color_cells_by="Type", label_cell_groups=FALSE, show_trajectory_graph = FALSE)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","Patient.png")) # 設定輸出圖檔
plot_cells(cds, color_cells_by="Patient", label_cell_groups=FALSE, show_trajectory_graph = FALSE)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CONDITION.png")) # 設定輸出圖檔
plot_cells(cds, color_cells_by="CONDITION", label_cell_groups=FALSE, show_trajectory_graph = FALSE)
dev.off() # 關閉輸出圖檔

plot_cells(cds, genes=c("TOP2A","TOP2B","TP53")) #error
cds@rowRanges@elementMetadata@listData$gene_short_name <- cds@assays@data@listData[["counts"]]@Dimnames[[1]]

plot_cells(cds, genes=c("TOP2A","TOP2B","TP53","CCNE1")) #ok

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_",Main,".png")) # 設定輸出圖檔
plot_cells(cds, genes=c(Main),cell_size=1,label_cell_groups = FALSE, show_trajectory_graph = FALSE)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_","Main_Group",".png")) # 設定輸出圖檔
plot_cells(cds, genes=c(Main_Group),cell_size=0.5,label_cell_groups = FALSE, show_trajectory_graph = FALSE)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_","Main_Group2",".png")) # 設定輸出圖檔
plot_cells(cds, genes=c(Main_Group2),cell_size=0.5,label_cell_groups = FALSE, show_trajectory_graph = FALSE)
dev.off() # 關閉輸出圖檔

# plot_cells(cds, genes=c("TOP2A"))
# plot_cells(cds, genes=c("TOP2A"), label_cell_groups=FALSE, show_trajectory_graph = FALSE)
# plot_cells(cds, genes=c("TOP2B"), label_cell_groups=FALSE, show_trajectory_graph = FALSE)
# plot_cells(cds, genes=c("EXO1"), label_cell_groups=FALSE, show_trajectory_graph = FALSE)
# plot_cells(cds, genes=c("CDK2"))
# plot_cells(cds, genes=c("CDK1"))
# plot_cells(cds, genes=c("CCNE1"))
# 
# #
# plot_cells(cds, genes=c("KRAS"))
# plot_cells(cds, genes=c("TP53"))
# plot_cells(cds, genes=c("NSUN2"))
# #
# plot_cells(cds, genes=c("MUC1"))
# plot_cells(cds, genes=c("AMBP"))
# plot_cells(cds, genes=c("FXYD2"))
# 
# #
# plot_cells(cds, genes=c("BRIP1","KIF23","FOSL1","FAM25A")) 
# plot_cells(cds, genes=c("ANLN","NCAPH","KRT9","MCM4")) 
# plot_cells(cds, genes=c("CKAP2L","CENPE","RACGAP1","DTL")) 
# plot_cells(cds, genes=c("RAD51AP1"))

cds <- cluster_cells(cds)
png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_","cluster_par",".png")) # 設定輸出圖檔
plot_cells(cds, color_cells_by = "partition", label_cell_groups=FALSE, show_trajectory_graph = FALSE)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_","cluster_clu",".png")) # 設定輸出圖檔
plot_cells(cds, color_cells_by = "cluster", label_cell_groups=FALSE, show_trajectory_graph = FALSE)
dev.off() # 關閉輸出圖檔

#########################################Cell cycle#################################################################
### Cell cycle

# Humans contain two distinct types of cyclin A: A1, the embryonic-specific form, and A2, the somatic form. Cyclin A1 is prevalently expressed during meiosis and early on in embryogenesis. Cyclin A2 is expressed in dividing somatic cells.
# https://en.wikipedia.org/wiki/Cyclin_A
# Yam CH, Fung TK, Poon RY (August 2002). "Cyclin A in cell cycle control and cancer". Cell. Mol. Life Sci. 59 (8): 1317–26. doi:10.1007/s00018-002-8510-y. PMID 12363035.

### Assign Cell-Cycle Scores
library(Seurat)

# https://github.com/rstudio/rstudio/issues/4741
library(SummarizedExperiment) 

library(AnnotationDbi)
library(org.Mm.eg.db)
library('org.Hs.eg.db')
library(Hmisc)

CCdata <- cds@assays@data@listData$counts
colnames(CCdata) = cds@assays@data@listData[["counts"]]@Dimnames[[2]]
rownames(CCdata) = cds@assays@data@listData[["counts"]]@Dimnames[[1]]


DataCellcycle <- CCdata
#rownames(DataCellcycle) <- make.names(DataCellcycle[,1], unique = TRUE)
#DataCellcycle <- DataCellcycle[1:length(data[,1]), 2:length(data[1,])]
DataCellcycle <- as(as.matrix(DataCellcycle), "dgCMatrix")

marrow <- CreateSeuratObject(counts = DataCellcycle)
marrow <- NormalizeData(marrow)
marrow <- FindVariableFeatures(marrow, selection.method = "vst")
marrow <- ScaleData(marrow, features = rownames(marrow))
#錯誤: 無法配置大小為 3.9 Gb 的向量
#https://d.cosx.org/d/413001-413001
#memory.limit(15000)

marrow@assays[["RNA"]]@counts@Dimnames[[1]] <- cds@assays@data@listData[["counts"]]@Dimnames[[1]]
marrow@assays[["RNA"]]@data@Dimnames[[1]] <- cds@assays@data@listData[["counts"]]@Dimnames[[1]]
marrow@commands[["ScaleData.RNA"]]@params[["features"]] <- cds@assays@data@listData[["counts"]]@Dimnames[[1]]

marrow <- RunPCA(marrow, features = VariableFeatures(marrow), ndims.print = 1:10, nfeatures.print = 10)

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DimHeatmap.png")) # 設定輸出圖檔
DimHeatmap(marrow, dims = c(1: 8))
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DimHeatmap2.png")) # 設定輸出圖檔
DimHeatmap(marrow, dims = c(1, 2))
dev.off() # 關閉輸出圖檔

# Also read in a list of cell cycle markers, from Tirosh et al, 2015
# cc.genes <- readLines(con = "~/Downloads/seurat_resources/regev_lab_cell_cycle_genes.txt")
cc.genes <- read.csv(paste0(PathName,"/Cell cycle/regev_lab_cell_cycle_genesCh.csv"))
# We can segregate this list into markers of G2/M phase and markers of S
# phase
s.genes <- cc.genes[1:43,]
s.genes2 <- capitalize(tolower(s.genes))
# https://www.360kuai.com/pc/9ae634cb511f34aee?cota=4&kuai_so=1&tj_url=so_rec&sign=360_7bc3b157
library(tidyverse)
G_listCCs.genes <- AnnotationDbi::select(org.Mm.eg.db, keys=s.genes2, columns='ENSEMBL', keytype='SYMBOL')
# G_listCCs.genes <- getBM(filters= "hgnc_symbol", attributes= c("ensembl_gene_id","hgnc_symbol"),values=s.genes,mart= martCCG,uniqueRows= TRUE,
#                          quote = "\"")
G_listCCs.genes2 = G_listCCs.genes[!duplicated(G_listCCs.genes[2]),]
G_listCCs.genes3 <- na.omit(G_listCCs.genes2[2])
G_listCCs.genes4 <- G_listCCs.genes3[,1]
# (Ori)s.genes <- cc.genes[,1:43]
g2m.genes <- cc.genes[44:97,]
g2m.genes2 <- capitalize(tolower(g2m.genes))
G_listCCg2m.genes <- AnnotationDbi::select(org.Mm.eg.db, keys=g2m.genes2, columns='ENSEMBL', keytype='SYMBOL')
G_listCCg2m.genes2 = G_listCCg2m.genes[!duplicated(G_listCCg2m.genes[2]),]
G_listCCg2m.genes3 <- na.omit(G_listCCg2m.genes2[2])
G_listCCg2m.genes4 <-G_listCCg2m.genes3[,1]
#marrow <- CellCycleScoring(marrow, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
#marrow <- CellCycleScoring(marrow, s.features = G_listCCs.genes2, g2m.features = G_listCCg2m.genes2, set.ident = TRUE)

marrow <- CellCycleScoring(marrow, s.features = G_listCCs.genes4, g2m.features = G_listCCg2m.genes4, set.ident = TRUE)
marrow <- CellCycleScoring(marrow, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(marrow[[]])

############################################################################################

colorsT <- c("#FF9912B3", "#32CD3299", "#4169E1B3") 

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_RidgePlot_Main.png")) # 設定輸出圖檔

RidgePlot(marrow,cols = colorsT, features = c(Main), ncol = 1) 

# RidgePlot(marrow,cols = colorsT, features = c(Main), ncol = 2,idents=c("G1","S","G2M")) 
# RidgePlot(marrow, features = c(MainGroup[,1]), ncol = 2,idents=c("G1","S","G2M"),sort= 'increasing' ,same.y.lims=TRUE) 
dev.off() # 關閉輸出圖檔

# https://datavizpyr.com/ridgeline-plot-with-ggridges-in-r/
# https://www.r-graph-gallery.com/294-basic-ridgeline-plot.html
RidgePlot(marrow,cols = colorsT, features = c(Main), ncol = 2,log=TRUE) 
RidgePlot(marrow,cols = colorsT, features = c(Main), ncol = 2,y.max = 100) 

#將Seurat跑出的Cell cycle結果寫入Monocle3的cds檔

#cell_cycle <- marrow@active.ident
cds@colData@listData$cell_cycle <- marrow@active.ident
#cds@colData@listData$cell_cycle <- marrow@meta.data[["Phase"]]
#TTT <- marrow@active.ident

##
Maingroup_ciliated_genes <- c(Main)
cds_marrow_subset <- cds[rowData(cds)$gene_short_name %in% Maingroup_ciliated_genes,]

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_Violin_Main.png")) # 設定輸出圖檔
plot_genes_violin(cds_marrow_subset, group_cells_by="cell_cycle", ncol=2, log_scale = FALSE) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
dev.off() # 關閉輸出圖檔

##
png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_CellCycle.png")) # 設定輸出圖檔

colorsT <- c("#FF9912B3", "#32CD3299", "#4169E1B3", "#4169E1B4") 

#plot_cells(cds  , color_cells_by="cell_cycle",cell_size=2, label_cell_groups=FALSE) + scale_color_manual(values = colorsT)
plot_cells(cds  , color_cells_by="cell_cycle",cell_size=1, label_cell_groups=FALSE) + scale_color_manual(values = colorsT)

#plot_cells(cds, color_cells_by="cell_cycle",cell_size=2, label_cell_groups=FALSE)

dev.off() # 關閉輸出圖檔
plot_cells(cds  , color_cells_by="cell_cycle", label_cell_groups=FALSE) + scale_color_manual(values = colorsT)


# ##################
# library(garnett)
# ## Install the monocle3 branch of garnett
# # BiocManager::install(c("org.Mm.eg.db", "org.Hs.eg.db"))
# # devtools::install_github("cole-trapnell-lab/garnett", ref="monocle3")
# # Install problem
# # https://stackoverflow.com/questions/42807247/installing-package-cannot-open-file-permission-denied
# 
# colData(cds)$garnett_cluster <- clusters(cds)
# Human_classifier <- train_cell_classifier(cds = cds,
#                                           marker_file = Marker_file,   # Import the marker_file
#                                           db=org.Hs.eg.db::org.Hs.eg.db,
#                                           cds_gene_id_type = "SYMBOL",
#                                           #num_unknown = 2215,
#                                           #max_training_samples = 10000,
#                                           marker_file_gene_id_type = "SYMBOL",
#                                           cores=8)
# 
# 
# cds_subset <- classify_cells(cds, Human_classifier,
#                              db = org.Hs.eg.db::org.Hs.eg.db,
#                              cluster_extend = TRUE,
#                              cds_gene_id_type = "SYMBOL")
# 
# 
# 
# png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_",Marker_file_Name,".png")) # 設定輸出圖檔
# plot_cells(cds,
#            group_cells_by="cluster",
#            cell_size=1.5,
#            color_cells_by="cluster_ext_type", show_trajectory_graph = FALSE)
# dev.off() # 關閉輸出圖檔
# 
# png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_",Marker_file_Name,"_2.png")) # 設定輸出圖檔
# plot_cells(cds,
#            group_cells_by="cluster",
#            cell_size=1.5,
#            color_cells_by="cluster_ext_type",
#            label_cell_groups=FALSE, show_trajectory_graph = FALSE)
# dev.off() # 關閉輸出圖檔


######################################  cds_subset ########################################
########################  DucT2 ##########################
cds_sub_DucT2 <- choose_cells(cds)
#cds_subset <- reduce_dimension(cds_subset)
plot_cells(cds_sub_DucT2, label_cell_groups=FALSE, show_trajectory_graph = FALSE)

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_","_sub_DucT2_cluster_clu_ori",".png")) # 設定輸出圖檔
plot_cells(cds_sub_DucT2, color_cells_by = "cluster", label_cell_groups=FALSE, show_trajectory_graph = FALSE)
dev.off() # 關閉輸出圖檔

cds_subset_K20 <- cluster_cells(cds_sub_DucT2,k = 20, resolution=1e-5)
png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_","_sub_DucT2_cluster_cluster_clu_K20",".png")) # 設定輸出圖檔
plot_cells(cds_subset_K20, color_cells_by = "cluster", label_cell_groups=FALSE, show_trajectory_graph = FALSE)
dev.off() # 關閉輸出圖檔

cds_subset_K50 <- cluster_cells(cds_sub_DucT2,k = 50, resolution=1e-5)
png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_","_sub_DucT2_cluster_clu_K50",".png")) # 設定輸出圖檔
plot_cells(cds_subset_K50, color_cells_by = "cluster", label_cell_groups=FALSE, show_trajectory_graph = FALSE)
dev.off() # 關閉輸出圖檔

cds_subset_K10 <- cluster_cells(cds_sub_DucT2,k = 10, resolution=1e-5)
png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_","_sub_DucT2_cluster_clu_K10",".png")) # 設定輸出圖檔
plot_cells(cds_subset_K10, color_cells_by = "cluster", label_cell_groups=FALSE, show_trajectory_graph = FALSE)
dev.off() # 關閉輸出圖檔

cds_subset_K100 <- cluster_cells(cds_sub_DucT2,k = 100, resolution=1e-5)
png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_","_sub_DucT2_cluster_clu_K100",".png")) # 設定輸出圖檔
plot_cells(cds_subset_K100, color_cells_by = "cluster", label_cell_groups=FALSE, show_trajectory_graph = FALSE)
dev.off() # 關閉輸出圖檔

cds_subset_K100 <- cluster_cells(cds_sub_DucT2,k = 100, resolution=1e-5)
png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_","_sub_DucT2_cluster_clu_K100_2",".png")) # 設定輸出圖檔
plot_cells(cds_subset_K100, color_cells_by = "cluster",group_label_size = 5, show_trajectory_graph = FALSE)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_","_sub_DucT2_cluster_clu_K100_",Main,".png")) # 設定輸出圖檔
plot_cells(cds_subset_K100, genes=c(Main),cell_size=1,label_cell_groups = FALSE, show_trajectory_graph = FALSE)
dev.off() # 關閉輸出圖檔
 png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_","_sub_DucT2_cluster_clu_K100_TP53.png")) # 設定輸出圖檔
 plot_cells(cds_subset_K100, genes=c("TP53"),cell_size=1,label_cell_groups = FALSE, show_trajectory_graph = FALSE)
 dev.off() # 關閉輸出圖檔
 png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_","_sub_DucT2_cluster_clu_K100_PTK2.png")) # 設定輸出圖檔
 plot_cells(cds_subset_K100, genes=c("PTK2"),cell_size=1,label_cell_groups = FALSE, show_trajectory_graph = FALSE)
 dev.off() # 關閉輸出圖檔
 png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_","_sub_DucT2_cluster_clu_K100_KRAS.png")) # 設定輸出圖檔
 plot_cells(cds_subset_K100, genes=c("KRAS"),cell_size=1,label_cell_groups = FALSE, show_trajectory_graph = FALSE)
 dev.off() # 關閉輸出圖檔
 png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_","_sub_DucT2_cluster_clu_K100_BRIP1.png")) # 設定輸出圖檔
 plot_cells(cds_subset_K100, genes=c("BRIP1"),cell_size=1,label_cell_groups = FALSE, show_trajectory_graph = FALSE)
 dev.off() # 關閉輸出圖檔
 png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_","_sub_DucT2_cluster_clu_K100_H2AX.png")) # 設定輸出圖檔
 plot_cells(cds_subset_K100, genes=c("H2AX"),cell_size=1,label_cell_groups = FALSE, show_trajectory_graph = FALSE)
 dev.off() # 關閉輸出圖檔
 png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_","_sub_DucT2_cluster_clu_K100_BRCA2.png")) # 設定輸出圖檔
 plot_cells(cds_subset_K100, genes=c("BRCA2"),cell_size=1,label_cell_groups = FALSE, show_trajectory_graph = FALSE)
 dev.off() # 關閉輸出圖檔
 png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_","_sub_DucT2_cluster_clu_K100_HMMR.png")) # 設定輸出圖檔
 plot_cells(cds_subset_K100, genes=c("HMMR"),cell_size=1,label_cell_groups = FALSE, show_trajectory_graph = FALSE)
 dev.off() # 關閉輸出圖檔
 png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_","_sub_DucT2_cluster_clu_K100_UBE2S.png")) # 設定輸出圖檔
 plot_cells(cds_subset_K100, genes=c("UBE2S"),cell_size=1,label_cell_groups = FALSE, show_trajectory_graph = FALSE)
 dev.off() # 關閉輸出圖檔
 png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_","_sub_DucT2_cluster_clu_K100_TOP2B.png")) # 設定輸出圖檔
 plot_cells(cds_subset_K100, genes=c("TOP2B"),cell_size=1,label_cell_groups = FALSE, show_trajectory_graph = FALSE)
 dev.off() # 關閉輸出圖檔
 png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_","_sub_DucT2_cluster_clu_K100_NDRG1.png")) # 設定輸出圖檔
 plot_cells(cds_subset_K100, genes=c("NDRG1"),cell_size=1,label_cell_groups = FALSE, show_trajectory_graph = FALSE)
 dev.off() # 關閉輸出圖檔
 png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_","_sub_DucT2_cluster_clu_K100_HLA-A.png")) # 設定輸出圖檔
 plot_cells(cds_subset_K100, genes=c("HLA-A"),cell_size=1,label_cell_groups = FALSE, show_trajectory_graph = FALSE)
 dev.off() # 關閉輸出圖檔
 png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_","_sub_DucT2_cluster_clu_K100_ADM.png")) # 設定輸出圖檔
 plot_cells(cds_subset_K100, genes=c("ADM"),cell_size=1,label_cell_groups = FALSE, show_trajectory_graph = FALSE)
 dev.off() # 關閉輸出圖檔
 png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_","_sub_DucT2_cluster_clu_K100_PAG1.png")) # 設定輸出圖檔
 plot_cells(cds_subset_K100, genes=c("PAG1"),cell_size=1,label_cell_groups = FALSE, show_trajectory_graph = FALSE)
 dev.off() # 關閉輸出圖檔
 
 
 png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_","_sub_DucT2_cluster_clu_K100_AGR2.png")) # 設定輸出圖檔
 plot_cells(cds_subset_K100, genes=c("AGR2"),cell_size=1,label_cell_groups = FALSE, show_trajectory_graph = FALSE)
 dev.off() # 關閉輸出圖檔
 
 png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_","_sub_DucT2_cluster_clu_K100_MUC1.png")) # 設定輸出圖檔
 plot_cells(cds_subset_K100, genes=c("MUC1"),cell_size=1,label_cell_groups = FALSE, show_trajectory_graph = FALSE)
 dev.off() # 關閉輸出圖檔
 png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_","_sub_DucT2_cluster_clu_K100_MUC6.png")) # 設定輸出圖檔
 plot_cells(cds_subset_K100, genes=c("MUC6"),cell_size=1,label_cell_groups = FALSE, show_trajectory_graph = FALSE)
 dev.off() # 關閉輸出圖檔
 png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_","_sub_DucT2_cluster_clu_K100_MUC4.png")) # 設定輸出圖檔
 plot_cells(cds_subset_K100, genes=c("MUC4"),cell_size=1,label_cell_groups = FALSE, show_trajectory_graph = FALSE)
 dev.off() # 關閉輸出圖檔
 png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_","_sub_DucT2_cluster_clu_K100_MUC5AC.png")) # 設定輸出圖檔
 plot_cells(cds_subset_K100, genes=c("MUC5AC"),cell_size=1,label_cell_groups = FALSE, show_trajectory_graph = FALSE)
 dev.off() # 關閉輸出圖檔
 
 png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_","_sub_DucT2_cluster_clu_K100_FXYD3.png")) # 設定輸出圖檔
 plot_cells(cds_subset_K100, genes=c("FXYD3"),cell_size=1,label_cell_groups = FALSE, show_trajectory_graph = FALSE)
 dev.off() # 關閉輸出圖檔
 
 png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_","_sub_DucT2_cluster_clu_K100_candidates14.png")) # 設定輸出圖檔
 plot_cells(cds_subset_K100, genes= candidates14,cell_size=1,label_cell_groups = FALSE, show_trajectory_graph = FALSE)
 dev.off() # 關閉輸出圖檔
 
 
##  Find marker genes expressed by each cluster
marker_test_res_DucT2 <- top_markers(cds_subset_K100, group_cells_by="cluster", 
                               reference_cells=1000, cores=8)

top_specific_markers_DucT2 <- marker_test_res_DucT2 %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(100, pseudo_R2)

top_specific_marker_ids_DucT2  <- unique(top_specific_markers_DucT2  %>% pull(gene_id))


plot_genes_by_group(cds_subset_K100,
                    top_specific_marker_ids_DucT2,
                    group_cells_by="cluster",
                    ordering_type="maximal_on_diag",
                    max.size=3)

top_specific_markers_DucT2 <- data.frame(top_specific_markers_DucT2)

top_specific_markers_DucT2_Sub1 <- top_specific_markers_DucT2[top_specific_markers_DucT2$cell_group =="1",]
top_specific_markers_DucT2_Sub2 <- top_specific_markers_DucT2[top_specific_markers_DucT2$cell_group =="2",]
top_specific_markers_DucT2_Sub3 <- top_specific_markers_DucT2[top_specific_markers_DucT2$cell_group =="3",]
top_specific_markers_DucT2_Sub4 <- top_specific_markers_DucT2[top_specific_markers_DucT2$cell_group =="4",]
top_specific_markers_DucT2_Sub5 <- top_specific_markers_DucT2[top_specific_markers_DucT2$cell_group =="5",]
top_specific_markers_DucT2_Sub6 <- top_specific_markers_DucT2[top_specific_markers_DucT2$cell_group =="6",]
top_specific_markers_DucT2_Sub7 <- top_specific_markers_DucT2[top_specific_markers_DucT2$cell_group =="7",]
top_specific_markers_DucT2_Sub8 <- top_specific_markers_DucT2[top_specific_markers_DucT2$cell_group =="8",]
top_specific_markers_DucT2_Sub9 <- top_specific_markers_DucT2[top_specific_markers_DucT2$cell_group =="9",]
top_specific_markers_DucT2_Sub10 <- top_specific_markers_DucT2[top_specific_markers_DucT2$cell_group =="10",]
top_specific_markers_DucT2_Sub11 <- top_specific_markers_DucT2[top_specific_markers_DucT2$cell_group =="11",]
top_specific_markers_DucT2_Sub12 <- top_specific_markers_DucT2[top_specific_markers_DucT2$cell_group =="12",]

marker_test_res_DucT2_Sub1 <- marker_test_res_DucT2[marker_test_res_DucT2$cell_group =="1",]
marker_test_res_DucT2_Sub2 <- marker_test_res_DucT2[marker_test_res_DucT2$cell_group =="2",]


#
assigned_DucT2_marker_test_res <- top_markers(cds_subset_K100,
                                             group_cells_by="cluster",
                                            # reference_cells=1000,
                                             cores=8)

# Require that markers have at least JS specificty score > 0.5 and
# be significant in the logistic test for identifying their cell type:
garnett_markers_DucT2 <- assigned_DucT2_marker_test_res %>%
  filter(marker_test_q_value < 0.01 & specificity >= 0.1) %>%
  group_by(cell_group) %>%
  top_n(100, marker_score)
# # Exclude genes that are good markers for more than one cell type:
# garnett_markers_DucT2 <- garnett_markers_DucT2 %>% 
#   group_by(gene_short_name) %>%
#   filter(n() == 2)
generate_garnett_marker_file(garnett_markers_DucT2,max_genes_per_group = 100, file="./DucT2_marker_file2.txt")



########################  DucT2_TOP2ACenter ##########################
cds_sub_DucT2_TOP2ACenter <- choose_cells(cds_subset_K100)
plot_cells(cds_sub_DucT2_TOP2ACenter, label_cell_groups=FALSE)
plot_cells(cds_sub_DucT2_TOP2ACenter, label_cell_groups=FALSE, show_trajectory_graph = FALSE, cell_size = 2)
plot_cells(cds_sub_DucT2_TOP2ACenter, label_cell_groups=FALSE, show_trajectory_graph = FALSE, cell_size = 2, color_cells_by="cell_cycle")+ scale_color_manual(values = colorsT)
#plot_cells(cds_sub_DucT2_TOP2ACenter, genes=c(Main),cell_size=1,label_cell_groups = FALSE, show_trajectory_graph = FALSE)

CCdata_sub_DucT2_TOP2ACenter <- cds_sub_DucT2_TOP2ACenter@assays@data@listData$counts
colnames(CCdata_sub_DucT2_TOP2ACenter) = cds_sub_DucT2_TOP2ACenter@assays@data@listData[["counts"]]@Dimnames[[2]]
rownames(CCdata_sub_DucT2_TOP2ACenter) = cds_sub_DucT2_TOP2ACenter@assays@data@listData[["counts"]]@Dimnames[[1]]


DataCellcycle_sub_DucT2_TOP2ACenter <- CCdata_sub_DucT2_TOP2ACenter
DataCellcycle_sub_DucT2_TOP2ACenter <- as(as.matrix(DataCellcycle_sub_DucT2_TOP2ACenter), "dgCMatrix")

marrow_sub_DucT2_TOP2ACenter <- CreateSeuratObject(counts = DataCellcycle_sub_DucT2_TOP2ACenter)
marrow_sub_DucT2_TOP2ACenter <- NormalizeData(marrow_sub_DucT2_TOP2ACenter)
#marrow_sub_DucT2_TOP2ACenter <- FindVariableFeatures(marrow_sub_DucT2_TOP2ACenter, selection.method = "vst")
marrow_sub_DucT2_TOP2ACenter <- FindVariableFeatures(marrow_sub_DucT2_TOP2ACenter, selection.method = "vst")
marrow_sub_DucT2_TOP2ACenter <- ScaleData(marrow_sub_DucT2_TOP2ACenter, features = rownames(marrow_sub_DucT2_TOP2ACenter))
#錯誤: 無法配置大小為 3.9 Gb 的向量
#https://d.cosx.org/d/413001-413001
#memory.limit(15000)


marrow_sub_DucT2_TOP2ACenter <- RunPCA(marrow_sub_DucT2_TOP2ACenter, features = VariableFeatures(marrow_sub_DucT2_TOP2ACenter), ndims.print = 6:10, nfeatures.print = 10)
marrow_sub_DucT2_TOP2ACenter@assays[["RNA"]]@counts@Dimnames[[1]] <- cds_sub_DucT2_TOP2ACenter@assays@data@listData[["counts"]]@Dimnames[[1]]
marrow_sub_DucT2_TOP2ACenter@assays[["RNA"]]@data@Dimnames[[1]] <- cds_sub_DucT2_TOP2ACenter@assays@data@listData[["counts"]]@Dimnames[[1]]
#png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DimHeatmap.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter, dims = c(3, 4))
#dev.off() # 關閉輸出圖檔

#png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DimHeatmap2.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2_TOP2ACenter, dims = c(1, 2))
#dev.off() # 關閉輸出圖檔

###########TTT Ori cell cycle
# marrow_sub_DucT2_TOP2ACenter@active.ident <- cds_sub_DucT2_TOP2ACenter@colData@listData$cell_cycle 
# marrow_sub_DucT2_TOP2ACenter@meta.data[["Phase"]] <- cds_sub_DucT2_TOP2ACenter@colData@listData[["cell_cycle"]]


# Also read in a list of cell cycle markers, from Tirosh et al, 2015
# cc.genes <- readLines(con = "~/Downloads/seurat_resources/regev_lab_cell_cycle_genes.txt")
cc.genes <- read.csv(paste0(PathName,"/Cell cycle/regev_lab_cell_cycle_genesCh.csv"))
# We can segregate this list into markers of G2/M phase and markers of S
# phase
s.genes <- cc.genes[1:43,]
s.genes2 <- capitalize(tolower(s.genes))
# https://www.360kuai.com/pc/9ae634cb511f34aee?cota=4&kuai_so=1&tj_url=so_rec&sign=360_7bc3b157
library(tidyverse)
G_listCCs.genes <- AnnotationDbi::select(org.Mm.eg.db, keys=s.genes2, columns='ENSEMBL', keytype='SYMBOL')
# G_listCCs.genes <- getBM(filters= "hgnc_symbol", attributes= c("ensembl_gene_id","hgnc_symbol"),values=s.genes,mart= martCCG,uniqueRows= TRUE,
#                          quote = "\"")
G_listCCs.genes2 = G_listCCs.genes[!duplicated(G_listCCs.genes[2]),]
G_listCCs.genes3 <- na.omit(G_listCCs.genes2[2])
G_listCCs.genes4 <- G_listCCs.genes3[,1]
# (Ori)s.genes <- cc.genes[,1:43]
g2m.genes <- cc.genes[44:97,]
g2m.genes2 <- capitalize(tolower(g2m.genes))
G_listCCg2m.genes <- AnnotationDbi::select(org.Mm.eg.db, keys=g2m.genes2, columns='ENSEMBL', keytype='SYMBOL')
G_listCCg2m.genes2 = G_listCCg2m.genes[!duplicated(G_listCCg2m.genes[2]),]
G_listCCg2m.genes3 <- na.omit(G_listCCg2m.genes2[2])
G_listCCg2m.genes4 <-G_listCCg2m.genes3[,1]
#marrow <- CellCycleScoring(marrow, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
#marrow <- CellCycleScoring(marrow, s.features = G_listCCs.genes2, g2m.features = G_listCCg2m.genes2, set.ident = TRUE)



marrow_sub_DucT2_TOP2ACenter <- CellCycleScoring(marrow_sub_DucT2_TOP2ACenter, s.features = G_listCCs.genes4, g2m.features = G_listCCg2m.genes4, set.ident = TRUE)
marrow_sub_DucT2_TOP2ACenter <- CellCycleScoring(marrow_sub_DucT2_TOP2ACenter, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(marrow_sub_DucT2_TOP2ACenter[[]])


##############
colorsT <- c("#FF9912B3", "#32CD3299", "#4169E1B3")
RidgePlot(marrow_sub_DucT2_TOP2ACenter,cols = colorsT, features = c("TOP2A"), ncol = 1)


png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_RidgePlot_Sub_DucT2_TOP2ACenter_V2.png")) # 設定輸出圖檔

RidgePlot(marrow_sub_DucT2_TOP2ACenter,cols = colorsT, features = c(Main), ncol = 1)
dev.off() # 關閉輸出圖檔

RidgePlot(marrow_sub_DucT2_TOP2ACenter,cols = colorsT, features = c("TOP2A"), ncol = 1)
# https://datavizpyr.com/ridgeline-plot-with-ggridges-in-r/
# https://www.r-graph-gallery.com/294-basic-ridgeline-plot.html
RidgePlot(marrow_sub_DucT2_TOP2ACenter,cols = colorsT, features = c(Main), ncol = 2,log=TRUE)
RidgePlot(marrow_sub_DucT2_TOP2ACenter,cols = colorsT, features = c(Main), ncol = 2,y.max = 100)

#將Seurat跑出的Cell cycle結果寫入Monocle3的cds檔
cds_sub_DucT2_TOP2ACenter@colData@listData$cell_cycle <- marrow_sub_DucT2_TOP2ACenter@active.ident


##
Maingroup_ciliated_genes <- c("TOP2A")
cds_sub_DucT2_TOP2ACenter <- cds_sub_DucT2_TOP2ACenter[rowData(cds_sub_DucT2_TOP2ACenter)$gene_short_name %in% Maingroup_ciliated_genes,]

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

plot_cells(cds_sub_DucT2_TOP2ACenter_T1, color_cells_by="cell_cycle",cell_size=2, label_cell_groups=FALSE) + scale_color_manual(values = colorsT)


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



########################  Check DucT2 cell cycle ##########################
#cds_sub_DucT2 <- choose_cells(cds)

CCdata_sub_DucT2 <- cds_sub_DucT2@assays@data@listData$counts
colnames(CCdata_sub_DucT2) = cds_sub_DucT2@assays@data@listData[["counts"]]@Dimnames[[2]]
rownames(CCdata_sub_DucT2) = cds_sub_DucT2@assays@data@listData[["counts"]]@Dimnames[[1]]


DataCellcycle_sub_DucT2 <- CCdata_sub_DucT2
DataCellcycle_sub_DucT2 <- as(as.matrix(DataCellcycle_sub_DucT2), "dgCMatrix")

marrow_sub_DucT2 <- CreateSeuratObject(counts = DataCellcycle_sub_DucT2)
marrow_sub_DucT2 <- NormalizeData(marrow_sub_DucT2)
marrow_sub_DucT2 <- FindVariableFeatures(marrow_sub_DucT2, selection.method = "vst")
marrow_sub_DucT2 <- ScaleData(marrow_sub_DucT2, features = rownames(marrow_sub_DucT2))
#錯誤: 無法配置大小為 3.9 Gb 的向量
#https://d.cosx.org/d/413001-413001
#memory.limit(15000)


marrow_sub_DucT2 <- RunPCA(marrow_sub_DucT2, features = VariableFeatures(marrow_sub_DucT2), ndims.print = 6:10, nfeatures.print = 10)
marrow_sub_DucT2@assays[["RNA"]]@counts@Dimnames[[1]] <- cds_sub_DucT2@assays@data@listData[["counts"]]@Dimnames[[1]]
marrow_sub_DucT2@assays[["RNA"]]@data@Dimnames[[1]] <- cds_sub_DucT2@assays@data@listData[["counts"]]@Dimnames[[1]]
#png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DimHeatmap.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2, dims = c(3, 4))
#dev.off() # 關閉輸出圖檔

#png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DimHeatmap2.png")) # 設定輸出圖檔
DimHeatmap(marrow_sub_DucT2, dims = c(1, 2))
#dev.off() # 關閉輸出圖檔

marrow_sub_DucT2 <- CellCycleScoring(marrow_sub_DucT2, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(marrow_sub_DucT2[[]])


##############
colorsT <- c("#FF9912B3", "#32CD3299", "#4169E1B3")

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_RidgePlot_Sub_DucT2.png")) # 設定輸出圖檔

RidgePlot(marrow_sub_DucT2,cols = colorsT, features = c(Main), ncol = 1)
dev.off() # 關閉輸出圖檔

# https://datavizpyr.com/ridgeline-plot-with-ggridges-in-r/
# https://www.r-graph-gallery.com/294-basic-ridgeline-plot.html
RidgePlot(marrow_sub_DucT2,cols = colorsT, features = c(Main), ncol = 2,log=TRUE)
RidgePlot(marrow_sub_DucT2,cols = colorsT, features = c(Main), ncol = 2,y.max = 100)

#將Seurat跑出的Cell cycle結果寫入Monocle3的cds檔
cds_sub_DucT2@colData@listData$cell_cycle_DuctT2 <- marrow_sub_DucT2@active.ident


##
Maingroup_ciliated_genes <- c("TOP2A")
cds_sub_DucT2 <- cds_sub_DucT2[rowData(cds_sub_DucT2)$gene_short_name %in% Maingroup_ciliated_genes,]

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_Violin_Main_Sub_DucT2.png")) # 設定輸出圖檔
plot_genes_violin(cds_sub_DucT2, group_cells_by="cell_cycle", ncol=2, log_scale = FALSE) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
dev.off() # 關閉輸出圖檔

##
png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_CellCycle_Sub_DucT2_Ori.png")) # 設定輸出圖檔
plot_cells(cds_sub_DucT2, color_cells_by="cell_cycle", label_cell_groups=FALSE, show_trajectory_graph = FALSE) + scale_color_manual(values = colorsT)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_CellCycle_Sub_DucT2_Re.png")) # 設定輸出圖檔
plot_cells(cds_sub_DucT2, color_cells_by="cell_cycle_DuctT2", label_cell_groups=FALSE, show_trajectory_graph = FALSE) + scale_color_manual(values = colorsT)
dev.off() # 關閉輸出圖檔

plot_cells(cds_sub_DucT2, color_cells_by="cell_cycle",cell_size=2, label_cell_groups=FALSE) + scale_color_manual(values = colorsT)



plot_cells(cds_sub_DucT2, genes=c("TOP2A"),cell_size=2, label_cell_groups=FALSE, show_trajectory_graph = FALSE)

AFD_lineage_cds_sub_DucT2 <- cds_sub_DucT2[rowData(cds_sub_DucT2)$gene_short_name %in% Main]
plot_genes_in_pseudotime(AFD_lineage_cds_sub_DucT2,
                         color_cells_by="cell_cycle",cell_size=2,
                         min_expr=0.5)+ scale_color_manual(values = colorsT)



########################  Heterogeneity center and Ori Ductal2 ##########################
cds_sub_HeteroCent_OriDucT2 <- choose_cells(cds)
#cds_subset <- reduce_dimension(cds_subset)
plot_cells(cds_sub_HeteroCent_OriDucT2, label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 2)

cds_sub_HeteroCent_K100 <- cluster_cells(cds_sub_HeteroCent_OriDucT2,k = 100, resolution=1e-5)
png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_","_sub_HeteroCent_OriDucT_cluster_clu_K100_1",".png")) # 設定輸出圖檔
plot_cells(cds_sub_HeteroCent_K100, label_cell_groups=FALSE, color_cells_by = "cluster", show_trajectory_graph = FALSE,cell_size = 2)
dev.off() # 關閉輸出圖檔

##  Find marker genes expressed by each cluster
marker_test_HeteroCent_OriDucT2 <- top_markers(cds_sub_HeteroCent_K100, group_cells_by="cluster", 
                                     reference_cells=1000, cores=8)

top_specific_markers_HeteroCent_OriDucT2 <- marker_test_HeteroCent_OriDucT2 %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(25, pseudo_R2)

top_specific_marker_ids_HeteroCent_OriDucT2  <- unique(top_specific_markers_HeteroCent_OriDucT2  %>% pull(gene_id))


plot_genes_by_group(cds_sub_HeteroCent_K100,
                    top_specific_marker_ids_HeteroCent_OriDucT2,
                    group_cells_by="cluster",
                    ordering_type="maximal_on_diag",
                    max.size=3)

top_specific_markers_HeteroCent_OriDucT2 <- data.frame(top_specific_markers_HeteroCent_OriDucT2)

top_marker_HeteroCent_OriDucT2_Sub1 <- top_specific_markers_HeteroCent_OriDucT2[top_specific_markers_HeteroCent_OriDucT2$cell_group =="1",]
top_marker_HeteroCent_OriDucT2_Sub2 <- top_specific_markers_HeteroCent_OriDucT2[top_specific_markers_HeteroCent_OriDucT2$cell_group =="2",]










########################  cds_subset ##########################


# ############  Annotate your cells according to type (Custom Marker)  ############
# Cell type gene expression markers https://panglaodb.se/markers.html
# https://www.ncbi.nlm.nih.gov/mesh?Db=mesh&Cmd=DetailsSearch&Term=%22Genetic+Markers%22%5BMeSH+Terms%5D
#
library(monocle3)
library(garnett)
## Install the monocle3 branch of garnett
# BiocManager::install(c("org.Mm.eg.db", "org.Hs.eg.db"))
# devtools::install_github("cole-trapnell-lab/garnett", ref="monocle3")
# Install problem
# https://stackoverflow.com/questions/42807247/installing-package-cannot-open-file-permission-denied

colData(cds_subset_K100)$garnett_cluster <- clusters(cds_subset_K100)

# Installing Cicero
# https://cole-trapnell-lab.github.io/cicero-release/docs_m3/#installing-cicero
library(cicero)
Human_classifier <- train_cell_classifier(cds = cds_subset_K100,
                                          marker_file = Marker_file,   # Import the marker_file
                                          db=org.Hs.eg.db::org.Hs.eg.db,
                                          cds_gene_id_type = "SYMBOL",
                                          #num_unknown = 2215,
                                          #max_training_samples = 10000,
                                          marker_file_gene_id_type = "SYMBOL",
                                          cores=8)


cds_subset <- classify_cells(cds_subset_K100, Human_classifier,
                      db = org.Hs.eg.db::org.Hs.eg.db,
                      cluster_extend = TRUE,
                      cds_gene_id_type = "SYMBOL")



png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_",Marker_file_Name,".png")) # 設定輸出圖檔
plot_cells(cds_subset,
           group_cells_by="cluster",
           cell_size=1.5,
           color_cells_by="cluster_ext_type", show_trajectory_graph = FALSE)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_",Marker_file_Name,"_2.png")) # 設定輸出圖檔
plot_cells(cds_subset,
           group_cells_by="cluster",
           cell_size=1.5,
           color_cells_by="cluster_ext_type",
           label_cell_groups=FALSE, show_trajectory_graph = FALSE)
dev.off() # 關閉輸出圖檔
############  Annotate your cells according to type (Custom Marker)  ############


########################  Constructing single-cell trajectories ##########################
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

cds <- order_cells(cds)

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)


AFD_genes <- c("TOP2A","TOP2B", "MUC1", "BRIP1","KRAS")

AFD_lineage_cds <- cds[rowData(cds)$gene_short_name %in% AFD_genes]

plot_genes_in_pseudotime(AFD_lineage_cds,
                         color_cells_by="cell_cycle",cell_size=2,
                         min_expr=0.5)

####################
cds_subTra <- choose_graph_segments(cds)
cds_subTra <- preprocess_cds(cds_subTra, num_dim = 100)
cds_subTra <- reduce_dimension(cds_subTra)
cds_subTra <- cluster_cells(cds_subTra)
cds_subTra <- learn_graph(cds_subTra)

cds_subTra <- order_cells(cds_subTra)
plot_cells(cds_subTra,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

# my_genes <- row.names(subset(fData(cds_subTra),
#                               gene_short_name %in% c("TOP2A", "MUC1", "BRIP1","KRAS")))
#
# mos_subset <- cds_subTra[my_genes,]
#
# heatmap <- plot_pseudotime_heatmap(mos_subset,
#                                    num_clusters = 1,
#                                    show_rownames = FALSE, return_heatmap = T)

CCdata2 <- cds_subTra@assays@data@listData$counts
colnames(CCdata2) = cds_subTra@assays@data@listData[["counts"]]@Dimnames[[2]]
rownames(CCdata2) = cds_subTra@assays@data@listData[["counts"]]@Dimnames[[1]]


DataCellcycle2 <- CCdata2
#rownames(DataCellcycle2) <- make.names(DataCellcycle2[,1], unique = TRUE)
#DataCellcycle <- DataCellcycle[1:length(data[,1]), 2:length(data[1,])]
DataCellcycle2 <- as(as.matrix(DataCellcycle2), "dgCMatrix")

marrow2 <- CreateSeuratObject(counts = DataCellcycle2)
marrow2 <- NormalizeData(marrow2)
marrow2 <- FindVariableFeatures(marrow2, selection.method = "vst")
marrow2 <- ScaleData(marrow2, features = rownames(marrow2))
#錯誤: 無法配置大小為 3.9 Gb 的向量
#https://d.cosx.org/d/413001-413001
#memory.limit(15000)


marrow2 <- RunPCA(marrow2, features = VariableFeatures(marrow2), ndims.print = 6:10, nfeatures.print = 10)
marrow2@assays[["RNA"]]@counts@Dimnames[[1]] <- cds@assays@data@listData[["counts"]]@Dimnames[[1]]
marrow2@assays[["RNA"]]@data@Dimnames[[1]] <- cds@assays@data@listData[["counts"]]@Dimnames[[1]]
#png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DimHeatmap.png")) # 設定輸出圖檔
DimHeatmap(marrow2, dims = c(3, 4))
#dev.off() # 關閉輸出圖檔

#png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DimHeatmap2.png")) # 設定輸出圖檔
DimHeatmap(marrow2, dims = c(1, 2))
#dev.off() # 關閉輸出圖檔

# Also read in a list of cell cycle markers, from Tirosh et al, 2015
# cc.genes <- readLines(con = "~/Downloads/seurat_resources/regev_lab_cell_cycle_genes.txt")
cc.genes <- read.csv(paste0(PathName,"/Cell cycle/regev_lab_cell_cycle_genesCh.csv"))
# We can segregate this list into markers of G2/M phase and markers of S
# phase
s.genes <- cc.genes[1:43,]
s.genes2 <- capitalize(tolower(s.genes))
# https://www.360kuai.com/pc/9ae634cb511f34aee?cota=4&kuai_so=1&tj_url=so_rec&sign=360_7bc3b157
library(tidyverse)
G_listCCs.genes <- AnnotationDbi::select(org.Mm.eg.db, keys=s.genes2, columns='ENSEMBL', keytype='SYMBOL')
# G_listCCs.genes <- getBM(filters= "hgnc_symbol", attributes= c("ensembl_gene_id","hgnc_symbol"),values=s.genes,mart= martCCG,uniqueRows= TRUE,
#                          quote = "\"")
G_listCCs.genes2 = G_listCCs.genes[!duplicated(G_listCCs.genes[2]),]
G_listCCs.genes3 <- na.omit(G_listCCs.genes2[2])
G_listCCs.genes4 <- G_listCCs.genes3[,1]
# (Ori)s.genes <- cc.genes[,1:43]
g2m.genes <- cc.genes[44:97,]
g2m.genes2 <- capitalize(tolower(g2m.genes))
G_listCCg2m.genes <- AnnotationDbi::select(org.Mm.eg.db, keys=g2m.genes2, columns='ENSEMBL', keytype='SYMBOL')
G_listCCg2m.genes2 = G_listCCg2m.genes[!duplicated(G_listCCg2m.genes[2]),]
G_listCCg2m.genes3 <- na.omit(G_listCCg2m.genes2[2])
G_listCCg2m.genes4 <-G_listCCg2m.genes3[,1]
#marrow <- CellCycleScoring(marrow, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
#marrow <- CellCycleScoring(marrow, s.features = G_listCCs.genes2, g2m.features = G_listCCg2m.genes2, set.ident = TRUE)



marrow2 <- CellCycleScoring(marrow2, s.features = G_listCCs.genes4, g2m.features = G_listCCg2m.genes4, set.ident = TRUE)
marrow2 <- CellCycleScoring(marrow2, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(marrow2[[]])


##############
colorsT <- c("#FF9912B3", "#32CD3299", "#4169E1B3")

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_RidgePlot_Main_Tr1.png")) # 設定輸出圖檔

RidgePlot(marrow2,cols = colorsT, features = c(Main), ncol = 1)
dev.off() # 關閉輸出圖檔

# https://datavizpyr.com/ridgeline-plot-with-ggridges-in-r/
# https://www.r-graph-gallery.com/294-basic-ridgeline-plot.html
RidgePlot(marrow2,cols = colorsT, features = c(Main), ncol = 2,log=TRUE)
RidgePlot(marro2w,cols = colorsT, features = c(Main), ncol = 2,y.max = 100)

#將Seurat跑出的Cell cycle結果寫入Monocle3的cds檔
cds_subTra@colData@listData$cell_cycle <- marrow2@active.ident


##
Main <- c("TOP2A","TOP2B","PHGR1","NDRG1","FABP1","MKI67","TUBA1B","BIRC5")
Maingroup_ciliated_genes <- c(Main)
cds_subTra2 <- cds_subTra[rowData(cds_subTra)$gene_short_name %in% Maingroup_ciliated_genes,]

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_Violin_Main_Tr1.png")) # 設定輸出圖檔
plot_genes_violin(cds_subTra2, group_cells_by="cell_cycle", ncol=2, log_scale = FALSE) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
dev.off() # 關閉輸出圖檔

##
png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_CellCycle_Tr1.png")) # 設定輸出圖檔
plot_cells(cds_subTra2, color_cells_by="cell_cycle",cell_size=1, label_cell_groups=FALSE) + scale_color_manual(values = colorsT)
dev.off() # 關閉輸出圖檔

plot_cells(cds_subTra2, color_cells_by="cell_cycle",cell_size=2, label_cell_groups=FALSE) + scale_color_manual(values = colorsT)



plot_cells(cds_subTra2, genes=c("TOP2A"),cell_size=2, label_cell_groups=FALSE, show_trajectory_graph = FALSE)

AFD_lineage_cds <- cds_subTra2[rowData(cds_subTra2)$gene_short_name %in% Main]
plot_genes_in_pseudotime(AFD_lineage_cds,
                         color_cells_by="cell_cycle",cell_size=2,
                         min_expr=0.5)+ scale_color_manual(values = colorsT)

############(ERROR) Finding modules of co-regulated genes ############
ciliated_cds_pr_test_res <- graph_test(cds_subTra2, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))
gene_module_df <- find_gene_modules(cds_subTra2[pr_deg_ids,], resolution=1e-2)

cell_group_df <- tibble::tibble(cell=row.names(colData(cds_subTra2)),
                                cell_group=partitions(cds)[colnames(cds_subTra2)])
agg_mat <- aggregate_gene_expression(cds_subTra2, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c("Partition ", colnames(agg_mat))

pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                   scale="column", clustering_method="ward.D2",
                   fontsize=6)



############ Working with 3D trajectories ############
cds_3d <- reduce_dimension(cds, max_components = 3)
cds_3d <- cluster_cells(cds_3d)
cds_3d <- learn_graph(cds_3d)
cds_3d <- order_cells(cds_3d, root_pr_nodes=get_earliest_principal_node(cds))

cds_3d_plot_obj <- plot_cells_3d(cds_3d, color_cells_by="partition")



