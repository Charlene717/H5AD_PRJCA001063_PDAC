## (scRNA-seq data analysis for H5AD files)

#############
rm(list = ls()) # Clean variable

memory.limit(150000)

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
# library(cicero)
# detach("package:monocle", unload = TRUE)
# 錯誤: package 'monocle' is required by 'cicero' so will not be detached

############# Import files settings #############
## General setting
PathName = setwd(getwd())
RVersion = "20210523V1"
dir.create(paste0(PathName,"/",RVersion))

## Marker genes file
Marker_file_Name <- c("NAKAMURA_METASTASIS_MODEL_UP")
Marker_file <- paste0(PathName,"/",Marker_file_Name,".txt")
Marker_List <- read.delim(Marker_file,header=F,sep= c("\t"))
Marker_List2 <- as.data.frame(Marker_List[-1:-2,])

Garnett_Marker_file_Name <- c("NAKAMURA_METASTASIS_MODEL_M18483")
Garnett_Marker_file <- paste0(PathName,"/marker_file_",Garnett_Marker_file_Name,".txt")

## Cell cycle genes file
cc.genes_list <- read.csv(paste0(PathName,"/Cell cycle/regev_lab_cell_cycle_genesCh.csv")) # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat. 
# cc: Cell-Cycle

############# Marker genes file (Old Version) #############
# Marker_file_Name <- c("NAKAMURA_METASTASIS_MODEL_M18483")
# Marker_file <- paste0(PathName,"/marker_file_",Marker_file_Name,".txt")
# Marker_List <- read.table(Marker_file,header=F,sep= c(","),stringsAsFactors = FALSE, fill = TRUE)
# library(stringr)
# Marker_List_1 <- Marker_List[2,1]
# Marker_List_2 <- str_replace_all(Marker_List_1,"expressed: ","")
# Marker_List <- str_trim(Marker_List[2,-1], side = c("both"))
# Marker_List <- c(Marker_List_2,Marker_List)
############# Marker genes file (Old Version) #############


############# Parameter setting #############
## Gene list of interest 
Main = c("TOP2A")
Main_Group = c("TOP2A","TOP2B","TP53","CCNE1")
Main_Group2 = c("KRAS","EXO1","NSUN2","MUC1","AMBP","FXYD2")
EMT_Meta = c("ANLN","APLP2","CD63","CDH2","CLIC4","CTSB","CX3CR1","DSG2","EDNRB")
candidates14 = c("BRIP1","KIF23","TOP2A","FOSL1","FAM25A","ANLN","NCAPH","KRT9","MCM4","CKAP2L","CENPE","RACGAP1","DTL","RAD51AP1")

DREAM_complex= c("RBL2","E2F4","E2F5","TFDP1","TFDP2")
Regulators= c("TP53","YBX1","E2F1")

## Color setting
colors_cc <- c("#FF9912B3", "#32CD3299", "#4169E1B3") ## Color for Cell-Cycle

## Format of data
GeneNAFMT <- c("HuGSymbol") # Gene names format of data: HuGSymbol,MouGSymbol,HuENSEMBL,MouENSEMBL

## Cluster cells setting
k_cds_sub_DucT2 <- c(100) # k for ductal cell type2

## Threshold of PCA scores
PCAThreshold_Pos <- 0.03
PCAThreshold_Neg <- -0.03

##################  Grab specific terms ################## 
## grepl Stroma
Stroma_cds <- cds[,grepl("Stromal", colData(cds)$Broad.cell.type, ignore.case=TRUE)]
plot_cells(Stroma_cds, reduction_method="tSNE", color_cells_by="partition")


##################  Function setting ################## 

## Call function
filePath <- ""
#匯入 同一個資料夾中的R檔案
getFilePath <- function(fileName) {
  # path <- setwd("~")  #專案資料夾絕對路徑
  path <- setwd(getwd()) 
  #字串合併無間隔
  # 「<<-」為全域變數給值的指派
  filePath <<- paste0(path ,"/" , fileName)  
  # 載入檔案
  sourceObj <- source(filePath)
  return(sourceObj)
}



#######################################################################################################


############# Import raw data #############
library(SummarizedExperiment)
library(Seurat)
library(SeuratDisk)

## Convert h5ad to h5seurat
Convert(paste0(PathName,"/StdWf1_PRJCA001063_CRC_besca2.raw.h5ad"), "PRJCA001063.h5seurat", assay = "RNA",) # This creates a copy of this .h5ad object reformatted into .h5seurat inside the example_dir directory
seuratObject <- LoadH5Seurat(paste0(PathName,"/PRJCA001063.h5seurat")) # This .d5seurat object can then be read in manually
 
## Convert Seurat Object to Monocle3 Object
library(SeuratWrappers)
cds <- as.cell_data_set(seuratObject) # Convert objects to Monocle3 'cell_data_set' objects

############# Run Monocle3 #############
library(monocle3)
cds <- estimate_size_factors(cds) # issues with cds object in monocle3 #54 # https://github.com/satijalab/seurat-wrappers/issues/54

###### Pre-process the data ######
cds <- preprocess_cds(cds, num_dim = 100)
plot_pc_variance_explained(cds)
# Error in plot_pc_variance_explained function
# https://github.com/cole-trapnell-lab/monocle3/issues/170
# package 'monocle' is required by 'cicero'


###### Reduce dimensionality and visualize the cells ######

##### UMAP #####
cds <- reduce_dimension(cds,preprocess_method = 'PCA')
plot_cells(cds)

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_Ori.png")) # Set output image file
plot_cells(cds)
dev.off() # Close the output image file

##### (UMAP) Plot different phenotype #####
png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_CellType.png")) # 設定輸出圖檔
plot_cells(cds, color_cells_by="Cell_type", label_cell_groups=FALSE, show_trajectory_graph = FALSE)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_Type.png")) # 設定輸出圖檔
plot_cells(cds, color_cells_by="Type", label_cell_groups=FALSE, show_trajectory_graph = FALSE)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_Patient.png")) # 設定輸出圖檔
plot_cells(cds, color_cells_by="Patient", label_cell_groups=FALSE, show_trajectory_graph = FALSE)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_CONDITION.png")) # 設定輸出圖檔
plot_cells(cds, color_cells_by="CONDITION", label_cell_groups=FALSE, show_trajectory_graph = FALSE)
dev.off() # 關閉輸出圖檔

##### (UMAP) Plot genes #####
plot_cells(cds, genes=c("TOP2A","TOP2B","TP53")) #error
cds@rowRanges@elementMetadata@listData$gene_short_name <- cds@assays@data@listData[["counts"]]@Dimnames[[1]]

plot_cells(cds, genes=c("TOP2A","TOP2B","TP53","CCNE1")) #ok
plot_cells(cds, genes=c("CGAS","TOP2A"),cell_size = 1)

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_",Main,".png")) # 設定輸出圖檔
plot_cells(cds, genes=c(Main),cell_size=1,label_cell_groups = FALSE, show_trajectory_graph = FALSE)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_","Main_Group",".png")) # 設定輸出圖檔
plot_cells(cds, genes=c(Main_Group),cell_size=0.5,label_cell_groups = FALSE, show_trajectory_graph = FALSE)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_","Main_Group2",".png")) # 設定輸出圖檔
plot_cells(cds, genes=c(Main_Group2),cell_size=0.5,label_cell_groups = FALSE, show_trajectory_graph = FALSE)
dev.off() # 關閉輸出圖檔


##### Group cells into clusters ######
cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "partition", label_cell_groups=FALSE, show_trajectory_graph = FALSE)
plot_cells(cds, color_cells_by = "cluster", label_cell_groups=FALSE, show_trajectory_graph = FALSE)
plot_cells(cds, label_cell_groups=FALSE, show_trajectory_graph = FALSE)

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_","cluster_par",".png")) # 設定輸出圖檔
plot_cells(cds, color_cells_by = "partition", label_cell_groups=FALSE, show_trajectory_graph = FALSE)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_","cluster_clu",".png")) # 設定輸出圖檔
plot_cells(cds, color_cells_by = "cluster", label_cell_groups=FALSE, show_trajectory_graph = FALSE)
dev.off() # 關閉輸出圖檔


############# Run Seurat #############

############ Cell-Cycle Scoring and Regression - Monocle3 & Seurat Mutual conversion #############

## Load package
library(Seurat)
library(SummarizedExperiment) 
library(AnnotationDbi)
library(org.Mm.eg.db)
library('org.Hs.eg.db')
library(Hmisc)

###### Convert Monocle3 Object to Seurat Object ######
getFilePath("Monocle3_To_Seurat.R")
marrow <- Monocle3_To_Seurat(cds,"cds") #這個function存在於Monocle3_To_Seurat.R裡面

###### Assign Cell-Cycle Scores ######
getFilePath("Cell-Cycle Scoring and Regression.R")
marrow <- CCScorReg(GeneNAFMT,marrow) #這個function存在於Cell-Cycle Scoring and Regression.R裡面

RidgePlot(marrow,cols = colors_cc, features = c(Main), ncol = 1)
RidgePlot(marrow,cols = colors_cc, features = c(Main_Group), ncol = 2,log=TRUE) 
RidgePlot(marrow,cols = colors_cc, features = c(Main_Group), ncol = 2,y.max = 100) 

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","cds_CC_RidgePlot_Main.png")) # 設定輸出圖檔
RidgePlot(marrow,cols = colors_cc, features = c(Main), ncol = 1) 
dev.off() # 關閉輸出圖檔

###### Insert the cell cycle results from Seurat into the  Monocle3 cds object ######
cds@colData@listData$cell_cycle <- marrow@active.ident
# cds@colData@listData$cell_cycle <- marrow@meta.data[["Phase"]]

plot_cells(cds, color_cells_by="cell_cycle", label_cell_groups=FALSE) + scale_color_manual(values = colors_cc)

## Plot the violin diagram
Maingroup_ciliated_genes <- c(Main_Group)
cds_marrow_cc <- cds[rowData(cds)$gene_short_name %in% Maingroup_ciliated_genes,]

plot_genes_violin(cds_marrow_cc, group_cells_by="cell_cycle", ncol=2, log_scale = FALSE)+ scale_fill_manual(values = colors_cc)
plot_genes_violin(cds_marrow_cc, group_cells_by="cell_cycle", ncol=2, log_scale = T)+ scale_fill_manual(values = colors_cc)
plot_genes_violin(cds_marrow_cc, group_cells_by="cell_cycle", ncol=2, log_scale = T)+ scale_fill_manual(values = colors_cc)+
  geom_boxplot(width=0.1, fill="white")

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CC_Violin_Main.png")) # 設定輸出圖檔
plot_genes_violin(cds_marrow_cc, group_cells_by="cell_cycle", ncol=2, log_scale = FALSE)+ 
  scale_fill_manual(values = colors_cc)+
  theme(axis.text.x=element_text(angle=45, hjust=1))
dev.off() # 關閉輸出圖檔

pdf(paste0(PathName,"/",RVersion,"/",RVersion,"_","CC_Violin_Main.pdf")) # 設定輸出圖檔
plot_genes_violin(cds_marrow_cc, group_cells_by="cell_cycle", ncol=2, log_scale = FALSE)+ 
  scale_fill_manual(values = colors_cc)+
  theme(axis.text.x=element_text(angle=45, hjust=1))
dev.off() # 關閉輸出圖檔



###### Cell discrimination ######
getFilePath("Monocle3_AddModuleScore.R")

PDAC_Marker_file_Name <- c("GRUETZMANN_PANCREATIC_CANCER_UP")
PDAC_Marker_Name <- c("PDAC_Marker")

cds <- Monocle3_AddModuleScore(PDAC_Marker_file_Name,PDAC_Marker_Name,marrow,cds)
plot_cells(cds, color_cells_by= Marker_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE)

cds_sub_DucT2 <- Monocle3_AddModuleScore(PDAC_Marker_file_Name,PDAC_Marker_Name,marrow_sub_DucT2,cds_sub_DucT2)
plot_cells(cds_sub_DucT2, color_cells_by= Marker_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE)


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


########### Constructing single-cell trajectories ###########
cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "cluster",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

MainGroup_lineage_cds <- cds[rowData(cds)$gene_short_name %in% Main_Group]

plot_genes_in_pseudotime(MainGroup_lineage_cds,
                         color_cells_by="cell_cycle",cell_size=2,
                         min_expr=0.5)+ scale_color_manual(values = colors_cc)


######################################  cds_subset ########################################
########################  DucT2 ##########################
cds_sub_DucT2 <- choose_cells(cds)
#cds_subset <- reduce_dimension(cds_subset)
plot_cells(cds_sub_DucT2, label_cell_groups=FALSE, show_trajectory_graph = FALSE)

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_","_sub_DucT2_cluster_clu_ori",".png")) # 設定輸出圖檔
plot_cells(cds_sub_DucT2, color_cells_by = "cluster", label_cell_groups=FALSE, show_trajectory_graph = FALSE)
dev.off() # 關閉輸出圖檔

cds_subset_NewK <- cluster_cells(cds_sub_DucT2,k = k_cds_sub_DucT2, resolution=1e-5)
plot_cells(cds_subset_NewK, color_cells_by = "cluster", label_cell_groups=FALSE, show_trajectory_graph = FALSE)
 png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_","_sub_DucT2_cluster_cluster_clu_NewK",".png")) # 設定輸出圖檔
 plot_cells(cds_subset_K20, color_cells_by = "cluster", label_cell_groups=FALSE, show_trajectory_graph = FALSE)
 dev.off() # 關閉輸出圖檔

 png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_","_sub_DucT2_cluster_clu_NewK_",Main,".png")) # 設定輸出圖檔
 plot_cells(cds_subset_NewK, genes=c(Main),cell_size=1,label_cell_groups = FALSE, show_trajectory_graph = FALSE)
 dev.off() # 關閉輸出圖檔

 
######  Find marker genes expressed by each cluster ######
marker_test_res_DucT2 <- top_markers(cds_subset_NewK, group_cells_by="cluster")

top_specific_markers_DucT2 <- marker_test_res_DucT2 %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(10, pseudo_R2)

top_specific_marker_ids_DucT2  <- unique(top_specific_markers_DucT2  %>% pull(gene_id))


plot_genes_by_group(cds_subset_NewK,
                    top_specific_marker_ids_DucT2,
                    group_cells_by="cluster",
                    ordering_type="maximal_on_diag",
                    max.size=3)
plot_genes_by_group(cds_subset_NewK,
                    top_specific_marker_ids_DucT2,
                    group_cells_by="cluster",
                    ordering_type="cluster_row_col",
                    max.size=3)


top_specific_markers_DucT2 <- data.frame(top_specific_markers_DucT2)

top_specific_markers_DucT2_Sub1 <- top_specific_markers_DucT2[top_specific_markers_DucT2$cell_group =="1",]
top_specific_markers_DucT2_Sub2 <- top_specific_markers_DucT2[top_specific_markers_DucT2$cell_group =="2",]
#...
marker_test_res_DucT2_Sub1 <- marker_test_res_DucT2[marker_test_res_DucT2$cell_group =="1",]
marker_test_res_DucT2_Sub2 <- marker_test_res_DucT2[marker_test_res_DucT2$cell_group =="2",]
#...

## Export a marker genes information file
write.table(marker_test_res_DucT2, file=paste0(PathName,"/",RVersion,"/",RVersion,"_",
                                                                     "marker_test_res_DucT2.txt"),  sep="\t", row.names=FALSE)

## Generate a Garnett file
assigned_DucT2_marker_test_res <- top_markers(cds_subset_NewK, group_cells_by="cluster")

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

generate_garnett_marker_file(garnett_markers_DucT2,max_genes_per_group = 100, file=paste0(PathName,"/",RVersion,"/",RVersion,"_","DucT2_marker_file2.txt"))
# generate_garnett_marker_file(garnett_markers_DucT2,max_genes_per_group = 100, file="./DucT2_marker_file2.txt")



########################  DucT2_TOP2ACenter ##########################
cds_sub_DucT2_TOP2ACenter <- choose_cells(cds_subset_NewK)
plot_cells(cds_sub_DucT2_TOP2ACenter, label_cell_groups=FALSE)
plot_cells(cds_sub_DucT2_TOP2ACenter, label_cell_groups=FALSE, show_trajectory_graph = FALSE, cell_size = 2)
plot_cells(cds_sub_DucT2_TOP2ACenter, label_cell_groups=FALSE, show_trajectory_graph = FALSE, cell_size = 2, color_cells_by="cell_cycle")+ scale_color_manual(values = colors_cc)
plot_cells(cds_sub_DucT2_TOP2ACenter, genes=c(Main), label_cell_groups=FALSE, show_trajectory_graph = FALSE, cell_size = 2)


###### Convert Monocle3 Object to Seurat Object ######
# getFilePath("Monocle3_To_Seurat.R")
marrow_sub_DucT2_TOP2ACenter <- Monocle3_To_Seurat(cds_sub_DucT2_TOP2ACenter,"sub_DT2TOP2ACTR") #sub_DT2TOP2ACTR:sub_DucT2_TOP2ACenter

###### Assign Cell-Cycle Scores ######
# getFilePath("Cell-Cycle Scoring and Regression.R")
marrow_sub_DucT2_TOP2ACenter <- CCScorReg(GeneNAFMT,marrow_sub_DucT2_TOP2ACenter) #這個function存在於Cell-Cycle Scoring and Regression.R裡面
# view cell cycle scores and phase assignments
head(marrow_sub_DucT2_TOP2ACenter[[]])

## Plot the RidgePlot
RidgePlot(marrow_sub_DucT2_TOP2ACenter,cols = colors_cc, features = c(Main), ncol = 1)
RidgePlot(marrow_sub_DucT2_TOP2ACenter,cols = colors_cc, features = c(Main_Group), ncol = 2,log=TRUE) 
RidgePlot(marrow_sub_DucT2_TOP2ACenter,cols = colors_cc, features = c(Main_Group), ncol = 2,y.max = 100) 

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_RidgePlot_sub_DT2TOP2ACTR_V2.png")) # 設定輸出圖檔
RidgePlot(marrow_sub_DucT2_TOP2ACenter,cols = colorsT, features = c(Main), ncol = 1)
dev.off() # 關閉輸出圖檔


###### Insert the cell cycle results from Seurat into the  Monocle3 cds object ######
cds_sub_DucT2_TOP2ACenter@colData@listData$cell_cycle <- marrow_sub_DucT2_TOP2ACenter@active.ident

plot_cells(cds_sub_DucT2_TOP2ACenter, color_cells_by="cell_cycle",cell_size=2, label_cell_groups=FALSE, show_trajectory_graph = FALSE) + scale_color_manual(values = colors_cc)
plot_cells(cds_sub_DucT2_TOP2ACenter, genes=c("TOP2A"),cell_size=2, label_cell_groups=FALSE, show_trajectory_graph = FALSE)

## Plot the Violin Plot 
cds_sub_DT2TOP2ACTR_Maingroup <- cds_sub_DucT2_TOP2ACenter[rowData(cds_sub_DucT2_TOP2ACenter)$gene_short_name %in% Main_Group,]
plot_genes_violin(cds_sub_DT2TOP2ACTR_Maingroup, group_cells_by="cell_cycle", ncol=2, log_scale = FALSE) +
                  scale_fill_manual(values = colors_cc) + 
                  theme(axis.text.x=element_text(angle=45, hjust=1))


png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_Violin_Main_sub_DT2TOP2ACTR_V2.png")) # 設定輸出圖檔
plot_genes_violin(cds_sub_DT2TOP2ACTR_Maingroup, group_cells_by="cell_cycle", ncol=2, log_scale = FALSE) +
  scale_fill_manual(values = colors_cc) + 
  theme(axis.text.x=element_text(angle=45, hjust=1))
dev.off() # 關閉輸出圖檔

##
png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_CellCycle_sub_DT2TOP2ACTR_V2.png")) # 設定輸出圖檔
plot_cells(cds_sub_DucT2_TOP2ACenter, color_cells_by="cell_cycle",cell_size=3, label_cell_groups=FALSE, show_trajectory_graph = FALSE) + scale_color_manual(values = colors_cc)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_",Main,"_sub_DT2TOP2ACTR_V2.png")) # 設定輸出圖檔
plot_cells(cds_sub_DucT2_TOP2ACenter, genes=c(Main),cell_size=3, label_cell_groups=FALSE, show_trajectory_graph = FALSE)
dev.off() # 關閉輸出圖檔



MainGroup_lineage_sub_DT2TOP2ACTR <- cds_sub_DucT2_TOP2ACenter[rowData(cds_sub_DucT2_TOP2ACenter)$gene_short_name %in% Main_Group]
plot_genes_in_pseudotime(MainGroup_lineage_sub_DT2TOP2ACTR,
                         color_cells_by="cell_cycle",cell_size=2,
                         min_expr=0.5)+ scale_color_manual(values = colors_cc)




########################  DucT2_TOP2ACenter trajectories ##########################
for (i in c(1:8)) {
 
  cds_sub_DucT2_TOP2ACenter_Tn <- choose_graph_segments(cds_sub_DucT2 ,clear_cds = FALSE)
  plot_cells(cds_sub_DucT2_TOP2ACenter_Tn, color_cells_by="cell_cycle",cell_size=2, 
             label_cell_groups=FALSE) + scale_color_manual(values = colors_cc)
  
  ###### Convert Monocle3 Object to Seurat Object ######
  # getFilePath("Monocle3_To_Seurat.R")
  marrow_sub_DucT2_TOP2ACenter_Tn <- Monocle3_To_Seurat(cds_sub_DucT2_TOP2ACenter_Tn,paste0("sub_DT2TOP2ACTR_T", i)) #sub_DT2TOP2ACTR:sub_DucT2_TOP2ACenter
  
  ###### Assign Cell-Cycle Scores ######
  # getFilePath("Cell-Cycle Scoring and Regression.R")
  marrow_sub_DucT2_TOP2ACenter_Tn <- CCScorReg(GeneNAFMT,marrow_sub_DucT2_TOP2ACenter_Tn) #這個function存在於Cell-Cycle Scoring and Regression.R裡面
  assign(paste0("marrow_sub_DucT2_TOP2ACenter_T", i),marrow_sub_DucT2_TOP2ACenter_Tn)
  
  ###### Insert the cell cycle results from Seurat into the  Monocle3 cds object ######
  cds_sub_DucT2_TOP2ACenter_Tn@colData@listData$cell_cycle <- marrow_sub_DucT2_TOP2ACenter_Tn@active.ident
  plot_cells(cds_sub_DucT2_TOP2ACenter_Tn, color_cells_by="cell_cycle", label_cell_groups=FALSE) + scale_color_manual(values = colors_cc)
  assign(paste0("cds_sub_DucT2_TOP2ACenter_T", i),cds_sub_DucT2_TOP2ACenter_Tn)
  
  ###### PCA Scores for finding significantly different genes at the endpoints ######
  getFilePath("PCA_Threshold.R")
  PCA_DT2TOP2ACTR_Tn <- assign(paste0("PCA_DT2TOP2ACTR_T", i),marrow_sub_DucT2_TOP2ACenter_Tn@reductions[["pca"]]@feature.loadings)
  assign(paste0("PCA_DT2TOP2ACTR_T", i,"_PC_Sum"),PCA_Threshold_Pos(PCA_DT2TOP2ACTR_Tn, i ,PCAThreshold_Pos))
  assign(paste0("PCA_DT2TOP2ACTR_T", i,"_NC_Sum"),PCA_Threshold_Neg(PCA_DT2TOP2ACTR_Tn, i ,PCAThreshold_Neg))

  rm(cds_sub_DucT2_TOP2ACenter_Tn,marrow_sub_DucT2_TOP2ACenter_Tn)
  }




########################  Heterogeneity center and Ori Ductal2 ##########################
cds_sub_HeteroCent_OriDucT2 <- choose_cells(cds)
#cds_subset <- reduce_dimension(cds_subset)
plot_cells(cds_sub_HeteroCent_OriDucT2, label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 2)

cds_sub_HeteroCent_K100 <- cluster_cells(cds_sub_HeteroCent_OriDucT2,k = 100, resolution=1e-5)
png(paste0(PathName,"/",RVersion,"/",RVersion,"_","UMAP_","_sub_HeteroCent_OriDucT_cluster_clu_K100_1",".png")) # 設定輸出圖檔
plot_cells(cds_sub_HeteroCent_K100, label_cell_groups=FALSE, color_cells_by = "cluster", show_trajectory_graph = FALSE,cell_size = 2)
dev.off() # 關閉輸出圖檔

##  Find marker genes expressed by each cluster
marker_test_HeteroCent_OriDucT2 <- top_markers(cds_sub_HeteroCent_K100, group_cells_by="cluster")

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
                                          marker_file = Garnett_Marker_file,   # Import the marker_file
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
cds_3d <- reduce_dimension(cds, max_components = 3,preprocess_method = 'PCA')
cds_3d <- cluster_cells(cds_3d)
cds_3d <- learn_graph(cds_3d)
cds_3d <- order_cells(cds_3d, root_pr_nodes=get_earliest_principal_node(cds))
# Error in get_earliest_principal_node(cds) : 
#   沒有這個函數 "get_earliest_principal_node"

cds_3d_plot_obj <- plot_cells_3d(cds_3d, color_cells_by="partition")
plot_cells_3d(cds_3d)
plot_cells_3d(cds_3d, color_cells_by="cluster", show_trajectory_graph = FALSE)
plot_cells_3d(cds_3d, color_cells_by="Cell_type", show_trajectory_graph = FALSE)
plot_cells_3d(cds_3d, color_cells_by="Type", show_trajectory_graph = FALSE)
plot_cells_3d(cds_3d, color_cells_by="Patient", show_trajectory_graph = FALSE)
# plot_cells_3d(cds_3d, color_cells_by="CONDITION", show_trajectory_graph = FALSE)
cds_3d@colData@listData$PDAC_Marker <- marrow_PDAC_Marker@meta.data[["PDAC_Marker1"]]
plot_cells_3d(cds_3d, color_cells_by="PDAC_Marker", show_trajectory_graph = FALSE)
plot_cells_3d(cds_3d, genes = Main, show_trajectory_graph = FALSE)
plot_cells_3d(cds_3d, color_cells_by="cell_cycle", show_trajectory_graph = FALSE)

cds_3d_sub_DucT2_TOP2ACenter <- reduce_dimension(cds_sub_DucT2_TOP2ACenter, max_components = 3,preprocess_method = 'PCA')
plot_cells_3d(cds_3d_sub_DucT2_TOP2ACenter, genes = Main, show_trajectory_graph = FALSE)
plot_cells_3d(cds_3d_sub_DucT2_TOP2ACenter, color_cells_by="cell_cycle", show_trajectory_graph = FALSE)

cds_3d_sub_DucT2 <- reduce_dimension(cds_sub_DucT2, max_components = 3,preprocess_method = 'PCA')
plot_cells_3d(cds_3d_sub_DucT2, genes = Main, show_trajectory_graph = FALSE)
plot_cells_3d(cds_3d_sub_DucT2, color_cells_by="cell_cycle", show_trajectory_graph = FALSE)

 
