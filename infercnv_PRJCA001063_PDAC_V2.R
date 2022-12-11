# rm(list = ls()) #清除變數

library(infercnv)



############# Read file settings #############

PathName = setwd(getwd())
RVersion = "20210416V1"
dir.create(paste0(PathName,"/",RVersion))



##### (Test2V2)　#################################################################
cds_sub_CPV2 <- choose_cells(cds)
plot_cells(cds_sub_CPV2, label_cell_groups=FALSE, show_trajectory_graph = FALSE)
cds_sub_CP2V2<- cds_sub_CPV2[,!grepl("Macrophage cell", cds_sub_CPV2@colData@listData[["Cell_type"]], ignore.case=TRUE)]
cds_sub_CP2V2<- cds_sub_CP2V2[,!grepl("B cell", cds_sub_CP2V2@colData@listData[["Cell_type"]], ignore.case=TRUE)]
cds_sub_CP2V2<- cds_sub_CP2V2[,!grepl("T cell", cds_sub_CP2V2@colData@listData[["Cell_type"]], ignore.case=TRUE)]
cds_sub_CP2V2<- cds_sub_CP2V2[,!grepl("Fibroblast cell", cds_sub_CP2V2@colData@listData[["Cell_type"]], ignore.case=TRUE)]
cds_sub_CP2V2<- cds_sub_CP2V2[,!grepl("Stellate cell", cds_sub_CP2V2@colData@listData[["Cell_type"]], ignore.case=TRUE)]
cds_sub_CP2V2<- cds_sub_CP2V2[,!grepl("Endothelial cell", cds_sub_CP2V2@colData@listData[["Cell_type"]], ignore.case=TRUE)]
cds_sub_CP2V2<- cds_sub_CP2V2[,!grepl("Endocrine cell", cds_sub_CP2V2@colData@listData[["Cell_type"]], ignore.case=TRUE)]


### !!!各種細胞至少要有2個

#c("Acinar cell","Ductal cell type 1")
cds_sub_CPV2 <- cds_sub_CP2V2

plot_cells(cds_sub_CPV2, color_cells_by="partition", group_cells_by="partition",
           label_groups_by_cluster = F,show_trajectory_graph = F, group_label_size = 5)

cds_sub_CPV2_2 <- cds_sub_CPV2
cds_sub_CPV2_2 = cluster_cells(cds_sub_CPV2_2, resolution=1e-5)
plot_cells(cds_sub_CPV2_2, color_cells_by="cluster", group_cells_by="cluster",
           label_groups_by_cluster = F,show_trajectory_graph = F, group_label_size = 5)


colData(cds_sub_CPV2_2)$assigned_cell_type <- as.character(clusters(cds_sub_CPV2_2))
colData(cds_sub_CPV2_2)$assigned_cell_type = dplyr::recode(colData(cds_sub_CPV2_2)$assigned_cell_type, 
                                                "1"="Ductal cell type 1 S1",
                                                "3"="Acinar cell",
                                                "4"="Ductal cell type 1 S2",
                                                "10"="Ductal cell type 1 S3",
                                                
                                                "2"="Ductal cell type 2 S0",                                                
                                                
                                                "5"="Ductal cell type 2 S1",
                                                "6"="Ductal cell type 2 S2",
                                                "7"="Ductal cell type 2 S3",
                                                "8"="Ductal cell type 2 S4",
                                                "9"="Ductal cell type 2 S5",
                                                "11"="Ductal cell type 2 S6",
                                                
                                                "12"="Ductal cell type 2 S7",
                                                "13"="Ductal cell type 2 S8",
                                                "14"="Ductal cell type 2 S9",
                                                "15"="Ductal cell type 2 S10",
                                                "16"="Ductal cell type 2 S11",
                                                "17"="Ductal cell type 2 S12",
                                                "18"="Ductal cell type 2 S13",
                                                "19"="Ductal cell type 2 S14"
                                                
                                                )
plot_cells(cds_sub_CPV2_2, color_cells_by="assigned_cell_type", group_cells_by="cluster",
           label_cell_groups = F,show_trajectory_graph = F, group_label_size = 5)
plot_cells(cds_sub_CPV2_2, color_cells_by="assigned_cell_type", group_cells_by="cluster",
           label_cell_groups = T,show_trajectory_graph = F, group_label_size = 3)


tttExpression <- as.data.frame(cds_sub_CPV2_2@assays@data@listData[["logcounts"]])
# 可以考慮拿掉
tttExpression2 <- tttExpression[which(rowSums( tttExpression) > 0),]
tttExpression2 <- tttExpression


write.table(tttExpression2,file=paste0(PathName,"/",RVersion,"/",RVersion,"_tttExpression2.csv"),row.names = TRUE,sep = ',')
write.table(tttExpression2,file=paste0(PathName,"/",RVersion,"/",RVersion,"_tttExpression2.txt"),row.names = TRUE,sep = '\t')

tttExpression3<- as.matrix(tttExpression2)

# tttGene <- as.data.frame(row.names(tttExpression))
# row.names(tttGene) <- tttGene [,1]
tttCT <- as.data.frame(cds_sub_CPV2_2@colData)
tttCT2 <- tttCT[,-1]

library(dplyr)
# https://yijutseng.github.io/DataScienceRBook/eda.html
tttCT3 <- mutate(as.data.frame(tttCT[,9]),tttCT)
colnames(tttCT3)[1] <- colnames(tttCT3)[8]
row.names(tttCT3) <- row.names(tttCT)
write.table(tttCT3,file=paste0(PathName,"/",RVersion,"/",RVersion,"_tttCT3.csv"),row.names = TRUE,col.names = FALSE,sep = ',')
write.table(tttCT3,file=paste0(PathName,"/",RVersion,"/",RVersion,"_tttCT3.txt"),row.names = TRUE,col.names = FALSE,sep = '\t')
#write.table(tttCT3,file=paste0(PathName,"/",RVersion,"/",RVersion,"_tttCT3.csv"),row.names = TRUE,sep = ',')

tttCT4 <- as.data.frame(tttCT3[,1])

tttCT4 <-tttCT3[,1,drop = FALSE]

colnames(tttCT4)[1] <- 1
row.names(tttCT4)<- row.names(tttCT)

write.table(tttCT4,file=paste0(PathName,"/",RVersion,"/",RVersion,"_tttCT4.csv"),row.names = TRUE,col.names = FALSE,sep = ',')
write.table(tttCT4,file=paste0(PathName,"/",RVersion,"/",RVersion,"_tttCT4.txt"),row.names = TRUE,col.names = FALSE,sep = '\t')
#write.table(tttCT4,file=paste0(PathName,"/",RVersion,"/",RVersion,"_tttCT4.csv"),row.names = TRUE,sep = ',')


tttCT5<- as.matrix(tttCT4)



##
# create the infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix= tttExpression3,
                                    annotations_file= tttCT5,
                                    #        delim="\t",
                                    gene_order_file=system.file("extdata", "gencode_downsampled.EXAMPLE_ONLY_DONT_REUSE.txt", package = "infercnv"),
                                    ref_group_names=c("Acinar cell","Ductal cell type 1 S1","Ductal cell type 1 S2","Ductal cell type 1 S3"))
# no inference
# infercnv_obj@reference_grouped_cell_indices[["Acinar cell"]] <- as.integer(infercnv_obj@reference_grouped_cell_indices[["Acinar cell"]])

# #### !!!!### !!!各種細胞至少要有2個，但在這裡改會改不到總體表達的表格
# infercnv_obj@observation_grouped_cell_indices[["Macrophage cell"]]<- NULL
#  待處理改寫下列指令來刪除infercnv_obj中只有單一細胞的種類
# cds_sub_CP2<- cds_sub_CP[,!grepl("Macrophage cell", cds_sub_CP@colData@listData[["Cell_type"]], ignore.case=TRUE)]



# perform infercnv operations to reveal cnv signal
# infercnv_obj = infercnv::run(infercnv_obj,
#                              cutoff=1,  # use 1 for smart-seq, 0.1 for 10x-genomics
#                              out_dir="output_dir",  # dir is auto-created for storing outputs
#                              cluster_by_groups=T,   # cluster
#                              denoise=T,
#                              resume_mode = FALSE,
#                              HMM=T
# )

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             #  out_dir= "output_dir",
                             out_dir=tempfile(), 
                             cluster_by_groups=TRUE,
                             plot_steps=FALSE,
                             no_plot=FALSE,
                             denoise=TRUE,
                             resume_mode = FALSE,
                             HMM=TRUE)



