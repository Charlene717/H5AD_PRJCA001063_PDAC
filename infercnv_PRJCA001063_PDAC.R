# rm(list = ls()) #清除變數

library(infercnv)



############# Read file settings #############

PathName = setwd(getwd())
RVersion = "20210415V1"
dir.create(paste0(PathName,"/",RVersion))



####### (Test All)　###############################################################
tttExpression <- as.data.frame(cds@assays@data@listData[["logcounts"]])
# tttGene <- as.data.frame(row.names(tttExpression))
# row.names(tttGene) <- tttGene [,1]
tttCT <- as.data.frame(cds@colData)
tttCT2 <- tttCT[,-1]

library(dplyr)
# https://yijutseng.github.io/DataScienceRBook/eda.html
tttCT3 <- mutate(as.data.frame(tttCT[,5]),tttCT)
colnames(tttCT3)[1] <- colnames(tttCT3)[6]
row.names(tttCT3) <- row.names(tttCT)


# create the infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix= tttExpression,
                                    annotations_file= tttCT3,
                                    delim="\t",
                                    gene_order_file=system.file("extdata", "gencode_downsampled.EXAMPLE_ONLY_DONT_REUSE.txt", package = "infercnv"),
                                    ref_group_names=c("Acinar cell","Ductal cell type 1"))

# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="output_dir",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T
)

####### (Test１)　###############################################################
plot_cells(cds_sub_DucT2, label_cell_groups=FALSE, show_trajectory_graph = FALSE)


tttExpression <- as.data.frame(cds_sub_DucT2@assays@data@listData[["logcounts"]])
# tttGene <- as.data.frame(row.names(tttExpression))
# row.names(tttGene) <- tttGene [,1]
tttCT <- as.data.frame(cds_sub_DucT2@colData)
tttCT2 <- tttCT[,-1]

library(dplyr)
# https://yijutseng.github.io/DataScienceRBook/eda.html
tttCT3 <- mutate(as.data.frame(tttCT[,5]),tttCT)
colnames(tttCT3)[1] <- colnames(tttCT3)[6]
row.names(tttCT3) <- row.names(tttCT)


# create the infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix= tttExpression,
                                    annotations_file= tttCT3,
                                    delim="\t",
                                    gene_order_file=system.file("extdata", "gencode_downsampled.EXAMPLE_ONLY_DONT_REUSE.txt", package = "infercnv"),
                                    ref_group_names=c("Acinar cell","Ductal cell type 1"))

# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="output_dir",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T
)

##### (Test2)　#################################################################
cds_sub_CP <- choose_cells(cds)
plot_cells(cds_sub_CP, label_cell_groups=FALSE, show_trajectory_graph = FALSE)
cds_sub_CP2<- cds_sub_CP[,!grepl("Macrophage cell", cds_sub_CP@colData@listData[["Cell_type"]], ignore.case=TRUE)]
### !!!各種細胞至少要有2個
cds_sub_CP <- cds_sub_CP2
#c("Acinar cell","Ductal cell type 1")

tttExpression <- as.data.frame(cds_sub_CP@assays@data@listData[["logcounts"]])
# 可以考慮拿掉
tttExpression2 <- tttExpression[which(rowSums( tttExpression) > 0),]

write.table(tttExpression2,file=paste0(PathName,"/",RVersion,"/",RVersion,"_tttExpression2.csv"),row.names = TRUE,sep = ',')
write.table(tttExpression2,file=paste0(PathName,"/",RVersion,"/",RVersion,"_tttExpression2.txt"),row.names = TRUE,sep = '\t')

tttExpression3<- as.matrix(tttExpression2)

# tttGene <- as.data.frame(row.names(tttExpression))
# row.names(tttGene) <- tttGene [,1]
tttCT <- as.data.frame(cds_sub_CP@colData)
tttCT2 <- tttCT[,-1]

library(dplyr)
# https://yijutseng.github.io/DataScienceRBook/eda.html
tttCT3 <- mutate(as.data.frame(tttCT[,5]),tttCT)
colnames(tttCT3)[1] <- colnames(tttCT3)[6]
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




# create the infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix= paste0(PathName,"/",RVersion,"/",RVersion,"_tttExpression2.txt"),
                                    annotations_file= paste0(PathName,"/",RVersion,"/",RVersion,"_tttCT4.txt"),
                                    delim="\t",
                                    gene_order_file=system.file("extdata", "gencode_downsampled.EXAMPLE_ONLY_DONT_REUSE.txt", package = "infercnv"),
                                    ref_group_names=c("Acinar cell","Ductal cell type 1"))

# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                         #   out_dir= "output_dir",
                             out_dir=tempfile(), 
                             cluster_by_groups=TRUE,
                             plot_steps=TRUE,
                             no_plot=FALSE,
                             denoise=TRUE,
                             resume_mode = FALSE,
                             HMM=TRUE)

##
# create the infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix= tttExpression3,
                                    annotations_file= tttCT5,
                           #        delim="\t",
                                    gene_order_file=system.file("extdata", "gencode_downsampled.EXAMPLE_ONLY_DONT_REUSE.txt", package = "infercnv"),
                                    ref_group_names=c("Acinar cell","Ductal cell type 1"))
# infercnv_obj@reference_grouped_cell_indices[["Acinar cell"]] <- as.integer(infercnv_obj@reference_grouped_cell_indices[["Acinar cell"]])
# infercnv_obj@reference_grouped_cell_indices[["Ductal cell type 1"]] <- as.integer(infercnv_obj@reference_grouped_cell_indices[["Ductal cell type 1"]])
# infercnv_obj@observation_grouped_cell_indices[["Ductal cell type 2"]] <- as.integer(infercnv_obj@observation_grouped_cell_indices[["Ductal cell type 2"]])
# infercnv_obj@observation_grouped_cell_indices[["T cell"]] <- as.integer(infercnv_obj@observation_grouped_cell_indices[["T cell"]] )
# infercnv_obj@observation_grouped_cell_indices[["Macrophage cell"]] <- as.integer(infercnv_obj@observation_grouped_cell_indices[["Macrophage cell"]])
# infercnv_obj@observation_grouped_cell_indices[["Endothelial cell"]] <- as.integer(infercnv_obj@observation_grouped_cell_indices[["Endothelial cell"]])
# infercnv_obj@observation_grouped_cell_indices[["Fibroblast cell"]] <- as.integer(infercnv_obj@observation_grouped_cell_indices[["Fibroblast cell"]] )
# infercnv_obj@observation_grouped_cell_indices[["Stellate cell"]] <- as.integer(infercnv_obj@observation_grouped_cell_indices[["Stellate cell"]] )
# infercnv_obj@observation_grouped_cell_indices[["Endocrine cell"]] <- as.integer(infercnv_obj@observation_grouped_cell_indices[["Endocrine cell"]])

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


# #########
# infercnv_obj = CreateInfercnvObject(raw_counts_matrix=system.file("extdata", "20210415V1_tttExpression2.txt", package = "infercnv"),
#                                     annotations_file=system.file("extdata", "20210415V1_tttCT3.txt", package = "infercnv"),
#                                     delim="\t",
#                                     gene_order_file=system.file("extdata", "gencode_downsampled.EXAMPLE_ONLY_DONT_REUSE.txt", package = "infercnv"),
#                                     ref_group_names=c("Microglia/Macrophage","Oligodendrocytes (non-malignant)")) 
# infercnv_obj = infercnv::run(infercnv_obj,
#                              cutoff=1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
#                              out_dir=tempfile(), 
#                              cluster_by_groups=TRUE,
#                              plot_steps=FALSE,
#                              no_plot=FALSE,
#                              denoise=TRUE,
#                              HMM=TRUE)
# 
