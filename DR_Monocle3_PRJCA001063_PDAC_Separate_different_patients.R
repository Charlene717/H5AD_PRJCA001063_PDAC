

ciliated_genes <- c("che-1",
                    "hlh-17",
                    "nhr-6",
                    "dmd-6",
                    "ceh-36",
                    "ham-1")


ciliated_genes <- c("TOP2A")

cds_subset <- cds[rowData(cds)$gene_short_name %in% ciliated_genes,]
plot_genes_violin(cds_subset, group_cells_by="Patient", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))


plot_genes_violin(cds_subset, group_cells_by="Cell_type", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))


##################  Grab specific terms in cds ################## 
## grepl Patient1
cds_Patient1 <- cds[,grepl("T1", colData(cds)$Patient, ignore.case=TRUE)]
plot_cells(cds_Patient1 , color_cells_by="Cell_type", show_trajectory_graph = F,label_cell_groups = FALSE)
plot_cells(cds_Patient1, genes=c(Main),cell_size=1,show_trajectory_graph = FALSE)

cds_Patient1_subset <- cds_Patient1 [rowData(cds_Patient1 )$gene_short_name %in% ciliated_genes,]
plot_genes_violin(cds_Patient1_subset, group_cells_by="Cell_type", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

## grepl Patient2
cds_Patient2 <- cds[,grepl("T2", colData(cds)$Patient, ignore.case=TRUE)]
plot_cells(cds_Patient2 , color_cells_by="Cell_type", show_trajectory_graph = F,label_cell_groups = FALSE)
plot_cells(cds_Patient2, genes=c(Main),cell_size=1,show_trajectory_graph = FALSE)

cds_Patient2_subset <- cds_Patient2[rowData(cds_Patient2)$gene_short_name %in% ciliated_genes,]
plot_genes_violin(cds_Patient2_subset, group_cells_by="Cell_type", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))


##################  Grab specific terms in AcinaDucT (ReCluster) ##################
plot_cells(cds_sub_AcinaDucT_NewK_ReCluster, color_cells_by="Patient",show_trajectory_graph = F,
           cell_size = 1.5,group_label_size = 4) 

cds_sub_AcinaDucT_NewK_ReCluster_subset <- cds_sub_AcinaDucT_NewK_ReCluster[rowData(cds_sub_AcinaDucT_NewK_ReCluster)$gene_short_name %in% ciliated_genes,]
plot_genes_violin(cds_sub_AcinaDucT_NewK_ReCluster_subset, group_cells_by="Cell_type", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

plot_genes_violin(cds_sub_AcinaDucT_NewK_ReCluster_subset, group_cells_by="ReCluster", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

# grab MDC00
cds_sub_AcinaDucT_NewK_ReCluster_MDC00 <- cds_sub_AcinaDucT_NewK_ReCluster[,grepl("MDC00", colData(cds_sub_AcinaDucT_NewK_ReCluster_subset)$ReCluster, ignore.case=TRUE)]
plot_cells(cds_sub_AcinaDucT_NewK_ReCluster_MDC00 , color_cells_by="cell_cycle", 
           show_trajectory_graph = F,label_cell_groups = FALSE,cell_size = 2)+ scale_color_manual(values = colors_cc)

plot_cells(cds_sub_AcinaDucT_NewK_ReCluster_MDC00 , genes = Main, 
           show_trajectory_graph = F,label_cell_groups = FALSE,cell_size = 2)

cds_sub_AcinaDucT_NewK_ReCluster_MDC00_subset <-cds_sub_AcinaDucT_NewK_ReCluster_MDC00[rowData(cds_sub_AcinaDucT_NewK_ReCluster_MDC00)$gene_short_name %in% ciliated_genes,]

plot_genes_violin(cds_sub_AcinaDucT_NewK_ReCluster_subset_MDC00, group_cells_by="cell_cycle", ncol=2, log_scale = FALSE)+ scale_fill_manual(values = colors_cc)
plot_genes_violin(cds_sub_AcinaDucT_NewK_ReCluster_subset_MDC00, group_cells_by="cell_cycle", ncol=2, log_scale = T)+ scale_fill_manual(values = colors_cc)

###### Convert Monocle3 Object to Seurat Object ######
getFilePath("Monocle3_To_Seurat.R")
getFilePath("Monocle3_To_Seurat2.R")
marrow_sub_AcinaDucT_NewK_ReCluster_MDC00 <- Monocle3_To_Seurat2(cds_sub_AcinaDucT_NewK_ReCluster_MDC00,"cds_sub_AcinaDucT_NewK_ReCluster_subset_MDC00") #這個function存在於Monocle3_To_Seurat.R裡面

###### Insert the cell cycle results from Monocle3 cds_sub into the  Seurat  marrow_sub ######
marrow_sub_AcinaDucT_NewK_ReCluster_MDC00@active.ident <- cds_sub_AcinaDucT_NewK_ReCluster_subset_MDC00@colData@listData$cell_cycle
RidgePlot(marrow_sub_AcinaDucT_NewK_ReCluster_MDC00,cols = colors_cc, features = c(Main), ncol = 1,log = T)


# grab MD00
cds_sub_AcinaDucT_NewK_ReCluster_subset_MD00 <- cds_sub_AcinaDucT_NewK_ReCluster_subset[,grepl("MD00", colData(cds_sub_AcinaDucT_NewK_ReCluster_subset)$ReCluster, ignore.case=TRUE)]
plot_cells(cds_sub_AcinaDucT_NewK_ReCluster_subset_MD00 , color_cells_by="cell_cycle", 
           show_trajectory_graph = F,label_cell_groups = FALSE,cell_size = 2)+ scale_color_manual(values = colors_cc)

plot_cells(cds_sub_AcinaDucT_NewK_ReCluster_subset_MD00 , genes = Main, 
           show_trajectory_graph = F,label_cell_groups = FALSE,cell_size = 2)

plot_genes_violin(cds_sub_AcinaDucT_NewK_ReCluster_subset_MD00, group_cells_by="cell_cycle", ncol=2, log_scale = FALSE)+ scale_fill_manual(values = colors_cc)


#Patient1
cds_sub_AcinaDucT_NewK_ReCluster_Patient1 <- cds_sub_AcinaDucT_NewK_ReCluster[,grepl("T1", colData(cds_sub_AcinaDucT_NewK_ReCluster)$Patient, ignore.case=TRUE)]
plot_cells(cds_sub_AcinaDucT_NewK_ReCluster_Patient1 , color_cells_by="Cell_type", show_trajectory_graph = F,label_cell_groups = FALSE)
plot_cells(cds_sub_AcinaDucT_NewK_ReCluster_Patient1, genes=c(Main),cell_size=1,show_trajectory_graph = FALSE)

cds_sub_AcinaDucT_NewK_ReCluster_Patient1_subset <- cds_sub_AcinaDucT_NewK_ReCluster_Patient1 [rowData(cds_sub_AcinaDucT_NewK_ReCluster_Patient1 )$gene_short_name %in% ciliated_genes,]
plot_genes_violin(cds_sub_AcinaDucT_NewK_ReCluster_Patient1_subset, group_cells_by="Cell_type", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

plot_genes_violin(cds_sub_AcinaDucT_NewK_ReCluster_Patient1_subset, group_cells_by="ReCluster", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))


########################  Try for loop for Separating differentpatients ##########################
##################  Grab specific terms in cds ################## 
for (i in c(1:24)) {
  
  ## grepl Patienti
  cds_Patienti <- cds[,grepl(paste0("T",i), colData(cds)$Patient, ignore.case=TRUE)]
  # 
  # ## Plot error
  # png(paste0(PathName,"/",RVersion,"/",RVersion,"_","cds_Patienti",i,"UMAP_Cell_type.png"), width = 640, height = 360) # 設定輸出圖檔
  # plot_cells(cds_Patienti , color_cells_by="Cell_type", show_trajectory_graph = F,label_cell_groups = FALSE)
  # dev.off() # 關閉輸出圖檔
  # 
  # png(paste0(PathName,"/",RVersion,"/",RVersion,"_","cds_Patienti",i,"_",ciliated_genes,"_UMAP.png"), width = 640, height = 360) # 設定輸出圖檔
  # plot_cells(cds_Patienti, genes=c(Main),cell_size=1,show_trajectory_graph = FALSE)
  # dev.off() # 關閉輸出圖檔
  
  cds_Patienti_subset <- cds_Patienti[rowData(cds_Patienti)$gene_short_name %in% ciliated_genes,]
  
  # png(paste0(PathName,"/",RVersion,"/",RVersion,"_","cds_Patienti",i,"_",ciliated_genes,"_UMAP_violin_type.png"), width = 640, height = 360) # 設定輸出圖檔
  # plot_genes_violin(cds_Patienti_subset, group_cells_by="Cell_type", ncol=2) +
  #   theme(axis.text.x=element_text(angle=45, hjust=1))
  # dev.off() # 關閉輸出圖檔
  
  assign(paste0("cds_Patient", i),cds_Patienti)
  assign(paste0("cds_Patient", i,"_subset"),cds_Patienti_subset)
#  rm(cds_Patienti,cds_Patienti_subset)
}

plot_genes_violin(cds_Patient1_subset, group_cells_by="Cell_type", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))


png(paste0(PathName,"/",RVersion,"/",RVersion,"_","cds_Patienti",i,"UMAP_Cell_type.png"), width = 640, height = 360) # 設定輸出圖檔
plot_cells(cds_Patienti , color_cells_by="Cell_type", show_trajectory_graph = F,label_cell_groups = FALSE)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","cds_Patienti",i,"_",ciliated_genes,"_UMAP.png"), width = 640, height = 360) # 設定輸出圖檔
plot_cells(cds_Patienti, genes=c(Main),cell_size=1,show_trajectory_graph = FALSE)
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_","cds_Patienti",i,"_",ciliated_genes,"_UMAP_violin_type.png"), width = 640, height = 360) # 設定輸出圖檔
plot_genes_violin(cds_Patienti_subset, group_cells_by="Cell_type", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
dev.off() # 關閉輸出圖檔


##################  Grab specific terms in AcinaDucT (ReCluster) ################## 
for (i in c(1:24)) {
  ## grepl Patienti
  cds_sub_AcinaDucT_NewK_ReCluster_Patienti <- cds_sub_AcinaDucT_NewK_ReCluster[,grepl(paste0("T",i), colData(cds_sub_AcinaDucT_NewK_ReCluster)$Patient, ignore.case=TRUE)]
  plot_cells(cds_sub_AcinaDucT_NewK_ReCluster_Patienti , color_cells_by="Cell_type", show_trajectory_graph = F,label_cell_groups = FALSE)
  plot_cells(cds_sub_AcinaDucT_NewK_ReCluster_Patienti, genes=c(Main),cell_size=1,show_trajectory_graph = FALSE)
  cds_sub_AcinaDucT_NewK_ReCluster_Patienti_subset <- cds_sub_AcinaDucT_NewK_ReCluster_Patienti[rowData(cds_sub_AcinaDucT_NewK_ReCluster_Patienti)$gene_short_name %in% ciliated_genes,]
  plot_genes_violin(cds_sub_AcinaDucT_NewK_ReCluster_Patienti_subset, group_cells_by="Cell_type", ncol=2) +
    theme(axis.text.x=element_text(angle=45, hjust=1))
  
  plot_genes_violin(cds_sub_AcinaDucT_NewK_ReCluster_Patienti_subset, group_cells_by="ReCluster", ncol=2) +
    theme(axis.text.x=element_text(angle=45, hjust=1))
  
  assign(paste0("cds_sub_AcinaDucT_NewK_ReCluster_Patient", i),cds_sub_AcinaDucT_NewK_ReCluster_Patienti)
  assign(paste0("cds_sub_AcinaDucT_NewK_ReCluster_Patient", i,"_subset"),cds_sub_AcinaDucT_NewK_ReCluster_Patienti_subset)
  
  #  rm(cds_sub_AcinaDucT_NewK_ReCluster_Patienti,cds_sub_AcinaDucT_NewK_ReCluster_Patienti_subset)
}

plot_genes_violin(cds_sub_AcinaDucT_NewK_ReCluster_Patient1_subset, group_cells_by="ReCluster", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

plot_cells(cds_sub_AcinaDucT_NewK_ReCluster_Patient1_subset, color_cells_by ="ReCluster",
           show_trajectory_graph = F,cell_size = 2,group_label_size = 4) 
