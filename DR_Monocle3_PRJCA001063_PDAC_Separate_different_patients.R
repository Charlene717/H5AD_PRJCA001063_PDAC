

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


##################  Grab specific terms ################## 
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


##################  AcinaDucT (ReCluster) ##################
cds_sub_AcinaDucT_NewK_ReCluster_Patient1 <- cds_sub_AcinaDucT_NewK_ReCluster[,grepl("T1", colData(cds_sub_AcinaDucT_NewK_ReCluster)$Patient, ignore.case=TRUE)]
plot_cells(cds_sub_AcinaDucT_NewK_ReCluster_Patient1 , color_cells_by="Cell_type", show_trajectory_graph = F,label_cell_groups = FALSE)
plot_cells(cds_sub_AcinaDucT_NewK_ReCluster_Patient1, genes=c(Main),cell_size=1,show_trajectory_graph = FALSE)

cds_sub_AcinaDucT_NewK_ReCluster_Patient1_subset <- cds_sub_AcinaDucT_NewK_ReCluster_Patient1 [rowData(cds_sub_AcinaDucT_NewK_ReCluster_Patient1 )$gene_short_name %in% ciliated_genes,]
plot_genes_violin(cds_sub_AcinaDucT_NewK_ReCluster_Patient1_subset, group_cells_by="Cell_type", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

plot_genes_violin(cds_sub_AcinaDucT_NewK_ReCluster_Patient1_subset, group_cells_by="ReCluster", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

