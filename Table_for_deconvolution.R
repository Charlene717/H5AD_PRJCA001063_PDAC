## Table for deconvolution

plot_cells(cds_sub_AcinaDucT_NewK_ReCluster, color_cells_by = "ReCluster",cell_size=2, 
           label_cell_groups=TRUE, show_trajectory_graph = FALSE, group_label_size =4)

cds@colData@listData[["ReCluster"]] <- as.character(cds@colData@listData[["Cell_type"]])
plot_cells(cds, color_cells_by = "ReCluster",cell_size=2, 
           label_cell_groups=TRUE, show_trajectory_graph = FALSE, group_label_size =4)

colData(cds)[colnames(cds_sub_AcinaDucT_NewK_ReCluster),]$ReCluster <- colData(cds_sub_AcinaDucT_NewK_ReCluster)$ReCluster 
plot_cells(cds, color_cells_by = "ReCluster",cell_size=2, 
           label_cell_groups=TRUE, show_trajectory_graph = FALSE, group_label_size =4)
plot_cells(cds, color_cells_by = "ReCluster",cell_size=2, 
           label_cell_groups=F, show_trajectory_graph = FALSE, group_label_size =4)


## Gene expression table
GeneExpMatrix <- as.data.frame(cds@assays@data@listData[["counts"]])
rownames(GeneExpMatrix)[[4137]] <- c("STR") # Unknown error record "<class 'str'>"
rownames(GeneExpMatrix)[[4818]] <- c("STR1") # Unknown error record "<class 'str'-1>"
rownames(GeneExpMatrix)[[5896]] <- c("STR2") # Unknown error record "<class 'str'-2>"
rownames(GeneExpMatrix)[[5897]] <- c("STR3") # Unknown error record "<class 'str'-3>"
rownames(GeneExpMatrix)[[15702]] <- c("STR4") # Unknown error record "<class 'str'-4>"
rownames(GeneExpMatrix)[[15860]] <- c("STR5") # Unknown error record "<class 'str'-5>"
rownames(GeneExpMatrix)[[16402]] <- c("STR6") # Unknown error record "<class 'str'-6>"
rownames(GeneExpMatrix)[[16863]] <- c("STR7") # Unknown error record "<class 'str'-7>"
rownames(GeneExpMatrix)[[17574]] <- c("STR8") # Unknown error record "<class 'str'-8>"


write.table(GeneExpMatrix, file = paste0(PathName,"/Deconvolution/scRNA/FCortex_PRJCA001063.txt"), quote = FALSE, 
            row.names = TRUE, col.names = TRUE, sep = " ")
# Error in paste(col.names, collapse = sep) : 
#   could not allocate memory (1 Mb) in C function 'R_AllocStringBuffer'


## Cell_type table
ReCluster <- as.data.frame(colData(cds)$ReCluster)
colnames(ReCluster) <- c("ReCluster")
Index <- as.data.frame(cds@colData@listData[["CELL"]])
colnames(Index) <- c("Index")

CellType <- cbind(Index,ReCluster)
write.table(CellType, file = paste0(PathName,"/Deconvolution/scRNA/CellType_PRJCA001063.txt"), quote = FALSE, 
            row.names = FALSE, col.names = TRUE, sep = "\t")


###################################################################
