library("ggpubr")

# Function to prepare data for scatter plots
prepare_data_for_scatter <- function(gene1, gene2) {
  as.data.frame(cbind(t(gene1), t(gene2)))
}

# Function to create scatter plots
create_scatter_plot <- function(data, x, y, xlab, ylab) {
  ggscatter(data, x = x, y = y, 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "pearson", 
            xlab = xlab, ylab = ylab)
}

# Extracting markers and genes from cds_sub_AcinaDucT_NewK_ReCluster
markers <- list(
  PDAC = cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["PDAC"]],
  EMT = cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["EMT"]],
  NPC = cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["NPC"]],
  NE = cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["NE"]],
  ATR = cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["ATR"]],
  Migration = cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["Migration"]],
  Metastasis = cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["Metastasis"]],
  ACST = cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["ACST"]]
)

genes <- list(
  PTK2 = as.data.frame(cds_sub_AcinaDucT_NewK_ReCluster[rowData(cds_sub_AcinaDucT_NewK_ReCluster)$gene_short_name %in% c("PTK2"),]@assays@data@listData[["counts"]]),
  TOP2A = as.data.frame(cds_sub_AcinaDucT_NewK_ReCluster[rowData(cds_sub_AcinaDucT_NewK_ReCluster)$gene_short_name %in% c("TOP2A"),]@assays@data@listData[["counts"]]),
  CGAS = as.data.frame(cds_sub_AcinaDucT_NewK_ReCluster[rowData(cds_sub_AcinaDucT_NewK_ReCluster)$gene_short_name %in% c("CGAS"),]@assays@data@listData[["counts"]]),
  TP53 = as.data.frame(cds_sub_AcinaDucT_NewK_ReCluster[rowData(cds_sub_AcinaDucT_NewK_ReCluster)$gene_short_name %in% c("TP53"),]@assays@data@listData[["counts"]])
)

# Generate scatter plots
for (marker_name in names(markers)) {
  for (gene_name in names(genes)) {
    data <- prepare_data_for_scatter(markers[[marker_name]], genes[[gene_name]])
    plot <- create_scatter_plot(data, x = gene_name, y = marker_name, xlab = gene_name, ylab = marker_name)
    print(plot)
  }
}



