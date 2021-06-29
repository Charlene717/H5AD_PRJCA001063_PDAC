memory.limit(150000)
# install.packages("ggpubr")
library("ggpubr")


## Markers
PDAC <- cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["PDAC"]]
EMT <- cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["EMT"]]
NPC <- cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["NPC"]]
NE <- cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["NE"]]
ATR <- cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["ATR"]]
Migration <- cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["Migration"]]
Metastasis <- cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["Metastasis"]]
ACST <- cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["ACST"]]

##
ReCluster <- as.character(cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["ReCluster"]])
Cell_cycle <- as.character(cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["cell_cycle"]])

## Genes
# PTK2
ciliated_genes_PTK2 <- c("PTK2")
cds_sub_AcinaDucT_NewK_ReCluster_PTK2 <- cds_sub_AcinaDucT_NewK_ReCluster[rowData(cds_sub_AcinaDucT_NewK_ReCluster)$gene_short_name %in% ciliated_genes_PTK2,]
plot_genes_violin(cds_sub_AcinaDucT_NewK_ReCluster_PTK2, group_cells_by="cell_cycle", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
PTK2 <- as.data.frame(cds_sub_AcinaDucT_NewK_ReCluster_PTK2@assays@data@listData[["counts"]])
# PTK2 <- as.data.frame(cds_sub_AcinaDucT_NewK_ReCluster_PTK2@assays@data@listData[["logcounts"]])

# TOP2A
ciliated_genes_TOP2A <- c("TOP2A")
cds_sub_AcinaDucT_NewK_ReCluster_TOP2A <- cds_sub_AcinaDucT_NewK_ReCluster[rowData(cds_sub_AcinaDucT_NewK_ReCluster)$gene_short_name %in% ciliated_genes_TOP2A,]
plot_genes_violin(cds_sub_AcinaDucT_NewK_ReCluster_TOP2A, group_cells_by="cell_cycle", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
TOP2A <- as.data.frame(cds_sub_AcinaDucT_NewK_ReCluster_TOP2A@assays@data@listData[["counts"]])
# TOP2A <- as.data.frame(cds_sub_AcinaDucT_NewK_ReCluster_TOP2A@assays@data@listData[["logcounts"]])

# CGAS
ciliated_genes_CGAS <- c("CGAS")
cds_sub_AcinaDucT_NewK_ReCluster_CGAS <- cds_sub_AcinaDucT_NewK_ReCluster[rowData(cds_sub_AcinaDucT_NewK_ReCluster)$gene_short_name %in% ciliated_genes_CGAS,]
plot_genes_violin(cds_sub_AcinaDucT_NewK_ReCluster_CGAS, group_cells_by="cell_cycle", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
CGAS <- as.data.frame(cds_sub_AcinaDucT_NewK_ReCluster_CGAS@assays@data@listData[["counts"]])

# TP53
ciliated_genes_TP53 <- c("TP53")
cds_sub_AcinaDucT_NewK_ReCluster_TP53 <- cds_sub_AcinaDucT_NewK_ReCluster[rowData(cds_sub_AcinaDucT_NewK_ReCluster)$gene_short_name %in% ciliated_genes_TP53,]
plot_genes_violin(cds_sub_AcinaDucT_NewK_ReCluster_TP53, group_cells_by="cell_cycle", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
TP53 <- as.data.frame(cds_sub_AcinaDucT_NewK_ReCluster_TP53@assays@data@listData[["counts"]])

## \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

library("ggpubr")

## NPC
## NPC vs PTK2
NP_Sum_PTK2 <- as.data.frame(cbind(NPC,t(PTK2)))
ggscatter(NP_Sum_PTK2, x = "PTK2", y = "NPC", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "PTK2", ylab = "NPC")

## NPC vs TOP2A
NP_Sum_TOP2A <- as.data.frame(cbind(NPC,t(TOP2A)))
ggscatter(NP_Sum_TOP2A, x = "TOP2A", y = "NPC", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "TOP2A", ylab = "NPC")

## NPC vs CGAS
NP_Sum_TOP2A <- as.data.frame(cbind(NPC,t(CGAS)))
ggscatter(NP_Sum_TOP2A, x = "CGAS", y = "NPC", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",  #"spearman"
          xlab = "CGAS", ylab = "NPC")

## NE
## NE vs PTK2
NP_Sum_PTK2 <- as.data.frame(cbind(NE,t(PTK2)))
ggscatter(NP_Sum_PTK2, x = "PTK2", y = "NE", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "PTK2", ylab = "NE")

## NE vs TOP2A
NP_Sum_TOP2A <- as.data.frame(cbind(NE,t(TOP2A)))
ggscatter(NP_Sum_TOP2A, x = "TOP2A", y = "NE", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "TOP2A", ylab = "NE")

## NE vs CGAS
NP_Sum_TOP2A <- as.data.frame(cbind(NE,t(CGAS)))
ggscatter(NP_Sum_TOP2A, x = "CGAS", y = "NE", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",  #"spearman"
          xlab = "CGAS", ylab = "NE")

## ACST
## ACST vs PTK2
NP_Sum_PTK2 <- as.data.frame(cbind(ACST,t(PTK2)))
ggscatter(NP_Sum_PTK2, x = "PTK2", y = "ACST", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "PTK2", ylab = "ACST")

## ACST vs TOP2A
NP_Sum_TOP2A <- as.data.frame(cbind(ACST,t(TOP2A)))
ggscatter(NP_Sum_TOP2A, x = "TOP2A", y = "ACST", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "TOP2A", ylab = "ACST")

## ACST vs CGAS
NP_Sum_TOP2A <- as.data.frame(cbind(ACST,t(CGAS)))
ggscatter(NP_Sum_TOP2A, x = "CGAS", y = "ACST", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",  #"spearman"
          xlab = "CGAS", ylab = "ACST")

## \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
######################################  cds_subset (CoreCD00) ########################################
## grepl CoreCD00
cds_sub_CoreCD00 <- cds_sub_AcinaDucT_NewK_ReCluster[,grepl("CoreCD00", colData(cds_sub_AcinaDucT_NewK_ReCluster)$ReCluster, ignore.case=TRUE)]
plot_cells(cds_sub_CoreCD00, color_cells_by="ReCluster", show_trajectory_graph = F,label_cell_groups = F,cell_size = 2)

##
## Markers
CoreCD00_PDAC <- cds_sub_CoreCD00@colData@listData[["PDAC"]]
CoreCD00_EMT <- cds_sub_CoreCD00@colData@listData[["EMT"]]
CoreCD00_NPC <- cds_sub_CoreCD00@colData@listData[["NPC"]]
CoreCD00_NE <- cds_sub_CoreCD00@colData@listData[["NE"]]
CoreCD00_ATR <- cds_sub_CoreCD00@colData@listData[["ATR"]]
CoreCD00_Migration <- cds_sub_CoreCD00@colData@listData[["Migration"]]
CoreCD00_Metastasis <- cds_sub_CoreCD00@colData@listData[["Metastasis"]]
CoreCD00_ACST <- cds_sub_CoreCD00@colData@listData[["ACST"]]

##
CoreCD00_ReCluster <- as.character(cds_sub_CoreCD00@colData@listData[["ReCluster"]])
CoreCD00_Cell_cycle <- as.character(cds_sub_CoreCD00@colData@listData[["cell_cycle"]])

## Genes
# PTK2
ciliated_genes_PTK2_CoreCD00 <- c("PTK2")
cds_sub_CoreCD00_PTK2 <- cds_sub_CoreCD00[rowData(cds_sub_CoreCD00)$gene_short_name %in% ciliated_genes_PTK2_CoreCD00,]
plot_genes_violin(cds_sub_CoreCD00_PTK2, group_cells_by="cell_cycle", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
CoreCD00_PTK2 <- as.data.frame(cds_sub_CoreCD00_PTK2@assays@data@listData[["counts"]])
# PTK2 <- as.data.frame(cds_sub_AcinaDucT_NewK_ReCluster_PTK2@assays@data@listData[["logcounts"]])

# TOP2A
ciliated_genes_TOP2A_CoreCD00 <- c("TOP2A")
cds_sub_CoreCD00_TOP2A <- cds_sub_CoreCD00[rowData(cds_sub_CoreCD00)$gene_short_name %in% ciliated_genes_TOP2A_CoreCD00,]
plot_genes_violin(cds_sub_CoreCD00_TOP2A, group_cells_by="cell_cycle", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
CoreCD00_TOP2A <- as.data.frame(cds_sub_CoreCD00_TOP2A@assays@data@listData[["counts"]])
# TOP2A <- as.data.frame(cds_sub_AcinaDucT_NewK_ReCluster_TOP2A@assays@data@listData[["logcounts"]])

# CGAS
ciliated_genes_CGAS_CoreCD00 <- c("CGAS")
cds_sub_CoreCD00_CGAS <- cds_sub_CoreCD00[rowData(cds_sub_CoreCD00)$gene_short_name %in% ciliated_genes_CGAS_CoreCD00,]
plot_genes_violin(cds_sub_CoreCD00_CGAS, group_cells_by="cell_cycle", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
CoreCD00_CGAS <- as.data.frame(cds_sub_CoreCD00_CGAS@assays@data@listData[["counts"]])

##--------------------------------------------------------------------------------------------------------------

## PTK2 vs CGAS
PTK2_Sum_CGAS_CoreCD00 <- as.data.frame(cbind(t(CoreCD00_CGAS),t(CoreCD00_PTK2)))
ggscatter(PTK2_Sum_CGAS_CoreCD00, x = "PTK2", y = "CGAS", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",  #"spearman"
          xlab = "PTK2", ylab = "CGAS")
## PTK2 vs TOP2A
PTK2_Sum_TOP2A_CoreCD00 <- as.data.frame(cbind(t(CoreCD00_TOP2A),t(CoreCD00_PTK2)))
ggscatter(PTK2_Sum_TOP2A_CoreCD00, x = "PTK2", y = "TOP2A", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",  #"spearman"
          xlab = "PTK2", ylab = "TOP2A")

## CGAS vs TOP2A
CGAS_Sum_TOP2A_CoreCD00 <- as.data.frame(cbind(t(CoreCD00_TOP2A),t(CoreCD00_CGAS)))
ggscatter(CGAS_Sum_TOP2A_CoreCD00, x = "CGAS", y = "TOP2A", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",  #"spearman"
          xlab = "CGAS", ylab = "TOP2A")


## NPC
## NPC vs PTK2
NP_Sum_PTK2_CoreCD00 <- as.data.frame(cbind(CoreCD00_NPC,t(CoreCD00_PTK2)))
ggscatter(NP_Sum_PTK2_CoreCD00, x = "PTK2", y = "CoreCD00_NPC", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",  #"spearman"
          xlab = "PTK2", ylab = "NPC")

## NPC vs TOP2A
NP_Sum_TOP2A_CoreCD00 <- as.data.frame(cbind(CoreCD00_NPC,t(CoreCD00_TOP2A)))
ggscatter(NP_Sum_TOP2A_CoreCD00, x = "TOP2A", y = "CoreCD00_NPC", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "TOP2A", ylab = "NPC")

## NPC vs CGAS
NP_Sum_CGAS_CoreCD00  <- as.data.frame(cbind(CoreCD00_NPC,t(CoreCD00_CGAS)))
ggscatter(NP_Sum_CGAS_CoreCD00, x = "CGAS", y = "CoreCD00_NPC", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",  #"spearman"
          xlab = "CGAS", ylab = "NPC")

## NE
## NE vs PTK2
NE_Sum_PTK2_CoreCD00 <- as.data.frame(cbind(CoreCD00_NE,t(CoreCD00_PTK2)))
ggscatter(NE_Sum_PTK2_CoreCD00, x = "PTK2", y = "CoreCD00_NE", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",  #"spearman"
          xlab = "PTK2", ylab = "NE")

## NE vs TOP2A
NE_Sum_TOP2A_CoreCD00 <- as.data.frame(cbind(CoreCD00_NE,t(CoreCD00_TOP2A)))
ggscatter(NE_Sum_TOP2A_CoreCD00, x = "TOP2A", y = "CoreCD00_NE", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "TOP2A", ylab = "NE")

## NE vs CGAS
NE_Sum_CGAS_CoreCD00  <- as.data.frame(cbind(CoreCD00_NE,t(CoreCD00_CGAS)))
ggscatter(NE_Sum_CGAS_CoreCD00, x = "CGAS", y = "CoreCD00_NE", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",  #"spearman"
          xlab = "CGAS", ylab = "NE")

## ACST
## ACST vs PTK2
ACST_Sum_PTK2_CoreCD00 <- as.data.frame(cbind(CoreCD00_ACST,t(CoreCD00_PTK2)))
ggscatter(ACST_Sum_PTK2_CoreCD00, x = "PTK2", y = "CoreCD00_ACST", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",  #"spearman"
          xlab = "PTK2", ylab = "ACST")

## ACST vs TOP2A
ACST_Sum_TOP2A_CoreCD00 <- as.data.frame(cbind(CoreCD00_ACST,t(CoreCD00_TOP2A)))
ggscatter(ACST_Sum_TOP2A_CoreCD00, x = "TOP2A", y = "CoreCD00_ACST", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "TOP2A", ylab = "ACST")

## ACST vs CGAS
ACST_Sum_CGAS_CoreCD00  <- as.data.frame(cbind(CoreCD00_ACST,t(CoreCD00_CGAS)))
ggscatter(ACST_Sum_CGAS_CoreCD00, x = "CGAS", y = "CoreCD00_ACST", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",  #"spearman"
          xlab = "CGAS", ylab = "ACST")


## \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

######################################  cds_subset (DucT2_TOP2A_SmallCTR)########################################
cds_sub_DucT2_TOP2A_SmallCTR <- choose_cells(cds_sub_AcinaDucT_NewK_ReCluster)
plot_cells(cds_sub_DucT2_TOP2A_SmallCTR, genes=c("TOP2A"),cell_size=2, label_cell_groups=FALSE, show_trajectory_graph = FALSE)

##
## Markers
DucT2_TOP2A_SmallCTR_PDAC <- cds_sub_DucT2_TOP2A_SmallCTR@colData@listData[["PDAC"]]
DucT2_TOP2A_SmallCTR_EMT <- cds_sub_DucT2_TOP2A_SmallCTR@colData@listData[["EMT"]]
DucT2_TOP2A_SmallCTR_NPC <- cds_sub_DucT2_TOP2A_SmallCTR@colData@listData[["NPC"]]
DucT2_TOP2A_SmallCTR_NE <- cds_sub_DucT2_TOP2A_SmallCTR@colData@listData[["NE"]]
DucT2_TOP2A_SmallCTR_ATR <- cds_sub_DucT2_TOP2A_SmallCTR@colData@listData[["ATR"]]
DucT2_TOP2A_SmallCTR_Migration <- cds_sub_DucT2_TOP2A_SmallCTR@colData@listData[["Migration"]]
DucT2_TOP2A_SmallCTR_Metastasis <- cds_sub_DucT2_TOP2A_SmallCTR@colData@listData[["Metastasis"]]
DucT2_TOP2A_SmallCTR_ACST <- cds_sub_DucT2_TOP2A_SmallCTR@colData@listData[["ACST"]]

##
DucT2_TOP2A_SmallCTR_ReCluster <- as.character(cds_sub_DucT2_TOP2A_SmallCTR@colData@listData[["ReCluster"]])
DucT2_TOP2A_SmallCTR_Cell_cycle <- as.character(cds_sub_DucT2_TOP2A_SmallCTR@colData@listData[["cell_cycle"]])

## Genes
# PTK2
ciliated_genes_PTK2_DucT2_TOP2A_SmallCTR <- c("PTK2")
cds_sub_DucT2_TOP2A_SmallCTR_PTK2 <- cds_sub_DucT2_TOP2A_SmallCTR[rowData(cds_sub_DucT2_TOP2A_SmallCTR)$gene_short_name %in% ciliated_genes_PTK2_DucT2_TOP2A_SmallCTR,]
plot_genes_violin(cds_sub_DucT2_TOP2A_SmallCTR_PTK2, group_cells_by="cell_cycle", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
DucT2_TOP2A_SmallCTR_PTK2 <- as.data.frame(cds_sub_DucT2_TOP2A_SmallCTR_PTK2@assays@data@listData[["counts"]])
# PTK2 <- as.data.frame(cds_sub_AcinaDucT_NewK_ReCluster_PTK2@assays@data@listData[["logcounts"]])

# TOP2A
ciliated_genes_TOP2A_DucT2_TOP2A_SmallCTR <- c("TOP2A")
cds_sub_DucT2_TOP2A_SmallCTR_TOP2A <- cds_sub_DucT2_TOP2A_SmallCTR[rowData(cds_sub_DucT2_TOP2A_SmallCTR)$gene_short_name %in% ciliated_genes_TOP2A_DucT2_TOP2A_SmallCTR,]
plot_genes_violin(cds_sub_DucT2_TOP2A_SmallCTR_TOP2A, group_cells_by="cell_cycle", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
DucT2_TOP2A_SmallCTR_TOP2A <- as.data.frame(cds_sub_DucT2_TOP2A_SmallCTR_TOP2A@assays@data@listData[["counts"]])
# TOP2A <- as.data.frame(cds_sub_AcinaDucT_NewK_ReCluster_TOP2A@assays@data@listData[["logcounts"]])

# CGAS
ciliated_genes_CGAS_DucT2_TOP2A_SmallCTR <- c("CGAS")
cds_sub_DucT2_TOP2A_SmallCTR_CGAS <- cds_sub_DucT2_TOP2A_SmallCTR[rowData(cds_sub_DucT2_TOP2A_SmallCTR)$gene_short_name %in% ciliated_genes_CGAS_DucT2_TOP2A_SmallCTR,]
plot_genes_violin(cds_sub_DucT2_TOP2A_SmallCTR_CGAS, group_cells_by="cell_cycle", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
DucT2_TOP2A_SmallCTR_CGAS <- as.data.frame(cds_sub_DucT2_TOP2A_SmallCTR_CGAS@assays@data@listData[["counts"]])

##--------------------------------------------------------------------------------------------------------------
## NPC
## NPC vs PTK2
NP_Sum_PTK2_DucT2_TOP2A_SmallCTR <- as.data.frame(cbind(DucT2_TOP2A_SmallCTR_NPC,t(DucT2_TOP2A_SmallCTR_PTK2)))
ggscatter(NP_Sum_PTK2_DucT2_TOP2A_SmallCTR, x = "PTK2", y = "DucT2_TOP2A_SmallCTR_NPC", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",  #"spearman"
          xlab = "PTK2", ylab = "NPC")

## NPC vs TOP2A
NP_Sum_TOP2A_DucT2_TOP2A_SmallCTR <- as.data.frame(cbind(DucT2_TOP2A_SmallCTR_NPC,t(DucT2_TOP2A_SmallCTR_TOP2A)))
ggscatter(NP_Sum_TOP2A_DucT2_TOP2A_SmallCTR, x = "TOP2A", y = "DucT2_TOP2A_SmallCTR_NPC", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "TOP2A", ylab = "NPC")

## NPC vs CGAS
NP_Sum_CGAS_DucT2_TOP2A_SmallCTR  <- as.data.frame(cbind(DucT2_TOP2A_SmallCTR_NPC,t(DucT2_TOP2A_SmallCTR_CGAS)))
ggscatter(NP_Sum_CGAS_DucT2_TOP2A_SmallCTR, x = "CGAS", y = "DucT2_TOP2A_SmallCTR_NPC", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",  #"spearman"
          xlab = "CGAS", ylab = "NPC")

## NE
## NE vs PTK2
NE_Sum_PTK2_DucT2_TOP2A_SmallCTR <- as.data.frame(cbind(DucT2_TOP2A_SmallCTR_NE,t(DucT2_TOP2A_SmallCTR_PTK2)))
ggscatter(NE_Sum_PTK2_DucT2_TOP2A_SmallCTR, x = "PTK2", y = "DucT2_TOP2A_SmallCTR_NE", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",  #"spearman"
          xlab = "PTK2", ylab = "NE")

## NE vs TOP2A
NE_Sum_TOP2A_DucT2_TOP2A_SmallCTR <- as.data.frame(cbind(DucT2_TOP2A_SmallCTR_NE,t(DucT2_TOP2A_SmallCTR_TOP2A)))
ggscatter(NE_Sum_TOP2A_DucT2_TOP2A_SmallCTR, x = "TOP2A", y = "DucT2_TOP2A_SmallCTR_NE", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "TOP2A", ylab = "NE")

## NE vs CGAS
NE_Sum_CGAS_DucT2_TOP2A_SmallCTR  <- as.data.frame(cbind(DucT2_TOP2A_SmallCTR_NE,t(DucT2_TOP2A_SmallCTR_CGAS)))
ggscatter(NE_Sum_CGAS_DucT2_TOP2A_SmallCTR, x = "CGAS", y = "DucT2_TOP2A_SmallCTR_NE", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",  #"spearman"
          xlab = "CGAS", ylab = "NE")


## ACST
## ACST vs PTK2
ACST_Sum_PTK2_DucT2_TOP2A_SmallCTR <- as.data.frame(cbind(DucT2_TOP2A_SmallCTR_ACST,t(DucT2_TOP2A_SmallCTR_PTK2)))
ggscatter(ACST_Sum_PTK2_DucT2_TOP2A_SmallCTR, x = "PTK2", y = "DucT2_TOP2A_SmallCTR_ACST", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",  #"spearman"
          xlab = "PTK2", ylab = "ACST")

## ACST vs TOP2A
ACST_Sum_TOP2A_DucT2_TOP2A_SmallCTR <- as.data.frame(cbind(DucT2_TOP2A_SmallCTR_ACST,t(DucT2_TOP2A_SmallCTR_TOP2A)))
ggscatter(ACST_Sum_TOP2A_DucT2_TOP2A_SmallCTR, x = "TOP2A", y = "DucT2_TOP2A_SmallCTR_ACST", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "TOP2A", ylab = "ACST")

## ACST vs CGAS
ACST_Sum_CGAS_DucT2_TOP2A_SmallCTR  <- as.data.frame(cbind(DucT2_TOP2A_SmallCTR_ACST,t(DucT2_TOP2A_SmallCTR_CGAS)))
ggscatter(ACST_Sum_CGAS_DucT2_TOP2A_SmallCTR, x = "CGAS", y = "DucT2_TOP2A_SmallCTR_ACST", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",  #"spearman"
          xlab = "CGAS", ylab = "ACST")

