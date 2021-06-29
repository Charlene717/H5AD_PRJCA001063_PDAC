rm(list=(ls()[ls()!="cds_sub_AcinaDucT_NewK_ReCluster"]))

GeneExpMatrix <- cds_sub_AcinaDucT_NewK_ReCluster@assays@data@listData[["counts"]]
GeneExpMatrix_DF <- as.data.frame(GeneExpMatrix)
GeneExpMatrix_DF2 <- t(GeneExpMatrix_DF)
# GeneExpMatrix_DF3 <- cbind(row.names(GeneExpMatrix_DF2),GeneExpMatrix_DF2)
# colnames(GeneExpMatrix_DF3)[1] <- c("Barcode")

# PhenoType
Cell <- as.character(cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["CELL"]])
Patient <- as.character(cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["Patient"]])
Type <- as.character(cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["Type"]])
Cell_type <- as.character(cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["Cell_type"]])
cell_cycle <- as.character(cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["cell_cycle"]])
ReCluster <- as.character(cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["ReCluster"]])

# Marker
PDAC <- cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["PDAC"]]
EMT <- cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["EMT"]]
Migration <- cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["Migration"]]
Metastasis <- cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["Metastasis"]]
ReCluster <- cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["ReCluster"]]
NE <- cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["NE"]]
NPC <- cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["NPC"]]
ACST <- cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["ACST"]]

# Combine
GeneExpMatrix_DF4 <- cbind(Cell,Patient,Type,Cell_type,ReCluster,PDAC,EMT,Migration,Metastasis,ReCluster,NE,NPC,ACST,cell_cycle,GeneExpMatrix_DF2)
GeneExpMatrix_DF4 <- cbind(row.names(GeneExpMatrix_DF4),GeneExpMatrix_DF4)
colnames(GeneExpMatrix_DF4)[1] <- c("Barcode")

PathName = setwd(getwd())
RVersion = "20210626V1_GEM"
dir.create(paste0(PathName,"/",RVersion))

# write.table(GeneExpMatrix_DF4, file = paste0(PathName,"/AcinaDucT_GeneExpMatrix_Pheno.txt"),sep = " ",header= T, quote = FALSE, na = "NA")
write.table(GeneExpMatrix_DF4, file=paste0(PathName,"/",RVersion,"/AcinaDucT_GeneExpMatrix_Pheno_Check_Marker.txt"),sep="\t", row.names=F)

# GeneExpMatrix_DF2 <- cbind(row.names(GeneExpMatrix_DF),GeneExpMatrix_DF)
# colnames(GeneExpMatrix_DF2)[1] <- c("gene_short_name")
