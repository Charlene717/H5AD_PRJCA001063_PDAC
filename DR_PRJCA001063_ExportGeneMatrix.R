rm(list=(ls()[ls()!="cds_sub_AcinaDucT_NewK_ReCluster"]))

GeneExpMatrix <- cds_sub_AcinaDucT_NewK_ReCluster@assays@data@listData[["counts"]]
GeneExpMatrix_DF <- as.data.frame(GeneExpMatrix)
GeneExpMatrix_DF2 <- t(GeneExpMatrix_DF)
# GeneExpMatrix_DF3 <- cbind(row.names(GeneExpMatrix_DF2),GeneExpMatrix_DF2)
# colnames(GeneExpMatrix_DF3)[1] <- c("Barcode")

Patient <- as.character(cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["Patient"]])
Type <- as.character(cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["Type"]])
Cell_type <- as.character(cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["Cell_type"]])
cell_cycle <- as.character(cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["cell_cycle"]])
ReCluster <- as.character(cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["ReCluster"]])

GeneExpMatrix_DF4 <- cbind(Patient,Type,Cell_type,cell_cycle,ReCluster,GeneExpMatrix_DF2)
GeneExpMatrix_DF4 <- cbind(row.names(GeneExpMatrix_DF4),GeneExpMatrix_DF4)
colnames(GeneExpMatrix_DF4)[1] <- c("Barcode")

PathName = setwd(getwd())
RVersion = "20210625V1_GEM"
dir.create(paste0(PathName,"/",RVersion))

# write.table(GeneExpMatrix_DF4, file = paste0(PathName,"/AcinaDucT_GeneExpMatrix_Pheno.txt"),sep = " ",header= T, quote = FALSE, na = "NA")
write.table(GeneExpMatrix_DF4, file=paste0(PathName,"/",RVersion,"/AcinaDucT_GeneExpMatrix_Pheno.txt"),sep="\t", row.names=F)

# GeneExpMatrix_DF2 <- cbind(row.names(GeneExpMatrix_DF),GeneExpMatrix_DF)
# colnames(GeneExpMatrix_DF2)[1] <- c("gene_short_name")
