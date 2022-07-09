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

##----------------------------Heatmap & bar-------------------------------##
## Bar
colors_cc <-c("#FF9912B3", "#2e6087", "#417034")  ## Color for Cell-Cycle

annotation_col <- as.data.frame(cbind(cell_cycle,ReCluster))
library(ggplot2)
ggplot(annotation_col, aes(x = ReCluster, fill = cell_cycle)) + 
  geom_bar(position = "fill")+
  theme(axis.text.x = element_text(face="bold", # color="#993333", 
                                   size=10, angle=75,vjust =0.55),
        axis.text.y = element_text(face="bold",hjust =-8),
        axis.title.x = element_text(size = 14,face="bold"),
        axis.title.y = element_text(size = 14,face="bold"),
        plot.title = element_text(color="black", size=20, 
                                  face="bold.italic",hjust = 0.1,vjust =-8), # margin = margin(t = 0.5, b = -7),
        legend.title = element_text(size=12, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"))+ scale_fill_manual(values = colors_cc)

# ## Heatmap & bar
# library(pheatmap)
# ##pheatmap(mat, annotation_col = annotation_col, annotation_row = annotation_row)
# pheatmap(GeneExpMatrix_DF2, annotation_col = annotation_col)
##----------------------------Heatmap & bar-------------------------------##

# Marker
PDAC <- cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["PDAC"]]
EMT <- cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["EMT"]]
Migration <- cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["Migration"]]
Metastasis <- cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["Metastasis"]]
NE <- cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["NE"]]
NPC <- cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["NPC"]]
ACST <- cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["ACST"]]
HYPOXIA <- cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["HYPOXIA"]]

# Combine
GeneExpMatrix_DF4 <- cbind(Cell,Patient,Type,Cell_type,ReCluster,PDAC,EMT,Migration,Metastasis,NE,NPC,ACST,HYPOXIA,cell_cycle,GeneExpMatrix_DF2)
# GeneExpMatrix_DF4 <- cbind(row.names(GeneExpMatrix_DF4),GeneExpMatrix_DF4)
GeneExpMatrix_DF4 <- data.frame(row.names(GeneExpMatrix_DF4),GeneExpMatrix_DF4)
colnames(GeneExpMatrix_DF4)[1] <- c("Barcode")

PathName = setwd(getwd())
RVersion = "20210815V1_GEM"
dir.create(paste0(PathName,"/",RVersion))

# write.table(GeneExpMatrix_DF4, file = paste0(PathName,"/AcinaDucT_GeneExpMatrix_Pheno.txt"),sep = " ",header= T, quote = FALSE, na = "NA")
write.table(GeneExpMatrix_DF4, file=paste0(PathName,"/",RVersion,"/AcinaDucT_GeneExpMatrix_Pheno_Check_Marker.txt"),sep="\t", row.names=F)

# GeneExpMatrix_DF2 <- cbind(row.names(GeneExpMatrix_DF),GeneExpMatrix_DF)
# colnames(GeneExpMatrix_DF2)[1] <- c("gene_short_name")


#-------------------------------------------------------------------------#

plot_genes_by_group(cds_sub_AcinaDucT_NewK_ReCluster,
                    c("TOP2A","PTK2","CGAS","TP53","COL17A1"),
                    group_cells_by="ReCluster",
                    ordering_type="cluster_row_col",
                    max.size=5)

library(ggplot2)
plot_genes_by_group(cds_sub_AcinaDucT_NewK_ReCluster,
                    c("TOP2A","PTK2","CGAS","TP53","COL17A1"),
                    group_cells_by="ReCluster",
                    ordering_type="maximal_on_diag",
                    max.size=5)+
  theme(axis.text.x=element_text(angle=80, hjust=1))

plot_genes_by_group(cds_sub_AcinaDucT_NewK_ReCluster,
                    c("KRAS","EGFR","BRAF","PIK3CA","MYC"),
                    group_cells_by="ReCluster",
                    ordering_type="maximal_on_diag",
                    max.size=5)+
theme(axis.text.x=element_text(angle=80, hjust=1))

# nuclear FAK Marker
plot_genes_by_group(cds_sub_AcinaDucT_NewK_ReCluster,
                    c("H19","NCAM1","MYOG","GATA4","VCAM1"),
                    group_cells_by="ReCluster",
                    ordering_type="maximal_on_diag",
                    max.size=5)+
  theme(axis.text.x=element_text(angle=80, hjust=1))

# nuclear FAK Marker 2
plot_genes_by_group(cds_sub_AcinaDucT_NewK_ReCluster,
                    # c("Has2","Has3","Loxl4","Lox","Col10a1","Col13a1","Col14a1","Col17a1","Col8a1",
                    #   "Fgf10","Fgf11","Fgf18"),
                    c("HAS2","HAS3","LOXL4","LOX","COL8A1","COL10A1","COL13A1","COL14A1","COL17A1",
                      "FGF10","FGF11","FGF18"),
                    group_cells_by="ReCluster",
                    ordering_type="maximal_on_diag",
                    max.size=5)+
  theme(axis.text.x=element_text(angle=80, hjust=1))

# nuclear FAK Marker 2
plot_genes_by_group(cds_sub_AcinaDucT_NewK_ReCluster,
                    # c("Flt3lg","M-CSF","GM-CSF","CX3CL1","CCL2","CCL11","CXCL2","G-CSF","DKK-1",
                    #   "IL-1α","IL-1β","E-Selectin","MMP-3","EPC-1","Gas6","CCL6","CXCL13","CCL20"),
                    c("FLT3LG","CSF1","CSF2","CX3CL1","CCL2","CCL11","CXCL2","CSF3","DKK1",
                      "IL1A","IL1B","SELE","MMP3","SERPINF1","GAS6","CCL6","CXCL13","CCL20"),
                    group_cells_by="ReCluster",
                    ordering_type="maximal_on_diag",
                    max.size=5)+
  theme(axis.text.x=element_text(angle=80, hjust=1))

# nuclear FAK Marker 4
plot_genes_by_group(cds_sub_AcinaDucT_NewK_ReCluster,
                    c("XIAP","PIAS1","TLN1","TLN2","ARHGEF28","PXN"),
                    group_cells_by="ReCluster",
                    ordering_type="maximal_on_diag",
                    max.size=5)+
  theme(axis.text.x=element_text(angle=80, hjust=1))


## Error
plot_genes_by_group(cds_sub_AcinaDucT_NewK_ReCluster,
                    c("PDAC","EMT","TP53"),
                    group_cells_by="ReCluster",
                    ordering_type="maximal_on_diag",
                    max.size=5)+
  theme(axis.text.x=element_text(angle=80, hjust=1))
#-------------------------------------------------------------------------#
GeneExpMatrix_DF4_Customize <- as.data.frame(GeneExpMatrix_DF4)

GeneExpMatrix_DF4_Customize2 <- GeneExpMatrix_DF4_Customize[ ,colnames(GeneExpMatrix_DF4_Customize)%in% 
                                                         c("ReCluster","TOP2A","PTK2","CGAS","TP53","GATA4","PDAC","EMT","Migration",
                                                           "Metastasis","NE","NPC","ACST","HYPOXIA")]
head(GeneExpMatrix_DF4_Customize2)

GeneExpMatrix_DF4_Customize_GEM <- GeneExpMatrix_DF4_Customize2[,2:14]
GeneExpMatrix_DF4_Customize_GEM <- as.data.frame(t(GeneExpMatrix_DF4_Customize_GEM))
GeneExpMatrix_DF4_Customize_GEM2 <- t( GeneExpMatrix_DF4_Customize2[,2:14])
GeneExpMatrix_DF4_Customize_Pheno <- as.data.frame(GeneExpMatrix_DF4_Customize2[,1])
GeneExpMatrix_DF4_Customize_Pheno <- cbind(row.names(GeneExpMatrix_DF4_Customize_Pheno),GeneExpMatrix_DF4_Customize_Pheno)
colnames(GeneExpMatrix_DF4_Customize_Pheno) <- c("samples","ReCluster")
GeneExpMatrix_DF4_Customize_GL <- as.data.frame(colnames(GeneExpMatrix_DF4_Customize2[,2:14]))
colnames(GeneExpMatrix_DF4_Customize_GL) <- c("gene_short_name")
row.names(GeneExpMatrix_DF4_Customize_GL) <- GeneExpMatrix_DF4_Customize_GL[,1]

cds_Customize <- new_cell_data_set(GeneExpMatrix_DF4_Customize_GEM2,
                         cell_metadata = GeneExpMatrix_DF4_Customize_Pheno,
                         gene_metadata = GeneExpMatrix_DF4_Customize_GL)

plot_genes_by_group(cds_Customize,
                    c("TOP2A","PTK2","CGAS","TP53","PDAC","GATA4","EMT","Migration",
                      "Metastasis","NE","NPC","ACST","HYPOXIA"),
                    group_cells_by="ReCluster",
                    ordering_type="maximal_on_diag",
                    max.size=5)+
  theme(axis.text.x=element_text(angle=80, hjust=1))

# Marker
plot_genes_by_group(cds_Customize,
                    c("PDAC","EMT","Migration",
                      "Metastasis","NE","NPC","ACST","HYPOXIA"),
                    group_cells_by="ReCluster",
                    ordering_type="maximal_on_diag",
                    max.size=10)+
#  scale_colour_gradient2(low = "#440075", mid = "#ffd261", high = "#4aff8c")+
  scale_colour_gradient2(low = "#440075", mid = "#e88a4f", high = "#e0ff63")+
  theme(axis.text.x=element_text(angle=80, hjust=1))+
  theme(axis.text.x=element_text(angle=80, hjust=1))+
  theme(axis.line.x = element_line(colour = "black", size = 0.7),
        axis.line.y = element_line(colour = "black", size = 0.7))+
  theme(axis.text.x = element_text(face="bold",color="black",  size=13),
        axis.text.y = element_text(face="bold",color="black",  size=13),
        axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
        axis.title = element_text(size = rel(1.5),face="bold"),
        plot.title = element_text(color="black", size=20, 
                                  face="bold.italic",hjust = 0.1,vjust =-8), # margin = margin(t = 0.5, b = -7),
        legend.title = element_text(size=12, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"))

# Genes #5:20 for PDF
plot_genes_by_group(cds_Customize,
                    c("TOP2A","PTK2","CGAS","TP53","GATA4"),
                    group_cells_by="ReCluster",
                    ordering_type="maximal_on_diag",
                    max.size=10)+
  scale_colour_gradient2(low = "black", mid = "#3528c7", high = "#ff85e9", 
                         guide = "colourbar",midpoint = 0, labs(fill ="Exp"))+
  theme(axis.text.x=element_text(angle=80, hjust=1))+
  theme(axis.line.x = element_line(colour = "black", size = 0.7),
        axis.line.y = element_line(colour = "black", size = 0.7))+
  theme(axis.text.x = element_text(face="bold",color="black",  size=13),
        axis.text.y = element_text(face="bold",color="black",  size=13),
        axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
        axis.title = element_text(size = rel(1.5),face="bold"),
        plot.title = element_text(color="black", size=20, 
                                  face="bold.italic",hjust = 0.1,vjust =-8), # margin = margin(t = 0.5, b = -7),
        legend.title = element_text(size=12, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"))
        
       #  + 
       #    scale_color_hue(l=40, c=35)+ 
       #    scale_fill_brewer(palette = "Spectral")+
       #    scale_fill_gradient(low = "blue", high = "red")+
       #    
       #    scale_colour_gradient2(low = "black", mid = "#3528c7", high = "#ff85e9", 
       #                           guide = "colourbar",midpoint = 0, labs(fill = "Exp"))
       #  #legend.position = c(0.1, 0.18),
       # # aspect.ratio=1
       #  ) 
#-------------------------------------------------------------------------#
# https://statisticsglobe.com/common-main-title-for-multiple-plots-in-r
# https://stackoverflow.com/questions/36458594/plotting-figures-using-parmfrow-c-in-r
ggp1 <- plot_genes_by_group(cds_Customize,
                    c("PDAC","EMT","Migration",
                      "Metastasis","NE","NPC","ACST"),
                    group_cells_by="ReCluster",
                    ordering_type="maximal_on_diag",
                    max.size=7)+
  #  scale_colour_gradient2(low = "#440075", mid = "#ffd261", high = "#4aff8c")+
  scale_colour_gradient2(low = "#440075", mid = "#e88a4f", high = "#e0ff63")+
  theme(axis.text.x=element_text(angle=80, hjust=1))+
  theme(axis.text.x=element_text(angle=80, hjust=1))+
  theme(axis.line.x = element_line(colour = "black", size = 0.7),
        axis.line.y = element_line(colour = "black", size = 0.7))+
  theme(axis.text.x = element_text(face="bold",color="black",  size=13),
        axis.text.y = element_text(face="bold",color="black",  size=13),
        axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
        axis.title = element_text(size = rel(1.5),face="bold"),
        plot.title = element_text(color="black", size=20, 
                                  face="bold.italic",hjust = 0.1,vjust =-8), # margin = margin(t = 0.5, b = -7),
        legend.title = element_text(size=12, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"))
ggp2 <- plot_genes_by_group(cds_Customize,
                    c("TOP2A","PTK2","CGAS","TP53"),
                    group_cells_by="ReCluster",
                    ordering_type="maximal_on_diag",
                    max.size=7)+
  scale_colour_gradient2(low = "black", mid = "#3528c7", high = "#ff85e9", 
                         guide = "colourbar",midpoint = 0, labs(fill ="Exp"))+
  theme(axis.text.x=element_text(angle=80, hjust=1))+
  theme(axis.line.x = element_line(colour = "black", size = 0.7),
        axis.line.y = element_line(colour = "black", size = 0.7))+
  theme(axis.text.x = element_text(face="bold",color="black",  size=13),
        axis.text.y = element_text(face="bold",color="black",  size=13),
        axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
        axis.title = element_text(size = rel(1.5),face="bold"),
        plot.title = element_text(color="black", size=20, 
                                  face="bold.italic",hjust = 0.1,vjust =-8), # margin = margin(t = 0.5, b = -7),
        legend.title = element_text(size=12, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"))


ggp3 <-  plot_genes_by_group(cds_sub_AcinaDucT_NewK_ReCluster,
                            c("KRAS","EGFR","BRAF","PIK3CA","MYC"),
                            group_cells_by="ReCluster",
                            ordering_type="maximal_on_diag",
                            max.size=7)+
  # scale_colour_gradient2(low = "black", mid = "#3528c7", high = "#ff85e9", 
  #                        guide = "colourbar",midpoint = 0, labs(fill ="Exp"))+
        theme(axis.text.x=element_text(angle=80, hjust=1))+
        theme(axis.line.x = element_line(colour = "black", size = 0.7),
              axis.line.y = element_line(colour = "black", size = 0.7))+
          theme(axis.text.x = element_text(face="bold",color="black",  size=13),
                axis.text.y = element_text(face="bold",color="black",  size=13),
                axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
                axis.title = element_text(size = rel(1.5),face="bold"),
                plot.title = element_text(color="black", size=20, 
                                          face="bold.italic",hjust = 0.1,vjust =-8), # margin = margin(t = 0.5, b = -7),
                legend.title = element_text(size=12, color = "black", face="bold"),
                legend.text = element_text(colour="black", size=12,face="bold"))
ggp3

library("patchwork")        
gp_all <- (ggp1) / (ggp2) / (ggp3)   # Create grid of plots with title
gp_all                                       # Draw grid of plots with title

# gp_all <- (ggp1 + ggp2) / (ggp3 + ggp4) +    # Create grid of plots with title
#   plot_annotation(title = "My Multiplot Title") & 
#   theme(plot.title = element_text(hjust = 0.5))
# ggp_all                                       # Draw grid of plots with title

#-------------------------------------------------------------------------#

# GeneExpMatrix_DF4_Bubble <- as.data.frame(GeneExpMatrix_DF4)
# 
# GeneExpMatrix_DF4_Bubble2 <- GeneExpMatrix_DF4_Bubble[ ,colnames(GeneExpMatrix_DF4_Bubble)%in% 
#                                                         c("ReCluster","TOP2A","PTK2","CGAS","TP53","PDAC","EMT","Migration",
#                                                           "Metastasis","NE","NPC","ACST")]
# head(GeneExpMatrix_DF4_Bubble2)
# #-------------------------------------------------------------------------#
# 
# 
# #-------------------------------------------------------------------------#
# 
# library(ggBubbles)
# library(tidyverse)
# library(ggplot2)
# 
# library(ggpubr)
# GeneExpMatrix_DF4_Bubble3 <- GeneExpMatrix_DF4_Bubble2
# row.names(GeneExpMatrix_DF4_Bubble3) <- GeneExpMatrix_DF4_Bubble3[,1]
# ggballoonplot(GeneExpMatrix_DF4_Bubble2, fill = "value")+
#   scale_fill_viridis_c(option = "C")
# 
# library(reshape2)
# data_melt <- melt(GeneExpMatrix_DF4_Bubble2)
# data_melt2 <- melt(GeneExpMatrix_DF4_Bubble2[,3:4])
# names(data_melt) = c('Gene', 'Cell', 'Value')
# p<-ggplot(GeneExpMatrix_DF4_Bubble2, aes(x = TOP2A, y = ReCluster, size = Value, color=Cell)) + geom_point()
# p<-ggplot(GeneExpMatrix_DF4_Bubble2, aes(x = TOP2A, y = ReCluster)) + geom_point()
# p
# 
# df <- GeneExpMatrix_DF4_Bubble %>%
#       group_by(ReCluster)%>%
#       summarize(Count = n(), AvgLe)