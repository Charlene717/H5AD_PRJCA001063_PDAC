rm(list=(ls()[ls()!="cds_sub_AcinaDucT_NewK_ReCluster"]))

#### AcinaDucT
## TOP2A
ciliated_genes_TOP2A <- c("TOP2A")
cds_sub_AcinaDucT_NewK_ReCluster_TOP2A <- cds_sub_AcinaDucT_NewK_ReCluster[rowData(cds_sub_AcinaDucT_NewK_ReCluster)$gene_short_name %in% ciliated_genes_TOP2A,]
plot_genes_violin(cds_sub_AcinaDucT_NewK_ReCluster_TOP2A, group_cells_by="cell_cycle", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
TOP2A <- as.data.frame(cds_sub_AcinaDucT_NewK_ReCluster_TOP2A@assays@data@listData[["counts"]])


# cds_sub_AcinaDucT_NewK_ReCluster_TOP2A <- cds_sub_AcinaDucT_NewK_ReCluster_TOP2A[,cds_sub_AcinaDucT_NewK_ReCluster_TOP2A@colData@listData[["ReCluster"]] %in% c("AC")]

## TOP2A G1
cds_sub_AcinaDucT_NewK_ReCluster_TOP2A_G1 <- cds_sub_AcinaDucT_NewK_ReCluster_TOP2A[,cds_sub_AcinaDucT_NewK_ReCluster_TOP2A@colData@listData[["cell_cycle"]] %in% c("G1")]
# plot_genes_violin(cds_sub_AcinaDucT_NewK_ReCluster_coreCD00_TOP2A_G1, group_cells_by="cell_cycle", ncol=2) +
#   theme(axis.text.x=element_text(angle=45, hjust=1))
TOP2A_G1 <- as.data.frame(cds_sub_AcinaDucT_NewK_ReCluster_TOP2A_G1@assays@data@listData[["counts"]])

## TOP2A G2M
cds_sub_AcinaDucT_NewK_ReCluster_TOP2A_G2M <- cds_sub_AcinaDucT_NewK_ReCluster_TOP2A[,cds_sub_AcinaDucT_NewK_ReCluster_TOP2A@colData@listData[["cell_cycle"]] %in% c("G2M")]
# plot_genes_violin(cds_sub_AcinaDucT_NewK_ReCluster_coreCD00_TOP2A_G2M, group_cells_by="cell_cycle", ncol=2) +
#   theme(axis.text.x=element_text(angle=45, hjust=1))
TOP2A_G2M <- as.data.frame(cds_sub_AcinaDucT_NewK_ReCluster_TOP2A_G2M@assays@data@listData[["counts"]])

## TOP2A S
cds_sub_AcinaDucT_NewK_ReCluster_TOP2A_S <- cds_sub_AcinaDucT_NewK_ReCluster_TOP2A[,cds_sub_AcinaDucT_NewK_ReCluster_TOP2A@colData@listData[["cell_cycle"]] %in% c("S")]
# plot_genes_violin(cds_sub_AcinaDucT_NewK_ReCluster_coreCD00_TOP2A_S, group_cells_by="cell_cycle", ncol=2) +
#   theme(axis.text.x=element_text(angle=45, hjust=1))
TOP2A_S <- as.data.frame(cds_sub_AcinaDucT_NewK_ReCluster_TOP2A_S@assays@data@listData[["counts"]])

##-------------------------- coreCD00 --------------------------##
cds_sub_AcinaDucT_NewK_ReCluster_coreCD00 <- cds_sub_AcinaDucT_NewK_ReCluster[,cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["ReCluster"]] %in% c("CoreCD00")]
## TOP2A
ciliated_genes_TOP2A <- c("TOP2A")
cds_sub_AcinaDucT_NewK_ReCluster_coreCD00_TOP2A <- cds_sub_AcinaDucT_NewK_ReCluster_coreCD00[rowData(cds_sub_AcinaDucT_NewK_ReCluster_coreCD00)$gene_short_name %in% ciliated_genes_TOP2A,]
plot_genes_violin(cds_sub_AcinaDucT_NewK_ReCluster_coreCD00_TOP2A, group_cells_by="cell_cycle", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
TOP2A_coreCD00 <- as.data.frame(cds_sub_AcinaDucT_NewK_ReCluster_coreCD00_TOP2A@assays@data@listData[["counts"]])

## TOP2A G1
cds_sub_AcinaDucT_NewK_ReCluster_coreCD00_TOP2A_G1 <- cds_sub_AcinaDucT_NewK_ReCluster_coreCD00_TOP2A[,cds_sub_AcinaDucT_NewK_ReCluster_coreCD00_TOP2A@colData@listData[["cell_cycle"]] %in% c("G1")]
# plot_genes_violin(cds_sub_AcinaDucT_NewK_ReCluster_coreCD00_TOP2A_G1, group_cells_by="cell_cycle", ncol=2) +
#   theme(axis.text.x=element_text(angle=45, hjust=1))
TOP2A_coreCD00_G1 <- as.data.frame(cds_sub_AcinaDucT_NewK_ReCluster_coreCD00_TOP2A_G1@assays@data@listData[["counts"]])


TOP2A_G1_2 <- t(TOP2A_G1)
TOP2A_G1_2_A <- rep("AcinaDucT", times=length(TOP2A_G1_2))
TOP2A_G1_3 <- as.data.frame(cbind(TOP2A_G1_2,TOP2A_G1_2_A))
colnames(TOP2A_G1_3)[2] <- c("Type")
  
  
TOP2A_coreCD00_G1_2 <- t(TOP2A_coreCD00_G1)
TOP2A_coreCD00_G1_2_A <- rep("CoreCD00", times=length(TOP2A_coreCD00_G1_2))
TOP2A_coreCD00_G1_3 <- as.data.frame(cbind(TOP2A_coreCD00_G1_2,TOP2A_coreCD00_G1_2_A))
colnames(TOP2A_coreCD00_G1_3)[2] <- c("Type")

# TOP2A_coreCD00_G1_3 <-TOP2A_coreCD00_G1_3+log10(TOP2A_coreCD00_G1_3) 


TOP2A_G1_Sum <- rbind(TOP2A_G1_3,TOP2A_coreCD00_G1_3 )
#!!! We need to trans the factor to numeric in here
TOP2A_G1_Sum$TOP2A <- as.numeric(TOP2A_G1_Sum$TOP2A)
# ####
# library(reshape2)
# TOP2A_G1_melt <- melt(TOP2A_G12,TOP2A_coreCD00_G12)


##-------------------------- Plot --------------------------##
ggplot_TOP2A_G1_Sum <- ggplot(data =TOP2A_G1_Sum, aes(x = Type , y = TOP2A , fill =Type ))
ggplot_TOP2A_G1_Sum 
ggplot_TOP2A_G1_Sum  + geom_violin()
ggplot_TOP2A_G1_Sum  + geom_boxplot(alpha = 0.5, show.legend = FALSE)
ggplot_TOP2A_G1_Sum + geom_violin(trim = FALSE) +
  stat_summary(fun= mean, geom = "point",
               shape = 23, size = 2, color = "blue")
library(ggpubr)
# https://www.r-bloggers.com/2017/06/add-p-values-and-significance-levels-to-ggplots/
ggplot_TOP2A_G1_Sum + geom_violin(trim = FALSE, size = 0.8) +
 # scale_y_continuous(trans = 'log10', limits = c(1, 5))+
  
   theme(panel.background = element_rect(fill = "white", colour = "black",
                                    size = 2, linetype = "solid"))+
  #theme_bw() + #  background
  stat_summary(fun= mean, geom = "point",
               shape = 18, size = 3, color = "black")+
  theme(axis.text.x = element_text(face="bold", # color="#993333", 
                                   size=14,color="black", angle=0,vjust=0.5),
        axis.text.y = element_text(face="bold",color="black",  size=14),
        axis.title.x = element_text(size = 14,face="bold"),
        axis.title.y = element_text(size = 14,face="bold"))+
  geom_hline(yintercept = mean(TOP2A_G1_Sum$TOP2A), linetype = 2)+ # Add horizontal line at base mean
  scale_fill_manual(values=c("#edd15f", "#db993b"))+
  stat_compare_means( aes(label = ..p.signif..), 
                      label.x = 1.5, label.y = 4,size = 7)+
#  theme(axis.line.x = element_line(colour = "black", size = 0.5),
#        axis.line.y = element_line(colour = "black", size = 0.5))+
  theme(legend.position="top",legend.text=element_text(size = 14),legend.title=element_text(size = 14,face="bold"))
  
# ggplot_TOP2A_G1_Sum + geom_boxplot(trim = FALSE, size = 0.8) + 
#   theme(panel.background = element_rect(fill = "white", colour = "black",
#                                         size = 2, linetype = "solid"))+
#   #theme_bw() + #  background
#   stat_summary(fun= mean, geom = "point",
#                shape = 18, size = 3, color = "black")+
#   theme(axis.text.x = element_text(face="bold", # color="#993333", 
#                                    size=14,color="black", angle=0,vjust=0.5),
#         axis.text.y = element_text(face="bold",color="black",  size=14),
#         axis.title.x = element_text(size = 14,face="bold"),
#         axis.title.y = element_text(size = 14,face="bold"))+
#   geom_hline(yintercept = mean(TOP2A_G1_Sum$TOP2A), linetype = 2)+ # Add horizontal line at base mean
#   scale_fill_manual(values=c("#edd15f", "#db993b"))+
#   stat_compare_means( aes(label = ..p.signif..), 
#                       label.x = 1.5, label.y = 5,size = 7)+
#   #  theme(axis.line.x = element_line(colour = "black", size = 0.5),
#   #        axis.line.y = element_line(colour = "black", size = 0.5))+
#   theme(legend.position="top",legend.text=element_text(size = 14),legend.title=element_text(size = 14,face="bold"))

