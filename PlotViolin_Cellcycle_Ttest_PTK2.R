rm(list=(ls()[ls()!="cds_sub_AcinaDucT_NewK_ReCluster"]))

#### AcinaDucT
## PTK2
ciliated_genes_PTK2 <- c("PTK2")
cds_sub_AcinaDucT_NewK_ReCluster_PTK2 <- cds_sub_AcinaDucT_NewK_ReCluster[rowData(cds_sub_AcinaDucT_NewK_ReCluster)$gene_short_name %in% ciliated_genes_PTK2,]
plot_genes_violin(cds_sub_AcinaDucT_NewK_ReCluster_PTK2, group_cells_by="cell_cycle", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
PTK2 <- as.data.frame(cds_sub_AcinaDucT_NewK_ReCluster_PTK2@assays@data@listData[["counts"]])


# cds_sub_AcinaDucT_NewK_ReCluster_PTK2 <- cds_sub_AcinaDucT_NewK_ReCluster_PTK2[,cds_sub_AcinaDucT_NewK_ReCluster_PTK2@colData@listData[["ReCluster"]] %in% c("AC")]

####### ------------------------------------- PTK2 G1 ------------------------------------- ####### 

## PTK2 G1
cds_sub_AcinaDucT_NewK_ReCluster_PTK2_G1 <- cds_sub_AcinaDucT_NewK_ReCluster_PTK2[,cds_sub_AcinaDucT_NewK_ReCluster_PTK2@colData@listData[["cell_cycle"]] %in% c("G1")]
# plot_genes_violin(cds_sub_AcinaDucT_NewK_ReCluster_coreCD00_PTK2_G1, group_cells_by="cell_cycle", ncol=2) +
#   theme(axis.text.x=element_text(angle=45, hjust=1))
PTK2_G1 <- as.data.frame(cds_sub_AcinaDucT_NewK_ReCluster_PTK2_G1@assays@data@listData[["counts"]])

## PTK2 G2M
cds_sub_AcinaDucT_NewK_ReCluster_PTK2_G2M <- cds_sub_AcinaDucT_NewK_ReCluster_PTK2[,cds_sub_AcinaDucT_NewK_ReCluster_PTK2@colData@listData[["cell_cycle"]] %in% c("G2M")]
# plot_genes_violin(cds_sub_AcinaDucT_NewK_ReCluster_coreCD00_PTK2_G2M, group_cells_by="cell_cycle", ncol=2) +
#   theme(axis.text.x=element_text(angle=45, hjust=1))
PTK2_G2M <- as.data.frame(cds_sub_AcinaDucT_NewK_ReCluster_PTK2_G2M@assays@data@listData[["counts"]])

## PTK2 S
cds_sub_AcinaDucT_NewK_ReCluster_PTK2_S <- cds_sub_AcinaDucT_NewK_ReCluster_PTK2[,cds_sub_AcinaDucT_NewK_ReCluster_PTK2@colData@listData[["cell_cycle"]] %in% c("S")]
# plot_genes_violin(cds_sub_AcinaDucT_NewK_ReCluster_coreCD00_PTK2_S, group_cells_by="cell_cycle", ncol=2) +
#   theme(axis.text.x=element_text(angle=45, hjust=1))
PTK2_S <- as.data.frame(cds_sub_AcinaDucT_NewK_ReCluster_PTK2_S@assays@data@listData[["counts"]])

##-------------------------- coreCD00 --------------------------##
cds_sub_AcinaDucT_NewK_ReCluster_coreCD00 <- cds_sub_AcinaDucT_NewK_ReCluster[,cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["ReCluster"]] %in% c("CoreCD00")]
## PTK2
ciliated_genes_PTK2 <- c("PTK2")
cds_sub_AcinaDucT_NewK_ReCluster_coreCD00_PTK2 <- cds_sub_AcinaDucT_NewK_ReCluster_coreCD00[rowData(cds_sub_AcinaDucT_NewK_ReCluster_coreCD00)$gene_short_name %in% ciliated_genes_PTK2,]
plot_genes_violin(cds_sub_AcinaDucT_NewK_ReCluster_coreCD00_PTK2, group_cells_by="cell_cycle", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
PTK2_coreCD00 <- as.data.frame(cds_sub_AcinaDucT_NewK_ReCluster_coreCD00_PTK2@assays@data@listData[["counts"]])

## PTK2 G1
cds_sub_AcinaDucT_NewK_ReCluster_coreCD00_PTK2_G1 <- cds_sub_AcinaDucT_NewK_ReCluster_coreCD00_PTK2[,cds_sub_AcinaDucT_NewK_ReCluster_coreCD00_PTK2@colData@listData[["cell_cycle"]] %in% c("G1")]
# plot_genes_violin(cds_sub_AcinaDucT_NewK_ReCluster_coreCD00_PTK2_G1, group_cells_by="cell_cycle", ncol=2) +
#   theme(axis.text.x=element_text(angle=45, hjust=1))
PTK2_coreCD00_G1 <- as.data.frame(cds_sub_AcinaDucT_NewK_ReCluster_coreCD00_PTK2_G1@assays@data@listData[["counts"]])


PTK2_G1_2 <- t(PTK2_G1)
PTK2_G1_2_A <- rep("AcinaDucT", times=length(PTK2_G1_2))
PTK2_G1_3 <- as.data.frame(cbind(PTK2_G1_2,PTK2_G1_2_A))
colnames(PTK2_G1_3)[2] <- c("Type")
  
  
PTK2_coreCD00_G1_2 <- t(PTK2_coreCD00_G1)
PTK2_coreCD00_G1_2_A <- rep("CoreCD00", times=length(PTK2_coreCD00_G1_2))
PTK2_coreCD00_G1_3 <- as.data.frame(cbind(PTK2_coreCD00_G1_2,PTK2_coreCD00_G1_2_A))
colnames(PTK2_coreCD00_G1_3)[2] <- c("Type")

# PTK2_coreCD00_G1_3 <-PTK2_coreCD00_G1_3+log10(PTK2_coreCD00_G1_3) 

###### ----------- Combine Two group ---------------- #####
PTK2_G1_Sum <- rbind(PTK2_G1_3,PTK2_coreCD00_G1_3 )
#!!! We need to trans the factor to numeric in here
PTK2_G1_Sum$PTK2 <- as.numeric(PTK2_G1_Sum$PTK2)
# ####
# library(reshape2)
# PTK2_G1_melt <- melt(PTK2_G12,PTK2_coreCD00_G12)


##-------------------------- Plot --------------------------##
ggplot_PTK2_G1_Sum <- ggplot(data =PTK2_G1_Sum, aes(x = Type , y = PTK2 , fill =Type ))
ggplot_PTK2_G1_Sum 
ggplot_PTK2_G1_Sum  + geom_violin()
ggplot_PTK2_G1_Sum  + geom_boxplot(alpha = 0.5, show.legend = FALSE)
ggplot_PTK2_G1_Sum + geom_violin(trim = FALSE) +
  stat_summary(fun= mean, geom = "point",
               shape = 23, size = 2, color = "blue")
library(ggpubr)
# https://www.r-bloggers.com/2017/06/add-p-values-and-significance-levels-to-ggplots/
ggplot_PTK2_G1_Sum_M <- ggplot_PTK2_G1_Sum + geom_violin(trim = FALSE, size = 0.8) +
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
  geom_hline(yintercept = mean(PTK2_G1_Sum$PTK2), linetype = 2)+ # Add horizontal line at base mean
  scale_fill_manual(values=c("#edd15f", "#db993b"))+
  stat_compare_means( aes(label = ..p.signif..), 
                      label.x = 1.5, label.y = 8,size = 7)+
#  theme(axis.line.x = element_line(colour = "black", size = 0.5),
#        axis.line.y = element_line(colour = "black", size = 0.5))+
  theme(legend.position="top",legend.text=element_text(size = 14),legend.title=element_text(size = 14,face="bold"))+ 
  ggtitle("G1")+ 
  theme(plot.title = element_text(size = 18, vjust = -5 ,hjust = 0,face="bold")
  )

ggplot_PTK2_G1_Sum_M
####### ------------------------------------- PTK2 S ------------------------------------- ####### 
cds_sub_AcinaDucT_NewK_ReCluster_coreCD00_PTK2_S <- cds_sub_AcinaDucT_NewK_ReCluster_coreCD00_PTK2[,cds_sub_AcinaDucT_NewK_ReCluster_coreCD00_PTK2@colData@listData[["cell_cycle"]] %in% c("S")]
PTK2_coreCD00_S <- as.data.frame(cds_sub_AcinaDucT_NewK_ReCluster_coreCD00_PTK2_S@assays@data@listData[["counts"]])

PTK2_S_2 <- t(PTK2_S)
PTK2_S_2_A <- rep("AcinaDucT", times=length(PTK2_S_2))
PTK2_S_3 <- as.data.frame(cbind(PTK2_S_2,PTK2_S_2_A))
colnames(PTK2_S_3)[2] <- c("Type")

PTK2_coreCD00_S_2 <- t(PTK2_coreCD00_S)
PTK2_coreCD00_S_2_A <- rep("CoreCD00", times=length(PTK2_coreCD00_S_2))
PTK2_coreCD00_S_3 <- as.data.frame(cbind(PTK2_coreCD00_S_2,PTK2_coreCD00_S_2_A))
colnames(PTK2_coreCD00_S_3)[2] <- c("Type")

# PTK2_coreCD00_S_3 <-PTK2_coreCD00_S_3+loS0(PTK2_coreCD00_S_3) 

##------ Combine Two group -----------##
PTK2_S_Sum <- rbind(PTK2_S_3,PTK2_coreCD00_S_3 )
#!!! We need to trans the factor to numeric in here
PTK2_S_Sum$PTK2 <- as.numeric(PTK2_S_Sum$PTK2)

##----- Plot -----##
ggplot_PTK2_S_Sum <- ggplot(data =PTK2_S_Sum, aes(x = Type , y = PTK2 , fill =Type ))
ggplot_PTK2_S_Sum 
ggplot_PTK2_S_Sum  + geom_violin()
ggplot_PTK2_S_Sum  + geom_boxplot(alpha = 0.5, show.legend = FALSE)
ggplot_PTK2_S_Sum + geom_violin(trim = FALSE) +
  stat_summary(fun= mean, geom = "point",
               shape = 23, size = 2, color = "blue")
library(ggpubr)
# https://www.r-bloggers.com/2017/06/add-p-values-and-significance-levels-to-ggplots/
ggplot_PTK2_S_Sum_M <- ggplot_PTK2_S_Sum + geom_violin(trim = FALSE, size = 0.8) +
  theme(panel.background = element_rect(fill = "white", colour = "black",
                                        size = 2, linetype = "solid"))+
  stat_summary(fun= mean, geom = "point",
               shape = 18, size = 3, color = "black")+
  theme(axis.text.x = element_text(face="bold", # color="#993333", 
                                   size=14,color="black", angle=0,vjust=0.5),
        axis.text.y = element_text(face="bold",color="black",  size=14),
        axis.title.x = element_text(size = 14,face="bold"),
        axis.title.y = element_text(size = 14,face="bold"))+
  geom_hline(yintercept = mean(PTK2_S_Sum$PTK2), linetype = 2)+ # Add horizontal line at base mean
  scale_fill_manual(values=c("#5fc244", "#417034"))+
  stat_compare_means( aes(label = ..p.signif..), 
                      label.x = 1.5, label.y = 8,size = 7)+
  theme(legend.position="top",legend.text=element_text(size = 14),legend.title=element_text(size = 14,face="bold"))+ 
  ggtitle("S")+ 
  theme(plot.title = element_text(size = 18, vjust = -5 ,hjust = 0,face="bold")
  )

ggplot_PTK2_S_Sum_M

####### ------------------------------------- PTK2 G2M ------------------------------------- #######
cds_sub_AcinaDucT_NewK_ReCluster_coreCD00_PTK2_G2M <- cds_sub_AcinaDucT_NewK_ReCluster_coreCD00_PTK2[,cds_sub_AcinaDucT_NewK_ReCluster_coreCD00_PTK2@colData@listData[["cell_cycle"]] %in% c("G2M")]
PTK2_coreCD00_G2M <- as.data.frame(cds_sub_AcinaDucT_NewK_ReCluster_coreCD00_PTK2_G2M@assays@data@listData[["counts"]])

PTK2_G2M_2 <- t(PTK2_G2M)
PTK2_G2M_2_A <- rep("AcinaDucT", times=length(PTK2_G2M_2))
PTK2_G2M_3 <- as.data.frame(cbind(PTK2_G2M_2,PTK2_G2M_2_A))
colnames(PTK2_G2M_3)[2] <- c("Type")

PTK2_coreCD00_G2M_2 <- t(PTK2_coreCD00_G2M)
PTK2_coreCD00_G2M_2_A <- rep("CoreCD00", times=length(PTK2_coreCD00_G2M_2))
PTK2_coreCD00_G2M_3 <- as.data.frame(cbind(PTK2_coreCD00_G2M_2,PTK2_coreCD00_G2M_2_A))
colnames(PTK2_coreCD00_G2M_3)[2] <- c("Type")

# PTK2_coreCD00_G2M_3 <-PTK2_coreCD00_G2M_3+loG2M0(PTK2_coreCD00_G2M_3) 

##------ Combine Two group -----------##
PTK2_G2M_Sum <- rbind(PTK2_G2M_3,PTK2_coreCD00_G2M_3 )
#!!! We need to trans the factor to numeric in here
PTK2_G2M_Sum$PTK2 <- as.numeric(PTK2_G2M_Sum$PTK2)

##----- Plot -----##
ggplot_PTK2_G2M_Sum <- ggplot(data =PTK2_G2M_Sum, aes(x = Type , y = PTK2 , fill =Type ))
ggplot_PTK2_G2M_Sum 
ggplot_PTK2_G2M_Sum  + geom_violin()
ggplot_PTK2_G2M_Sum  + geom_boxplot(alpha = 0.5, show.legend = FALSE)
ggplot_PTK2_G2M_Sum + geom_violin(trim = FALSE) +
  stat_summary(fun= mean, geom = "point",
               shape = 23, size = 2, color = "blue")
library(ggpubr)
# https://www.r-bloggers.com/2017/06/add-p-values-and-significance-levels-to-ggplots/
ggplot_PTK2_G2M_Sum_M <- ggplot_PTK2_G2M_Sum + geom_violin(trim = FALSE, size = 0.8) +
  theme(panel.background = element_rect(fill = "white", colour = "black",
                                        size = 2, linetype = "solid"))+
  stat_summary(fun= mean, geom = "point",
               shape = 18, size = 3, color = "black")+
  theme(axis.text.x = element_text(face="bold", # color="#993333", 
                                   size=14,color="black", angle=0,vjust=0.5),
        axis.text.y = element_text(face="bold",color="black",  size=14),
        axis.title.x = element_text(size = 14,face="bold"),
        axis.title.y = element_text(size = 14,face="bold"))+
  geom_hline(yintercept = mean(PTK2_G2M_Sum$PTK2), linetype = 2)+ # Add horizontal line at base mean
  scale_fill_manual(values=c("#538ebd", "#2e6087"))+
  stat_compare_means( aes(label = ..p.signif..), 
                      label.x = 1.5, label.y = 8,size = 7)+
  theme(legend.position="top",legend.text=element_text(size = 14),legend.title=element_text(size = 14,face="bold"))+ 
  ggtitle("G2M")+ 
  theme(plot.title = element_text(size = 18, vjust = -5 ,hjust = 0,face="bold")
  )

ggplot_PTK2_G2M_Sum_M


ggplot_PTK2_G1_Sum_M + ggplot_PTK2_S_Sum_M + ggplot_PTK2_G2M_Sum_M  
# ggplot_PTK2_G1_Sum + geom_boxplot(trim = FALSE, size = 0.8) + 
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
#   geom_hline(yintercept = mean(PTK2_G1_Sum$PTK2), linetype = 2)+ # Add horizontal line at base mean
#   scale_fill_manual(values=c("#edd15f", "#db993b"))+
#   stat_compare_means( aes(label = ..p.signif..), 
#                       label.x = 1.5, label.y = 5,size = 7)+
#   #  theme(axis.line.x = element_line(colour = "black", size = 0.5),
#   #        axis.line.y = element_line(colour = "black", size = 0.5))+
#   theme(legend.position="top",legend.text=element_text(size = 14),legend.title=element_text(size = 14,face="bold"))

