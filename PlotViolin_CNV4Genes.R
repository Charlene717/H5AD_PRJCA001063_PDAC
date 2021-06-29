library(ggpubr)
##### Genes ######
# CDT1
ciliated_genes_CDT1 <- c("CDT1")
cds_sub_AcinaDucT_NewK_ReCluster_CDT1 <- cds_sub_AcinaDucT_NewK_ReCluster[rowData(cds_sub_AcinaDucT_NewK_ReCluster)$gene_short_name %in% ciliated_genes_CDT1,]
plot_genes_violin(cds_sub_AcinaDucT_NewK_ReCluster_CDT1, group_cells_by="cell_cycle", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
CDT1 <- as.data.frame(cds_sub_AcinaDucT_NewK_ReCluster_CDT1@assays@data@listData[["counts"]])

# CDC6
ciliated_genes_CDC6 <- c("CDC6")
cds_sub_AcinaDucT_NewK_ReCluster_CDC6 <- cds_sub_AcinaDucT_NewK_ReCluster[rowData(cds_sub_AcinaDucT_NewK_ReCluster)$gene_short_name %in% ciliated_genes_CDC6,]
plot_genes_violin(cds_sub_AcinaDucT_NewK_ReCluster_CDC6, group_cells_by="cell_cycle", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
CDC6 <- as.data.frame(cds_sub_AcinaDucT_NewK_ReCluster_CDC6@assays@data@listData[["counts"]])

# RAD17
ciliated_genes_RAD17 <- c("RAD17")
cds_sub_AcinaDucT_NewK_ReCluster_RAD17 <- cds_sub_AcinaDucT_NewK_ReCluster[rowData(cds_sub_AcinaDucT_NewK_ReCluster)$gene_short_name %in% ciliated_genes_RAD17,]
plot_genes_violin(cds_sub_AcinaDucT_NewK_ReCluster_RAD17, group_cells_by="cell_cycle", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
RAD17 <- as.data.frame(cds_sub_AcinaDucT_NewK_ReCluster_RAD17@assays@data@listData[["counts"]])

# GMNN
ciliated_genes_GMNN <- c("GMNN")
cds_sub_AcinaDucT_NewK_ReCluster_GMNN <- cds_sub_AcinaDucT_NewK_ReCluster[rowData(cds_sub_AcinaDucT_NewK_ReCluster)$gene_short_name %in% ciliated_genes_GMNN,]
plot_genes_violin(cds_sub_AcinaDucT_NewK_ReCluster_GMNN, group_cells_by="cell_cycle", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
GMNN <- as.data.frame(cds_sub_AcinaDucT_NewK_ReCluster_GMNN@assays@data@listData[["counts"]])


###### Pheno ######
TNType <- as.character(cds_sub_AcinaDucT@colData@listData[["Type"]])


## \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
## Markers ##******************************************************************************************************************************
## CDT1 ##--------------------------------------------------------------------------------------------------------------------
CDT1_Sum <- as.data.frame(cbind(t(CDT1),ReCluster))
#!!! We need to trans the factor to numeric in here
CDT1_Sum$CDT1 <- as.numeric(CDT1_Sum$CDT1)

ggplot_CDT1 <- ggplot(data =CDT1_Sum, aes(x = ReCluster, y = CDT1, fill =ReCluster))
ggplot_CDT1 
ggplot_CDT1  + geom_violin()
ggplot_CDT1  + geom_boxplot(alpha = 0.5, show.legend = FALSE)
ggplot_CDT1 + geom_violin(trim = FALSE) +
  stat_summary(fun= mean, geom = "point",
               shape = 23, size = 2, color = "blue")

ggplot_CDT1 + geom_violin(trim = FALSE) +
  stat_summary(fun= mean, geom = "point",
               shape = 23, size = 2, color = "blue")+
  theme(axis.text.x = element_text(face="bold", # color="#993333", 
                                   size=10, angle=45),
        axis.text.y = element_text(face="bold"),
        axis.title.x = element_text(size = 14,face="bold"),
        axis.title.y = element_text(size = 14,face="bold"))+
  geom_hline(yintercept = mean(CDT1_Sum$CDT1), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(method = "anova", label.y = 10, size = 5)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", label.y = 10.5, size = 5)      # Pairwise comparison against all

## CDC6 ##--------------------------------------------------------------------------------------------------------------------
CDC6_Sum <- as.data.frame(cbind(t(CDC6),ReCluster))
#!!! We need to trans the factor to numeric in here
CDC6_Sum$CDC6 <- as.numeric(CDC6_Sum$CDC6)

ggplot_CDC6 <- ggplot(data =CDC6_Sum, aes(x = ReCluster, y = CDC6, fill =ReCluster))
ggplot_CDC6 
ggplot_CDC6  + geom_violin()
ggplot_CDC6  + geom_boxplot(alpha = 0.5, show.legend = FALSE)
ggplot_CDC6 + geom_violin(trim = FALSE) +
  stat_summary(fun= mean, geom = "point",
               shape = 23, size = 2, color = "blue")

ggplot_CDC6 + geom_violin(trim = FALSE) +
  stat_summary(fun= mean, geom = "point",
               shape = 23, size = 2, color = "blue")+
  theme(axis.text.x = element_text(face="bold", # color="#993333", 
                                   size=10, angle=45),
        axis.text.y = element_text(face="bold"),
        axis.title.x = element_text(size = 14,face="bold"),
        axis.title.y = element_text(size = 14,face="bold"))+
  geom_hline(yintercept = mean(CDC6_Sum$CDC6), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(method = "anova", label.y = 10, size = 5)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", label.y = 10.5, size = 5)      # Pairwise comparison against all

## RAD17 ##--------------------------------------------------------------------------------------------------------------------
RAD17_Sum <- as.data.frame(cbind(t(RAD17),ReCluster))
#!!! We need to trans the factor to numeric in here
RAD17_Sum$RAD17 <- as.numeric(RAD17_Sum$RAD17)

ggplot_RAD17 <- ggplot(data =RAD17_Sum, aes(x = ReCluster, y = RAD17, fill =ReCluster))
ggplot_RAD17 
ggplot_RAD17  + geom_violin()
ggplot_RAD17  + geom_boxplot(alpha = 0.5, show.legend = FALSE)
ggplot_RAD17 + geom_violin(trim = FALSE) +
  stat_summary(fun= mean, geom = "point",
               shape = 23, size = 2, color = "blue")

ggplot_RAD17 + geom_violin(trim = FALSE) +
  stat_summary(fun= mean, geom = "point",
               shape = 23, size = 2, color = "blue")+
  theme(axis.text.x = element_text(face="bold", # color="#993333", 
                                   size=10, angle=45),
        axis.text.y = element_text(face="bold"),
        axis.title.x = element_text(size = 14,face="bold"),
        axis.title.y = element_text(size = 14,face="bold"))+
  geom_hline(yintercept = mean(RAD17_Sum$RAD17), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(method = "anova", label.y = 10, size = 5)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", label.y = 10.5, size = 5)      # Pairwise comparison against all

## GMNN ##--------------------------------------------------------------------------------------------------------------------
GMNN_Sum <- as.data.frame(cbind(t(GMNN),ReCluster))
#!!! We need to trans the factor to numeric in here
GMNN_Sum$GMNN <- as.numeric(GMNN_Sum$GMNN)

ggplot_GMNN <- ggplot(data =GMNN_Sum, aes(x = ReCluster, y = GMNN, fill =ReCluster))
ggplot_GMNN 
ggplot_GMNN  + geom_violin()
ggplot_GMNN  + geom_boxplot(alpha = 0.5, show.legend = FALSE)
ggplot_GMNN + geom_violin(trim = FALSE) +
  stat_summary(fun= mean, geom = "point",
               shape = 23, size = 2, color = "blue")

ggplot_GMNN + geom_violin(trim = FALSE) +
  stat_summary(fun= mean, geom = "point",
               shape = 23, size = 2, color = "blue")+
  theme(axis.text.x = element_text(face="bold", # color="#993333", 
                                   size=10, angle=45),
        axis.text.y = element_text(face="bold"),
        axis.title.x = element_text(size = 14,face="bold"),
        axis.title.y = element_text(size = 14,face="bold"))+
  geom_hline(yintercept = mean(GMNN_Sum$GMNN), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(method = "anova", label.y = 62, size = 5)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
    
                                      ref.group = ".all.", label.y = 64.5, size = 5)      # Pairwise comparison against all


######################################  cds_sub_AcinaDucT_subset_AcinaDucT ########################################
##################  Grab specific terms ################## 
## grepl Normal
cds_sub_AcinaDucT_N <- cds_sub_AcinaDucT[,grepl("N", colData(cds_sub_AcinaDucT)$Type, ignore.case=TRUE)]
plot_cells(cds_sub_AcinaDucT_N, color_cells_by="partition")
plot_cells(cds_sub_AcinaDucT_N, color_cells_by="partition", show_trajectory_graph = F)
plot_cells(cds_sub_AcinaDucT_N, genes=c(Main),cell_size=1,label_cell_groups = FALSE, show_trajectory_graph = FALSE)
plot_cells(cds_sub_AcinaDucT_N, genes=c(Main_Group),cell_size=0.5,label_cell_groups = FALSE, show_trajectory_graph = FALSE)

## grepl Tumor
cds_sub_AcinaDucT_T <- cds_sub_AcinaDucT[,grepl("T", colData(cds_sub_AcinaDucT)$Type, ignore.case=TRUE)]
plot_cells(cds_sub_AcinaDucT_T, color_cells_by="partition", show_trajectory_graph = F)
plot_cells(cds_sub_AcinaDucT_T, color_cells_by="cluster", show_trajectory_graph = F)

###
## \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
###### Pheno ######*****************************************************************************************************************************
## TNType ##--------------------------------------------------------------------------------------------------------------------
## CDT1 ## ---------------------------------------------------------------------------------------------------------------------
TNType_Sum_CDT1 <- as.data.frame(cbind(t(CDT1),TNType))
#!!! We need to trans the factor to numeric in here
TNType_Sum_CDT1$CDT1 <- as.numeric(TNType_Sum_CDT1$CDT1)

ggplot_TNType <- ggplot(data =TNType_Sum_CDT1, aes(x =  TNType, y =CDT1 , fill =TNType))
ggplot_TNType 
ggplot_TNType  + geom_violin()
ggplot_TNType  + geom_boxplot(alpha = 0.5, show.legend = FALSE)
ggplot_TNType + geom_violin(trim = FALSE) +
  stat_summary(fun= mean, geom = "point",
               shape = 23, size = 2, color = "blue")

# https://www.r-bloggers.com/2017/06/add-p-values-and-significance-levels-to-ggplots/
ggplot_TNType + geom_violin(trim = FALSE) +
  stat_summary(fun= mean, geom = "point",
               shape = 23, size = 2, color = "blue")+
  theme(axis.text.x = element_text(face="bold", # color="#993333", 
                                   size=15, angle=45),
        axis.text.y = element_text(face="bold",size=15),
        axis.title.x = element_text(size = 14,face="bold"),
        axis.title.y = element_text(size = 14,face="bold"))+
  geom_hline(yintercept = mean(TNType_Sum_CDT1$TNType), linetype = 2)+ # Add horizontal line at base mean
   stat_compare_means( aes(label = ..p.signif..), 
                        label.x = 1.5, label.y = 8)

## CDC6 ## ---------------------------------------------------------------------------------------------------------------------
TNType_Sum_CDC6 <- as.data.frame(cbind(t(CDC6),TNType))
#!!! We need to trans the factor to numeric in here
TNType_Sum_CDC6$CDC6 <- as.numeric(TNType_Sum_CDC6$CDC6)

ggplot_TNType <- ggplot(data =TNType_Sum_CDC6, aes(x =  TNType, y =CDC6 , fill =TNType))
ggplot_TNType 
ggplot_TNType  + geom_violin()
ggplot_TNType  + geom_boxplot(alpha = 0.5, show.legend = FALSE)
ggplot_TNType + geom_violin(trim = FALSE) +
  stat_summary(fun= mean, geom = "point",
               shape = 23, size = 2, color = "blue")

# https://www.r-bloggers.com/2017/06/add-p-values-and-significance-levels-to-ggplots/
ggplot_TNType + geom_violin(trim = FALSE) +
  stat_summary(fun= mean, geom = "point",
               shape = 23, size = 2, color = "blue")+
  theme(axis.text.x = element_text(face="bold", # color="#993333", 
                                   size=15, angle=45),
        axis.text.y = element_text(face="bold",size=15),
        axis.title.x = element_text(size = 14,face="bold"),
        axis.title.y = element_text(size = 14,face="bold"))+
  geom_hline(yintercept = mean(TNType_Sum_CDC6$TNType), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means( aes(label = ..p.signif..), 
                      label.x = 1.5, label.y = 9)

## RAD17 ## ---------------------------------------------------------------------------------------------------------------------
TNType_Sum_RAD17 <- as.data.frame(cbind(t(RAD17),TNType))
#!!! We need to trans the factor to numeric in here
TNType_Sum_RAD17$RAD17 <- as.numeric(TNType_Sum_RAD17$RAD17)

ggplot_TNType <- ggplot(data =TNType_Sum_RAD17, aes(x =  TNType, y =RAD17 , fill =TNType))
ggplot_TNType 
ggplot_TNType  + geom_violin()
ggplot_TNType  + geom_boxplot(alpha = 0.5, show.legend = FALSE)
ggplot_TNType + geom_violin(trim = FALSE) +
  stat_summary(fun= mean, geom = "point",
               shape = 23, size = 2, color = "blue")

# https://www.r-bloggers.com/2017/06/add-p-values-and-significance-levels-to-ggplots/
ggplot_TNType + geom_violin(trim = FALSE) +
  stat_summary(fun= mean, geom = "point",
               shape = 23, size = 2, color = "blue")+
  theme(axis.text.x = element_text(face="bold", # color="#993333", 
                                   size=15, angle=45),
        axis.text.y = element_text(face="bold",size=15),
        axis.title.x = element_text(size = 14,face="bold"),
        axis.title.y = element_text(size = 14,face="bold"))+
  geom_hline(yintercept = mean(TNType_Sum_RAD17$TNType), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means( aes(label = ..p.signif..), 
                      label.x = 1.5, label.y = 9)

## GMNN ## ---------------------------------------------------------------------------------------------------------------------
TNType_Sum_GMNN <- as.data.frame(cbind(t(GMNN),TNType))
#!!! We need to trans the factor to numeric in here
TNType_Sum_GMNN$GMNN <- as.numeric(TNType_Sum_GMNN$GMNN)

ggplot_TNType <- ggplot(data =TNType_Sum_GMNN, aes(x =  TNType, y =GMNN , fill =TNType))
ggplot_TNType 
ggplot_TNType  + geom_violin()
ggplot_TNType  + geom_boxplot(alpha = 0.5, show.legend = FALSE)
ggplot_TNType + geom_violin(trim = FALSE) +
  stat_summary(fun= mean, geom = "point",
               shape = 23, size = 2, color = "blue")

# https://www.r-bloggers.com/2017/06/add-p-values-and-significance-levels-to-ggplots/
ggplot_TNType + geom_violin(trim = FALSE) +
  stat_summary(fun= mean, geom = "point",
               shape = 23, size = 2, color = "blue")+
  theme(axis.text.x = element_text(face="bold", # color="#993333", 
                                   size=15, angle=45),
        axis.text.y = element_text(face="bold",size=15),
        axis.title.x = element_text(size = 14,face="bold"),
        axis.title.y = element_text(size = 14,face="bold"))+
  geom_hline(yintercept = mean(TNType_Sum_GMNN$TNType), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means( aes(label = ..p.signif..), 
                      label.x = 1.5, label.y = 55)


