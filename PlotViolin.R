##### geom_violin()
##### https://zhuanlan.zhihu.com/p/148189818

memory.limit(150000)

colors_cc <-c("#FF9912B3", "#2e6087", "#417034")  ## Color for Cell-Cycle


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

## \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
## Markers ##******************************************************************************************************************************
## PDAC ##--------------------------------------------------------------------------------------------------------------------
PDAC_Sum <- as.data.frame(cbind(PDAC,ReCluster))
#!!! We need to trans the factor to numeric in here
PDAC_Sum$PDAC <- as.numeric(PDAC_Sum$PDAC)

ggplot_PDAC <- ggplot(data =PDAC_Sum, aes(x = ReCluster, y = PDAC, fill =ReCluster))
ggplot_PDAC 
ggplot_PDAC  + geom_violin()
ggplot_PDAC  + geom_boxplot(alpha = 0.5, show.legend = FALSE)
ggplot_PDAC + geom_violin(trim = FALSE) +
              stat_summary(fun= mean, geom = "point",
              shape = 23, size = 2, color = "blue")

ggplot_PDAC + geom_violin(trim = FALSE) +
  stat_summary(fun= mean, geom = "point",
               shape = 23, size = 2, color = "blue")+
  theme(axis.text.x = element_text(face="bold", # color="#993333", 
                                   size=10, angle=45),
        axis.text.y = element_text(face="bold"),
        axis.title.x = element_text(size = 14,face="bold"),
        axis.title.y = element_text(size = 14,face="bold"))+
  geom_hline(yintercept = mean(PDAC_Sum$PDAC), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(method = "anova", label.y = 34000, size = 5)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", label.y = 32000, size = 5)      # Pairwise comparison against all

## PDAC Cell cycle
PDAC_Sum_CC <- as.data.frame(cbind(PDAC,ReCluster,Cell_cycle))
#!!! We need to trans the factor to numeric in here
PDAC_Sum_CC$PDAC <- as.numeric(PDAC_Sum_CC$PDAC)

ggplot_PDAC_CC <- ggplot(data =PDAC_Sum_CC, aes(x = ReCluster, y = PDAC, fill =Cell_cycle))
ggplot_PDAC_CC  + geom_violin()+ scale_fill_manual(values = colors_cc) + 
  #      labs(title = "NP")+
  theme(axis.text.x = element_text(face="bold", # color="#993333", 
                                   size=10, angle=45),
        axis.text.y = element_text(face="bold"),
        axis.title.x = element_text(size = 14,face="bold"),
        axis.title.y = element_text(size = 14,face="bold"))

## EMT ##--------------------------------------------------------------------------------------------------------------------
EMT_Sum <- as.data.frame(cbind(EMT,ReCluster))
#!!! We need to trans the factor to numeric in here
EMT_Sum$EMT <- as.numeric(EMT_Sum$EMT)

ggplot_EMT <- ggplot(data =EMT_Sum, aes(x = ReCluster, y = EMT, fill =ReCluster))
ggplot_EMT 
ggplot_EMT  + geom_violin()
ggplot_EMT  + geom_boxplot(alpha = 0.5, show.legend = FALSE)
ggplot_EMT + geom_violin(trim = FALSE) +
             stat_summary(fun= mean, geom = "point",
             shape = 23, size = 2, color = "blue")

ggplot_EMT + geom_violin(trim = FALSE) +
  stat_summary(fun= mean, geom = "point",
               shape = 23, size = 2, color = "blue")+
  theme(axis.text.x = element_text(face="bold", # color="#993333", 
                                   size=10, angle=45),
        axis.text.y = element_text(face="bold"),
        axis.title.x = element_text(size = 14,face="bold"),
        axis.title.y = element_text(size = 14,face="bold"))+
  geom_hline(yintercept = mean(EMT_Sum$EMT), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(method = "anova", label.y = 34000, size = 5)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", label.y = 32000, size = 5)      # Pairwise comparison against all

## EMT Cell cycle
EMT_Sum_CC <- as.data.frame(cbind(EMT,ReCluster,Cell_cycle))
#!!! We need to trans the factor to numeric in here
EMT_Sum_CC$EMT <- as.numeric(EMT_Sum_CC$EMT)

ggplot_EMT_CC <- ggplot(data =EMT_Sum_CC, aes(x = ReCluster, y = EMT, fill =Cell_cycle))
ggplot_EMT_CC  + geom_violin()+ scale_fill_manual(values = colors_cc) + 
  #      labs(title = "NP")+
  theme(axis.text.x = element_text(face="bold", # color="#993333", 
                                   size=10, angle=45),
        axis.text.y = element_text(face="bold"),
        axis.title.x = element_text(size = 14,face="bold"),
        axis.title.y = element_text(size = 14,face="bold"))

## NP ##--------------------------------------------------------------------------------------------------------------------
NP_Sum <- as.data.frame(cbind(NPC,ReCluster))
#!!! We need to trans the factor to numeric in here
NP_Sum$NPC <- as.numeric(NP_Sum$NPC)

ggplot_NP <- ggplot(data =NP_Sum, aes(x = ReCluster, y = NPC, fill =ReCluster))
ggplot_NP 
ggplot_NP  + geom_violin()
ggplot_NP  + geom_boxplot(alpha = 0.5, show.legend = FALSE)
ggplot_NP + geom_violin(trim = FALSE) +
            stat_summary(fun= mean, geom = "point",
            shape = 23, size = 2, color = "blue")
ggplot_NP + geom_violin(trim = FALSE) +
  stat_summary(fun= mean, geom = "point",
               shape = 23, size = 2, color = "blue")+
  theme(axis.text.x = element_text(face="bold", # color="#993333", 
                                   size=10, angle=45),
        axis.text.y = element_text(face="bold"),
        axis.title.x = element_text(size = 14,face="bold"),
        axis.title.y = element_text(size = 14,face="bold"))+
  geom_hline(yintercept = mean(NP_Sum$NPC), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(method = "anova", label.y = 34000, size = 5)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", label.y = 32000, size = 5)      # Pairwise comparison against all

## NP PTK2 ## error
NP_Sum_PTK2 <- as.data.frame(cbind(NPC,ReCluster,t(PTK2)))
#!!! We need to trans the factor to numeric in here
NP_Sum_PTK2$NPC <- as.numeric(NP_Sum_PTK2$NPC)
#NP_Sum_PTK2$PTK2 <- as.factor(NP_Sum_PTK2$PTK2)
NP_Sum_PTK2$PTK2 <- as.numeric(NP_Sum_PTK2$PTK2)

ggplot_NP <- ggplot(data =NP_Sum_PTK2, aes(x = ReCluster, y = NPC, fill =PTK2))
ggplot_NP 
ggplot_NP  + geom_violin(aes(fill=PTK2))
ggplot_NP  + geom_boxplot(alpha = 0.5, show.legend = FALSE)
ggplot_NP + geom_violin(trim = FALSE) +
  stat_summary(fun= mean, geom = "point",
               shape = 23, size = 2, color = "blue")
## NPC Cell cycle
NP_Sum_CC <- as.data.frame(cbind(NPC,ReCluster,Cell_cycle))
#!!! We need to trans the factor to numeric in here
NP_Sum_CC$NPC <- as.numeric(NP_Sum_CC$NPC)

ggplot_NP_CC <- ggplot(data =NP_Sum_CC, aes(x = ReCluster, y = NPC, fill =Cell_cycle))
ggplot_NP_CC 

ggplot_NP_CC  + geom_violin()+ scale_fill_manual(values = colors_cc) + 
          #      labs(title = "NPC")+
                theme(axis.text.x = element_text(face="bold", # color="#993333", 
                     size=10, angle=45),
                     axis.text.y = element_text(face="bold"),
                     axis.title.x = element_text(size = 14,face="bold"),
                     axis.title.y = element_text(size = 14,face="bold"))

ggplot_NP_CC  + geom_boxplot(alpha = 0.5, show.legend = FALSE)
ggplot_NP_CC + geom_violin(trim = FALSE) +
  stat_summary(fun= mean, geom = "point",
               shape = 23, size = 2, color = "blue")


## NE  ##--------------------------------------------------------------------------------------------------------------------
NE_Sum <- as.data.frame(cbind(NE,ReCluster))
#!!! We need to trans the factor to numeric in here
NE_Sum$NE <- as.numeric(NE_Sum$NE)

ggplot_NE <- ggplot(data =NE_Sum, aes(x = ReCluster, y = NE, fill =ReCluster))
ggplot_NE 
ggplot_NE  + geom_violin()
ggplot_NE  + geom_boxplot(alpha = 0.5, show.legend = FALSE)
ggplot_NE + geom_violin(trim = FALSE) +
  stat_summary(fun= mean, geom = "point",
               shape = 23, size = 2, color = "blue")

ggplot_NE + geom_violin(trim = FALSE) +
  stat_summary(fun= mean, geom = "point",
               shape = 23, size = 2, color = "blue")+
                theme(axis.text.x = element_text(face="bold", # color="#993333", 
                                                 size=10, angle=45),
                      axis.text.y = element_text(face="bold"),
                      axis.title.x = element_text(size = 14,face="bold"),
                      axis.title.y = element_text(size = 14,face="bold"))+
                geom_hline(yintercept = mean(NE_Sum$NE), linetype = 2)+ # Add horizontal line at base mean
                stat_compare_means(method = "anova", label.y = 34000, size = 5)+        # Add global annova p-value
                stat_compare_means(label = "p.signif", method = "t.test",
                                   ref.group = ".all.", label.y = 32000, size = 5)      # Pairwise comparison against all

## NE Cell cycle
NE_Sum_CC <- as.data.frame(cbind(NE,ReCluster,Cell_cycle))
#!!! We need to trans the factor to numeric in here
NE_Sum_CC$NE <- as.numeric(NE_Sum_CC$NE)

ggplot_NE_CC <- ggplot(data =NE_Sum_CC, aes(x = ReCluster, y = NE, fill =Cell_cycle))
ggplot_NE_CC  + geom_violin()+ scale_fill_manual(values = colors_cc) + 
  #      labs(title = "NP")+
  theme(axis.text.x = element_text(face="bold", # color="#993333", 
                                   size=10, angle=45),
        axis.text.y = element_text(face="bold"),
        axis.title.x = element_text(size = 14,face="bold"),
        axis.title.y = element_text(size = 14,face="bold"))

library(ggpubr)
ggplot_NE_CC  + geom_violin()+ scale_fill_manual(values = colors_cc) + 
  #      labs(title = "NP")+
  theme(axis.text.x = element_text(face="bold", # color="#993333", 
                                   size=10, angle=45),
        axis.text.y = element_text(face="bold"),
        axis.title.x = element_text(size = 14,face="bold"),
        axis.title.y = element_text(size = 14,face="bold"))+
  geom_hline(yintercept = mean(NE_Sum_CC$NE), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(method = "anova", label.y = 1600)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.")                      # Pairwise comparison against all
## ATR ##--------------------------------------------------------------------------------------------------------------------
ATR_Sum <- as.data.frame(cbind(ATR,ReCluster))
#!!! We need to trans the factor to numeric in here
ATR_Sum$ATR <- as.numeric(ATR_Sum$ATR)

ggplot_ATR <- ggplot(data =ATR_Sum, aes(x = ReCluster, y = ATR, fill =ReCluster))
ggplot_ATR 
ggplot_ATR  + geom_violin()
ggplot_ATR  + geom_boxplot(alpha = 0.5, show.legend = FALSE)
ggplot_ATR + geom_violin(trim = FALSE) +
  stat_summary(fun= mean, geom = "point",
               shape = 23, size = 2, color = "blue")

## Migration ##--------------------------------------------------------------------------------------------------------------------
Migration_Sum <- as.data.frame(cbind(Migration,ReCluster))
#!!! We need to trans the factor to numeric in here
Migration_Sum$Migration <- as.numeric(Migration_Sum$Migration)

ggplot_Migration <- ggplot(data =Migration_Sum, aes(x = ReCluster, y = Migration, fill =ReCluster))
ggplot_Migration 
ggplot_Migration  + geom_violin()
ggplot_Migration  + geom_boxplot(alpha = 0.5, show.legend = FALSE)
ggplot_Migration + geom_violin(trim = FALSE) +
  stat_summary(fun= mean, geom = "point",
               shape = 23, size = 2, color = "blue")

## Migration Cell cycle
Migration_Sum_CC <- as.data.frame(cbind(Migration,ReCluster,Cell_cycle))
#!!! We need to trans the factor to numeric in here
Migration_Sum_CC$Migration <- as.numeric(Migration_Sum_CC$Migration)

ggplot_Migration_CC <- ggplot(data =Migration_Sum_CC, aes(x = ReCluster, y = Migration, fill =Cell_cycle))
ggplot_Migration_CC  + geom_violin()+ scale_fill_manual(values = colors_cc) + 
  #      labs(title = "NP")+
  theme(axis.text.x = element_text(face="bold", # color="#993333", 
                                   size=10, angle=45),
        axis.text.y = element_text(face="bold"),
        axis.title.x = element_text(size = 14,face="bold"),
        axis.title.y = element_text(size = 14,face="bold"))

ggplot_Migration + geom_violin(trim = FALSE) +
  stat_summary(fun= mean, geom = "point",
               shape = 23, size = 2, color = "blue")+
  theme(axis.text.x = element_text(face="bold", # color="#993333", 
                                   size=10, angle=45),
        axis.text.y = element_text(face="bold"),
        axis.title.x = element_text(size = 14,face="bold"),
        axis.title.y = element_text(size = 14,face="bold"))+
  geom_hline(yintercept = mean(Migration_Sum$Migration), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(method = "anova", label.y = 34000, size = 5)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", label.y = 32000, size = 5)      # Pairwise comparison against all


## Metastasis ##--------------------------------------------------------------------------------------------------------------------
Metastasis_Sum <- as.data.frame(cbind(Metastasis,ReCluster))
#!!! We need to trans the factor to numeric in here
Metastasis_Sum$Metastasis <- as.numeric(Metastasis_Sum$Metastasis)

ggplot_Metastasis <- ggplot(data = Metastasis_Sum, aes(x = ReCluster, y = Metastasis, fill =ReCluster))
ggplot_Metastasis 
ggplot_Metastasis  + geom_violin()
ggplot_Metastasis  + geom_boxplot(alpha = 0.5, show.legend = FALSE)
ggplot_Metastasis + geom_violin(trim = FALSE) +
  stat_summary(fun= mean, geom = "point",
               shape = 23, size = 2, color = "blue")
ggplot_Metastasis + geom_violin(trim = FALSE) +
  stat_summary(fun= mean, geom = "point",
               shape = 23, size = 2, color = "blue")+
  theme(axis.text.x = element_text(face="bold", # color="#993333", 
                                   size=10, angle=45),
        axis.text.y = element_text(face="bold"),
        axis.title.x = element_text(size = 14,face="bold"),
        axis.title.y = element_text(size = 14,face="bold"))+
  geom_hline(yintercept = mean(Metastasis_Sum$Metastasis), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(method = "anova", label.y = 34000, size = 5)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", label.y = 32000, size = 5)      # Pairwise comparison against all


## Metastasis Cell cycle
Metastasis_Sum_CC <- as.data.frame(cbind(Metastasis,ReCluster,Cell_cycle))
#!!! We need to trans the factor to numeric in here
Metastasis_Sum_CC$Metastasis <- as.numeric(Metastasis_Sum_CC$Metastasis)
#Metastasis_Sum_CC$Metastasis <- as.numeric(log10(Metastasis_Sum_CC$Metastasis))

ggplot_Metastasis_CC <- ggplot(data =Metastasis_Sum_CC, aes(x = ReCluster, y = Metastasis, fill =Cell_cycle))
ggplot_Metastasis_CC  + geom_violin()+ scale_fill_manual(values = colors_cc) + 
  #      labs(title = "NP")+
  theme(axis.text.x = element_text(face="bold", # color="#993333", 
                                   size=10, angle=45),
        axis.text.y = element_text(face="bold"),
        axis.title.x = element_text(size = 14,face="bold"),
        axis.title.y = element_text(size = 14,face="bold"))

## ACST ##--------------------------------------------------------------------------------------------------------------------
ACST_Sum <- as.data.frame(cbind(ACST,ReCluster))
#!!! We need to trans the factor to numeric in here
ACST_Sum$ACST <- as.numeric(ACST_Sum$ACST)

ggplot_ACST <- ggplot(data =ACST_Sum, aes(x = ReCluster, y = ACST, fill =ReCluster))
ggplot_ACST 
ggplot_ACST  + geom_violin()
ggplot_ACST  + geom_boxplot(alpha = 0.5, show.legend = FALSE)
ggplot_ACST + geom_violin(trim = FALSE) +
  stat_summary(fun= mean, geom = "point",
               shape = 23, size = 2, color = "blue")
ggplot_ACST + geom_violin(trim = FALSE) +
  stat_summary(fun= mean, geom = "point",
               shape = 23, size = 2, color = "blue")+
  theme(axis.text.x = element_text(face="bold", # color="#993333", 
                                   size=10, angle=45),
        axis.text.y = element_text(face="bold"),
        axis.title.x = element_text(size = 14,face="bold"),
        axis.title.y = element_text(size = 14,face="bold"))+
  geom_hline(yintercept = mean(ACST_Sum$ACST), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(method = "anova", label.y = 34000, size = 5)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", label.y = 32000, size = 5)      # Pairwise comparison against all


## ACST Cell cycle
ACST_Sum_CC <- as.data.frame(cbind(ACST,ReCluster,Cell_cycle))
#!!! We need to trans the factor to numeric in here
ACST_Sum_CC$ACST <- as.numeric(ACST_Sum_CC$ACST)
ACST_Sum_CC$ACST <- as.numeric(log2(ACST_Sum_CC$ACST))

ggplot_ACST_CC <- ggplot(data =ACST_Sum_CC, aes(x = ReCluster, y = ACST, fill =Cell_cycle))
ggplot_ACST_CC  + geom_violin()+ scale_fill_manual(values = colors_cc) + 
  #      labs(title = "ACST")+
  theme(axis.text.x = element_text(face="bold", # color="#993333", 
                                   size=10, angle=45),
        axis.text.y = element_text(face="bold"),
        axis.title.x = element_text(size = 14,face="bold"),
        axis.title.y = element_text(size = 14,face="bold"))


## Genes ##******************************************************************************************************************************
## PTK2 Cell cycle  ##--------------------------------------------------------------------------------------------------------------------
PTK2_Sum_CC <- as.data.frame(cbind(t(PTK2),ReCluster,Cell_cycle))

PTK2_Sum_CC$PTK2 <- as.numeric(PTK2_Sum_CC$PTK2) #!!! We need to trans the factor to numeric in here
PTK2_Sum_CC$PTK2 <- as.numeric(log10(PTK2_Sum_CC$PTK2)) #!!! We need to trans the factor to numeric in here


ggplot_PTK2_CC <- ggplot(data =PTK2_Sum_CC, aes(x = ReCluster, y = PTK2, fill =Cell_cycle))
ggplot_PTK2_CC  + geom_violin()+ scale_fill_manual(values = colors_cc) +
                  #      labs(title = "NP")+
                  theme(axis.text.x = element_text(face="bold", # color="#993333", 
                        size=10, angle=45),
                        axis.text.y = element_text(face="bold"),
                        axis.title.x = element_text(size = 14,face="bold"),
                        axis.title.y = element_text(size = 14,face="bold"))

ggplot_PTK2_CC  + geom_boxplot(alpha = 0.5, show.legend = FALSE)+ 
                  scale_fill_manual(values = colors_cc)
ggplot_PTK2_CC + geom_violin(trim = FALSE) +
  stat_summary(fun= mean, geom = "point",
               shape = 23, size = 2, color = "blue")

###
## TOP2A Cell cycle  ##--------------------------------------------------------------------------------------------------------------------
TOP2A_Sum_CC <- as.data.frame(cbind(t(TOP2A),ReCluster,Cell_cycle))

TOP2A_Sum_CC$TOP2A <- as.numeric(TOP2A_Sum_CC$TOP2A) #!!! We need to trans the factor to numeric in here
TOP2A_Sum_CC$TOP2A <- as.numeric(log2(TOP2A_Sum_CC$TOP2A)) #!!! We need to trans the factor to numeric in here

ggplot_TOP2A_CC <- ggplot(data =TOP2A_Sum_CC, aes(x = ReCluster, y = TOP2A, fill =Cell_cycle))
ggplot_TOP2A_CC  + geom_violin()+ scale_fill_manual(values = colors_cc) +
  #      labs(title = "NP")+
  theme(axis.text.x = element_text(face="bold", # color="#993333", 
                                   size=10, angle=45),
        axis.text.y = element_text(face="bold"),
        axis.title.x = element_text(size = 14,face="bold"),
        axis.title.y = element_text(size = 14,face="bold"))

ggplot_TOP2A_CC  + geom_boxplot(alpha = 0.8, show.legend = TRUE)+ 
                   scale_fill_manual(values = colors_cc)+
                  theme(axis.text.x = element_text(face="bold", # color="#993333", 
                                                   size=10, angle=45),
                        axis.text.y = element_text(face="bold"),
                        axis.title.x = element_text(size = 14,face="bold"),
                        axis.title.y = element_text(size = 14,face="bold"))


## CGAS Cell cycle  ##--------------------------------------------------------------------------------------------------------------------
CGAS_Sum_CC <- as.data.frame(cbind(t(CGAS),ReCluster,Cell_cycle))

CGAS_Sum_CC$CGAS<- as.numeric(CGAS_Sum_CC$CGAS) #!!! We need to trans the factor to numeric in here
CGAS_Sum_CC$CGAS<- as.numeric(log10(CGAS_Sum_CC$CGAS)) #!!! We need to trans the factor to numeric in here

ggplot_CGAS_CC <- ggplot(data =CGAS_Sum_CC, aes(x = ReCluster, y = CGAS, fill =Cell_cycle))
ggplot_CGAS_CC  + geom_violin()+ scale_fill_manual(values = colors_cc) +
  #      labs(title = "NP")+
  theme(axis.text.x = element_text(face="bold", # color="#993333", 
                                   size=10, angle=45),
        axis.text.y = element_text(face="bold"),
        axis.title.x = element_text(size = 14,face="bold"),
        axis.title.y = element_text(size = 14,face="bold"))

ggplot_CGAS_CC   + geom_boxplot(alpha = 0.8, show.legend = TRUE)+ 
                   scale_fill_manual(values = colors_cc)+
                  theme(axis.text.x = element_text(face="bold", # color="#993333", 
                                      size=10, angle=45),
                        axis.text.y = element_text(face="bold"),
                        axis.title.x = element_text(size = 14,face="bold"),
                        axis.title.y = element_text(size = 14,face="bold"))


######################################################################################


PDAC_Test1 <- cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["PDAC"]]
PDAC_Test1p100 <- 1*PDAC_Test1
# # We need to trans the factor to numeric
# PDAC_Test1p100 <- as.numeric(round(PDAC_Test1p100,1))

# We need to trans the factor to character
PDAC_Test2 <- as.character(cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["ReCluster"]])
PDAC_Test3 <- as.character(cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["cell_cycle"]])

PDAC_TestSum <- as.data.frame(cbind(PDAC_Test1p100,PDAC_Test2,PDAC_Test3))
#!!! We need to trans the factor to numeric in here
PDAC_TestSum$PDAC_Test1p100 <- as.numeric(PDAC_TestSum$PDAC_Test1p100)
# e <- ggplot(data =PDAC_TestSum, aes(x = PDAC_Test3, y = PDAC_Test1, fill =PDAC_Test2))
e <- ggplot(data =PDAC_TestSum, aes(x = PDAC_Test3, y = PDAC_Test1p100, fill =PDAC_Test2))
e
e + geom_violin()
e + geom_boxplot(alpha = 0.5, show.legend = FALSE)
e + geom_violin(trim = FALSE) +
  stat_summary(fun= mean, geom = "point",
               shape = 23, size = 2, color = "blue")

###
PDAC_TestSum2 <- as.data.frame(cbind(PDAC_Test1p100,PDAC_Test3))
PDAC_TestSum2$PDAC_Test1p100 <- as.numeric(PDAC_TestSum2$PDAC_Test1p100)

e2 <- ggplot(PDAC_TestSum2, aes(x = PDAC_Test3, y = PDAC_Test1p100))
e2
e2 + geom_violin()
e2 + geom_violin(trim = FALSE) +
  stat_summary(fun= mean, geom = "point",
               shape = 23, size = 2, color = "blue")

###
PDAC_TestSum3 <- as.data.frame(cbind(PDAC_Test1p100,PDAC_Test2, fill =PDAC_Test2))
PDAC_TestSum3$PDAC_Test1p100 <- as.numeric(PDAC_TestSum3$PDAC_Test1p100)

e3 <- ggplot(PDAC_TestSum3, aes(x = PDAC_Test2, y = PDAC_Test1p100))
e3
e3 + geom_violin(fill = "steelblue")
e3 + geom_violin(trim = FALSE) +
  stat_summary(fun= mean, geom = "point",
               shape = 23, size = 2, color = "blue")



# #### Test small1 ###
# PDAC_TestSum2_2 <- PDAC_TestSum2[1:9,]
# PDAC_TestSum2_3 <-as.matrix(PDAC_TestSum2_2)
# e3 <- ggplot(PDAC_TestSum2_2, aes(x = PDAC_Test3, y = PDAC_Test1p100))
# e3
# e3 + geom_violin()

# 
# #### Test small2 ###
# ##### https://zhuanlan.zhihu.com/p/148189818
# PDAC_Test1 <- cds_sub_DucT2_TOP2ACenter@colData@listData[["PDAC"]]
# PDAC_Test1p100 <- 100*PDAC_Test1
# PDAC_Test1p100 <- round(PDAC_Test1p100,1)
# PDAC_Test2 <- as.character(cds_sub_DucT2_TOP2ACenter@colData@listData[["ReCluster"]])
# PDAC_Test3 <- as.character(cds_sub_DucT2_TOP2ACenter@colData@listData[["cell_cycle"]])
# 
# PDAC_TestSum <- as.data.frame(cbind(PDAC_Test1p100,PDAC_Test2,PDAC_Test3))
# # e <- ggplot(data =PDAC_TestSum, aes(x = PDAC_Test3, y = PDAC_Test1, fill =PDAC_Test2))
# e <- ggplot(data =PDAC_TestSum, aes(x = PDAC_Test3, y = PDAC_Test1p100, fill =PDAC_Test2))
# e
# e + geom_violin()
# e + geom_violin(trim = FALSE) +
#   stat_summary(fun= mean, geom = "point",
#                shape = 23, size = 2, color = "blue")
# 
# 
# PDAC_TestSum2 <- as.data.frame(cbind(PDAC_Test1p100,PDAC_Test3))
# e2 <- ggplot(PDAC_TestSum2, aes(x = PDAC_Test3, y = PDAC_Test1p100))
# e2
# e2 + geom_violin()
# e2 + geom_violin(trim = FALSE) +
#   stat_summary(fun= mean, geom = "point",
#                shape = 23, size = 2, color = "blue")
