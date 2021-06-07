##### geom_violin()
##### https://zhuanlan.zhihu.com/p/148189818
PDAC <- cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["PDAC"]]
EMT <- cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["EMT"]]
NP <- cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["NP"]]
NE <- cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["NE"]]
ATR <- cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["ATR"]]
Migration <- cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["Migration"]]
Metastasis <- cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["Metastasis"]]
ACST <- cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["ACST"]]

ReCluster <- as.character(cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["ReCluster"]])
Cell_cycle <- as.character(cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["cell_cycle"]])

## PDAC
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

## EMT
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

## NP
NP_Sum <- as.data.frame(cbind(NP,ReCluster))
#!!! We need to trans the factor to numeric in here
NP_Sum$NP <- as.numeric(NP_Sum$NP)

ggplot_NP <- ggplot(data =NP_Sum, aes(x = ReCluster, y = NP, fill =ReCluster))
ggplot_NP 
ggplot_NP  + geom_violin()
ggplot_NP  + geom_boxplot(alpha = 0.5, show.legend = FALSE)
ggplot_NP + geom_violin(trim = FALSE) +
            stat_summary(fun= mean, geom = "point",
            shape = 23, size = 2, color = "blue")

## NE
NE_Sum <- as.data.frame(cbind(NE,ReCluster))
#!!! We need to trans the factor to numeric in here
NE_Sum$NE <- as.numeric(NE_Sum$NE)

ggplot_NE <- ggplot(data =NE_Sum, aes(x = ReCluster, y = NE, fill =ReCluster))
ggplot_NE 
ggplot_NE  + geom_violin()
ggplot_NE  + geom_boxplot(alpha = 0.5, show.legend = FALSE)
ggplot_Ne + geom_violin(trim = FALSE) +
            stat_summary(fun= mean, geom = "point",
            shape = 23, size = 2, color = "blue")

## ATR
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

## Migration
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

## Metastasis
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

## ACST
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
