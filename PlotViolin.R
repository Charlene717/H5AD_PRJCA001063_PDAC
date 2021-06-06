# #載入gcookbook(取得男女生身高體重數據) 
# library(gcookbook) 
# # Base plot 
# p <- ggplot(heightweight, aes(x = sex, y = heightIn)) 
# p + geom_violin() 

##### geom_violin()
##### https://zhuanlan.zhihu.com/p/148189818
PDAC_Test1 <- cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["PDAC"]]
PDAC_Test1p100 <- 100*PDAC_Test1
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
