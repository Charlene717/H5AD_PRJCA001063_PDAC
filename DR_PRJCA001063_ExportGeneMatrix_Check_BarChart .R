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
## https://ithelp.ithome.com.tw/articles/10209002
colors_cc <-c("#FF9912B3", "#2e6087", "#417034")  ## Color for Cell-Cycle

annotation_col <- as.data.frame(cbind(cell_cycle,ReCluster))
ggplot(annotation_col, aes(x = ReCluster, fill = cell_cycle)) + 
  geom_bar(position = "fill")+
  theme(axis.text.x = element_text(face="bold", # color="#993333", 
                                   size=10, angle=75,vjust =0.55),
        axis.text.y = element_text(face="bold",hjust =-8),
        axis.title.x = element_text(size = 16,face="bold"),
        axis.title.y = element_text(size = 16,face="bold"),
        plot.title = element_text(color="black", size=20, 
                                  face="bold.italic",hjust = 0.1,vjust =-8), # margin = margin(t = 0.5, b = -7),
        legend.title = element_text(size=13, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"))+ 
  scale_fill_manual(values = colors_cc)+ 
  theme(plot.title = element_text(size = 18,vjust = 0.5,hjust = 0.5)
  )

ggplot(annotation_col, aes(x = ReCluster, fill = cell_cycle))+ 
  geom_bar(position = "dodge")+
  theme(axis.text.x = element_text(face="bold", # color="#993333", 
                                   size=10, angle=75,vjust =0.55),
        axis.text.y = element_text(face="bold"),
        axis.title.x = element_text(size = 16,face="bold"),
        axis.title.y = element_text(size = 16,face="bold"),
        plot.title = element_text(color="black", size=20, 
                                  face="bold.italic",hjust = 0.1,vjust =-8), # margin = margin(t = 0.5, b = -7),
        legend.title = element_text(size=13, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"))+ 
  scale_fill_manual(values = colors_cc)+ 
  theme(plot.title = element_text(size = 18, vjust = 0.5,hjust = 0.5)
  )

ggplot(annotation_col, aes(x = ReCluster, fill = cell_cycle))+ 
  geom_bar(position = "stack")+
  theme(axis.text.x = element_text(face="bold", # color="#993333", 
                                   size=10, angle=75,vjust =0.55),
        axis.text.y = element_text(face="bold"),
        axis.title.x = element_text(size = 16,face="bold"),
        axis.title.y = element_text(size = 16,face="bold"),
        plot.title = element_text(color="black", size=20, 
                                  face="bold.italic",hjust = 0.1,vjust =-8), # margin = margin(t = 0.5, b = -7),
        legend.title = element_text(size=13, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"))+ 
  scale_fill_manual(values = colors_cc)+ 
  theme(plot.title = element_text(size = 18, vjust = 0.5,hjust = 0.5)
  )


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

# Combine
GeneExpMatrix_DF4 <- cbind(Cell,Patient,Type,Cell_type,ReCluster,PDAC,EMT,Migration,Metastasis,NE,NPC,ACST,cell_cycle,GeneExpMatrix_DF2)
GeneExpMatrix_DF4 <- cbind(row.names(GeneExpMatrix_DF4),GeneExpMatrix_DF4)
colnames(GeneExpMatrix_DF4)[1] <- c("Barcode")

PathName = setwd(getwd())
RVersion = "20210626V1_GEM"
dir.create(paste0(PathName,"/",RVersion))

# write.table(GeneExpMatrix_DF4, file = paste0(PathName,"/AcinaDucT_GeneExpMatrix_Pheno.txt"),sep = " ",header= T, quote = FALSE, na = "NA")
write.table(GeneExpMatrix_DF4, file=paste0(PathName,"/",RVersion,"/AcinaDucT_GeneExpMatrix_Pheno_Check_Marker.txt"),sep="\t", row.names=F)

# GeneExpMatrix_DF2 <- cbind(row.names(GeneExpMatrix_DF),GeneExpMatrix_DF)
# colnames(GeneExpMatrix_DF2)[1] <- c("gene_short_name")


#-------------------------------------------------------------------------#
##----------------------------TOP2A high bar-------------------------------##
# https://ithelp.ithome.com.tw/articles/10209002 # R語言_ggplot2長條圖

GeneExpMatrix_DF5 <- t(as.data.frame(GeneExpMatrix_DF4))
# https://weitinglin.com/2016/01/08/r-package-dplyr-%E7%84%A1%E7%97%9B%E5%88%86%E6%9E%90dataframe/ # R package: dplyr 無痛分析dataframe
# devtools::install_github('cole-trapnell-lab/monocle3')
library(dplyr)

##---------------- TOP2A > 0 ----------------##
GeneExpMatrix_DF5_TOP2A_T <- GeneExpMatrix_DF5[,GeneExpMatrix_DF5["TOP2A",] > 0]
GeneExpMatrix_DF5_TOP2A_T_t <- as.data.frame(t(GeneExpMatrix_DF5_TOP2A_T))

colors_cc <-c("#FF9912B3", "#2e6087", "#417034")  ## Color for Cell-Cycle
# https://ithelp.ithome.com.tw/articles/10209002
ggplot(GeneExpMatrix_DF5_TOP2A_T_t, aes(x = ReCluster, fill = cell_cycle)) + 
  geom_bar(position = "fill")+
  theme(axis.text.x = element_text(face="bold", # color="#993333", 
                                   size=10, angle=75,vjust =0.55),
        axis.text.y = element_text(face="bold",hjust =-8),
        axis.title.x = element_text(size = 16,face="bold"),
        axis.title.y = element_text(size = 16,face="bold"),
        plot.title = element_text(color="black", size=20, 
                                  face="bold.italic",hjust = 0.1,vjust =-8), # margin = margin(t = 0.5, b = -7),
        legend.title = element_text(size=13, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"))+ 
        scale_fill_manual(values = colors_cc)+ 
        ggtitle("TOP2A > 0")+ 
        theme(plot.title = element_text(size = 18,vjust = 0.5,hjust = 0.5)
        )

ggplot(GeneExpMatrix_DF5_TOP2A_T_t, aes(x = ReCluster, fill = cell_cycle))+ 
  geom_bar(position = "dodge")+
  theme(axis.text.x = element_text(face="bold", # color="#993333", 
                                   size=10, angle=75,vjust =0.55),
        axis.text.y = element_text(face="bold"),
        axis.title.x = element_text(size = 16,face="bold"),
        axis.title.y = element_text(size = 16,face="bold"),
        plot.title = element_text(color="black", size=20, 
                                  face="bold.italic",hjust = 0.1,vjust =-8), # margin = margin(t = 0.5, b = -7),
        legend.title = element_text(size=13, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"))+ 
  scale_fill_manual(values = colors_cc)+ 
  ggtitle("TOP2A > 0")+ 
  theme(plot.title = element_text(size = 18, vjust = 0.5,hjust = 0.5)
  )

ggplot(GeneExpMatrix_DF5_TOP2A_T_t, aes(x = ReCluster, fill = cell_cycle))+ 
  geom_bar(position = "stack")+
  theme(axis.text.x = element_text(face="bold", # color="#993333", 
                                   size=10, angle=75,vjust =0.55),
        axis.text.y = element_text(face="bold"),
        axis.title.x = element_text(size = 16,face="bold"),
        axis.title.y = element_text(size = 16,face="bold"),
        plot.title = element_text(color="black", size=20, 
                                  face="bold.italic",hjust = 0.1,vjust =-8), # margin = margin(t = 0.5, b = -7),
        legend.title = element_text(size=13, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"))+ 
  scale_fill_manual(values = colors_cc)+ 
  ggtitle("TOP2A > 0")+ 
  theme(plot.title = element_text(size = 18, vjust = 0.5,hjust = 0.5)
  )


##---------------- TOP2A = 0 ----------------##
GeneExpMatrix_DF5_TOP2A_F <- GeneExpMatrix_DF5[,GeneExpMatrix_DF5["TOP2A",] == 0]
GeneExpMatrix_DF5_TOP2A_F_t <- as.data.frame(t(GeneExpMatrix_DF5_TOP2A_F))

colors_cc <-c("#FF9912B3", "#2e6087", "#417034")  ## Color for Cell-Cycle
ggplot(GeneExpMatrix_DF5_TOP2A_F_t, aes(x = ReCluster, fill = cell_cycle)) + 
  geom_bar(position = "fill")+
  theme(axis.text.x = element_text(face="bold", # color="#993333", 
                                   size=10, angle=75,vjust =0.55),
        axis.text.y = element_text(face="bold",hjust =-8),
        axis.title.x = element_text(size = 16,face="bold"),
        axis.title.y = element_text(size = 16,face="bold"),
        plot.title = element_text(color="black", size=20, 
                                  face="bold.italic",hjust = 0.1,vjust =-8), # margin = margin(t = 0.5, b = -7),
        legend.title = element_text(size=13, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"))+ 
  scale_fill_manual(values = colors_cc)+ 
  ggtitle("TOP2A > 0")+ 
  theme(plot.title = element_text(size = 18,vjust = 0.5,hjust = 0.5)
  )


ggplot(GeneExpMatrix_DF5_TOP2A_F_t, aes(x = ReCluster, fill = cell_cycle))+ 
  geom_bar(position = "dodge")+
  theme(axis.text.x = element_text(face="bold", # color="#993333", 
                                   size=10, angle=75,vjust =0.55),
        axis.text.y = element_text(face="bold"),
        axis.title.x = element_text(size = 16,face="bold"),
        axis.title.y = element_text(size = 16,face="bold"),
        plot.title = element_text(color="black", size=20, 
                                  face="bold.italic",hjust = 0.1,vjust =-8), # margin = margin(t = 0.5, b = -7),
        legend.title = element_text(size=13, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"))+ 
  scale_fill_manual(values = colors_cc)+ 
  ggtitle("TOP2A > 0")+ 
  theme(plot.title = element_text(size = 18, vjust = 0.5,hjust = 0.5)
  )

ggplot(GeneExpMatrix_DF5_TOP2A_F_t, aes(x = ReCluster, fill = cell_cycle))+ 
  geom_bar(position = "stack")+
  theme(axis.text.x = element_text(face="bold", # color="#993333", 
                                   size=10, angle=75,vjust =0.55),
        axis.text.y = element_text(face="bold"),
        axis.title.x = element_text(size = 16,face="bold"),
        axis.title.y = element_text(size = 16,face="bold"),
        plot.title = element_text(color="black", size=20, 
                                  face="bold.italic",hjust = 0.1,vjust =-8), # margin = margin(t = 0.5, b = -7),
        legend.title = element_text(size=13, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"))+ 
  scale_fill_manual(values = colors_cc)+ 
  ggtitle("TOP2A > 0")+ 
  theme(plot.title = element_text(size = 18, vjust = 0.5,hjust = 0.5)
  )



##---------------- TOP2A < 0 ----------------##
GeneExpMatrix_DF5_TOP2A_S <- GeneExpMatrix_DF5[,GeneExpMatrix_DF5["TOP2A",] < 0]
GeneExpMatrix_DF5_TOP2A_S_t <- as.data.frame(t(GeneExpMatrix_DF5_TOP2A_S))

colors_cc <-c("#FF9912B3", "#2e6087", "#417034")  ## Color for Cell-Cycle
ggplot(GeneExpMatrix_DF5_TOP2A_S_t, aes(x = ReCluster, fill = cell_cycle)) + 
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



# #------#----統計起來為了做長條圖----#
# # https://ithelp.ithome.com.tw/articles/10209002
# library(plyr)
# ess2 = ddply(GeneExpMatrix_DF5_TOP2A_T_t,.(ReCluster),function(.){
#   res = prop.table(table(factor(.$cell_cycle)))
#   res2 = table(factor(.$cell_cycle))
#   data.frame(lab=names(res), y=c(res),yy =c(res2))
# })
# detach("package:plyr", unload=TRUE)
# 
# for_show = GeneExpMatrix_DF5_TOP2A_T_t %>% group_by(ReCluster) %>% summarise(length(ReCluster))
# 
#  windowsFonts(A=windowsFont("微軟正黑體"))
# ggplot(ess2,aes(x = ReCluster,y=y,fill = lab))+
#   geom_bar(stat = "identity") + 
#   theme_classic(base_size = 16)+
#   geom_text(mapping = aes(label = sprintf("%.2f%%",y*100)),
#             size = 5, colour = 'black', vjust = 2, hjust = .5, position = position_stack())+
#   annotate('text',x = 1:nrow(for_show),y=0.1,label= unlist(for_show[,2]),family = "A",size = 5.5,fontface =2)+
#   
#   theme(axis.text.y=element_text(face="bold",size=15,color="#333333"))+##調整y軸字型
#   theme(axis.text.x=element_text(face="bold",size=20,angle=360,color="#333333"))+#x軸字型(筆數的數字)
#   scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1) ,labels =c("0%","25%","50%","75%","100%"))
