

# # How to define a function
# FUNCTION_NAME <- function(INPUT1, INPUT2, ..., PARAM1, PARAM2, ...) {
#   # BODY
#   return(OUTPUT)
# }


PCA_Threshold_Pos <- function(PCA_F_T3,TN){

PCA_F_T3_PC_SUM <- list()
for(i in 1:50){
  LIstPC <- as.data.frame(PCA_F_T3[,i])
  LIstPC <-cbind(row.names(LIstPC),LIstPC)
  LIstPC2 <-   as.data.frame(LIstPC[LIstPC[,2]>=0.04,])
  colnames(LIstPC2) <- c("gene",paste0("PCA_F_T",TN,"_PC",i))
  #  assign(paste0("PCA_F_T3_PC",i),as.data.frame(PCA_F_T3[,i])) 
  #Error  assign(paste0("PCA_F_T3_PC",i),PCA_F_T3[,PCA_F_T3[,i]>=0.09]) 
  #Too much#  assign(paste0("PCA_F_T3_PC",i),LIstPC2)
  
  PCA_F_T3_PC_SUM[[i]] <- as.matrix(LIstPC2[,1])
  #  PCA_F_T3_PC_SUM[[i]] <- LIstPC2[,1]
  
}

# Format convert  
PCA_F_T3_PC_SUM2 <- data.frame()

for(i in 1:length(PCA_F_T3_PC_SUM)){
  String1_1 <- paste0(">", "PCA_F_T",TN,"_PC",i,sep = " ")
  String1_2 <- paste("expressed:", PCA_F_T3_PC_SUM[[i]][1], sep = " ")
  
  for (j in 2:length(PCA_F_T3_PC_SUM[[i]])) {
    String1_2 <- paste(String1_2,PCA_F_T3_PC_SUM[[i]][j], sep = ", ")  
  }
  
  String1_3 <- ""
  
  a <- 3*i
  PCA_F_T3_PC_SUM2[a-2,1] <- String1_1
  PCA_F_T3_PC_SUM2[a-1,1] <- String1_2
  PCA_F_T3_PC_SUM2[a,1] <- String1_3
  
}

# Write
out_file2 <- paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_PCA_T",TN,"_PT.txt")
write.table(PCA_F_T3_PC_SUM2,file = out_file2, quote = F, row.names = F, col.names = F)

PCA_Threshold_Output <- PCA_F_T3_PC_SUM


return(PCA_Threshold_Output)
}


#########################################################################################

PCA_P_T1 <- marrow_sub_DucT2_TOP2ACenter_T1@reductions[["pca"]]@feature.loadings
PCA_P_T1_PC_Sum <- PCA_Threshold_Pos(PCA_P_T1,1)

PCA_P_T2 <- marrow_sub_DucT2_TOP2ACenter_T2@reductions[["pca"]]@feature.loadings
PCA_P_T2_PC_Sum <- PCA_Threshold_Pos(PCA_P_T2,2)

PCA_P_T3 <- marrow_sub_DucT2_TOP2ACenter_T3@reductions[["pca"]]@feature.loadings
PCA_P_T3_PC_Sum <- PCA_Threshold_Pos(PCA_P_T3,3)

PCA_P_T4 <- marrow_sub_DucT2_TOP2ACenter_T4@reductions[["pca"]]@feature.loadings
PCA_P_T4_PC_Sum <- PCA_Threshold_Pos(PCA_P_T4,4)

PCA_P_T5 <- marrow_sub_DucT2_TOP2ACenter_T5@reductions[["pca"]]@feature.loadings
PCA_P_T5_PC_Sum <- PCA_Threshold_Pos(PCA_P_T5,5)

PCA_P_T6 <- marrow_sub_DucT2_TOP2ACenter_T6@reductions[["pca"]]@feature.loadings
PCA_P_T6_PC_Sum <- PCA_Threshold_Pos(PCA_P_T6,6)

PCA_P_T7 <- marrow_sub_DucT2_TOP2ACenter_T7@reductions[["pca"]]@feature.loadings
PCA_P_T7_PC_Sum <- PCA_Threshold_Pos(PCA_P_T7,7)

PCA_P_T8 <- marrow_sub_DucT2_TOP2ACenter_T8@reductions[["pca"]]@feature.loadings
PCA_P_T8_PC_Sum <- PCA_Threshold_Pos(PCA_P_T8,8)

#########################################################################################



PCA_Threshold_Neg <- function(PCA_F_T3,TN){
  
  PCA_F_T3_PC_SUM <- list()
  for(i in 1:50){
    LIstPC <- as.data.frame(PCA_F_T3[,i])
    LIstPC <-cbind(row.names(LIstPC),LIstPC)
    LIstPC2 <-   as.data.frame(LIstPC[LIstPC[,2]<=-0.04,])
    colnames(LIstPC2) <- c("gene",paste0("PCA_N_T",TN,"_PC",i))
    #  assign(paste0("PCA_F_T3_PC",i),as.data.frame(PCA_F_T3[,i])) 
    #Error  assign(paste0("PCA_F_T3_PC",i),PCA_F_T3[,PCA_F_T3[,i]>=0.09]) 
    #Too much#  assign(paste0("PCA_F_T3_PC",i),LIstPC2)
    
    PCA_F_T3_PC_SUM[[i]] <- as.matrix(LIstPC2[,1])
    #  PCA_F_T3_PC_SUM[[i]] <- LIstPC2[,1]
    
  }
  
  # Format convert  
  PCA_F_T3_PC_SUM2 <- data.frame()
  
  for(i in 1:length(PCA_F_T3_PC_SUM)){
    String1_1 <- paste0(">", "PCA_N_T",TN,"_PC",i,sep = " ")
    String1_2 <- paste("expressed:", PCA_F_T3_PC_SUM[[i]][1], sep = " ")
    
    for (j in 2:length(PCA_F_T3_PC_SUM[[i]])) {
      String1_2 <- paste(String1_2,PCA_F_T3_PC_SUM[[i]][j], sep = ", ")  
    }
    
    String1_3 <- ""
    
    a <- 3*i
    PCA_F_T3_PC_SUM2[a-2,1] <- String1_1
    PCA_F_T3_PC_SUM2[a-1,1] <- String1_2
    PCA_F_T3_PC_SUM2[a,1] <- String1_3
    
  }
  
  # Write
  out_file2 <- paste0(PathName,"/",RVersion,"/",RVersion,"_","CellCycle_DucT2_TOP2ACenter_PCA_T",TN,"_NT.txt")
  write.table(PCA_F_T3_PC_SUM2,file = out_file2, quote = F, row.names = F, col.names = F)
  
  PCA_Threshold_Output <- PCA_F_T3_PC_SUM
  
  
  return(PCA_Threshold_Output)
}


########
PCA_P_T4_NC_Sum <- PCA_Threshold_Neg(PCA_P_T4,4)



#########################################################################################

PCA_P_T1 <- marrow_sub_DucT2_TOP2ACenter_T1@reductions[["pca"]]@feature.loadings
PCA_P_T1_NC_Sum <- PCA_Threshold_Neg(PCA_P_T1,1)

PCA_P_T2 <- marrow_sub_DucT2_TOP2ACenter_T2@reductions[["pca"]]@feature.loadings
PCA_P_T2_NC_Sum <- PCA_Threshold_Neg(PCA_P_T2,2)

PCA_P_T3 <- marrow_sub_DucT2_TOP2ACenter_T3@reductions[["pca"]]@feature.loadings
PCA_P_T3_NC_Sum <- PCA_Threshold_Neg(PCA_P_T3,3)

PCA_P_T4 <- marrow_sub_DucT2_TOP2ACenter_T4@reductions[["pca"]]@feature.loadings
PCA_P_T4_NC_Sum <- PCA_Threshold_Neg(PCA_P_T4,4)

PCA_P_T5 <- marrow_sub_DucT2_TOP2ACenter_T5@reductions[["pca"]]@feature.loadings
PCA_P_T5_NC_Sum <- PCA_Threshold_Neg(PCA_P_T5,5)

PCA_P_T6 <- marrow_sub_DucT2_TOP2ACenter_T6@reductions[["pca"]]@feature.loadings
PCA_P_T6_NC_Sum <- PCA_Threshold_Neg(PCA_P_T6,6)

PCA_P_T7 <- marrow_sub_DucT2_TOP2ACenter_T7@reductions[["pca"]]@feature.loadings
PCA_P_T7_NC_Sum <- PCA_Threshold_Neg(PCA_P_T7,7)

PCA_P_T8 <- marrow_sub_DucT2_TOP2ACenter_T8@reductions[["pca"]]@feature.loadings
PCA_P_T8_NC_Sum <- PCA_Threshold_Neg(PCA_P_T8,8)

#########################################################################################

