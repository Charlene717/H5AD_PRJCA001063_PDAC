Test4_2 <- as.data.frame(DeconOutputMatrix)

## Link to TCGA PAAD Phenotype information
Bulk_TCGA_PAAD_Pheno_clinicalMatrix <- read.delim(paste0(PathName,"/Deconvolution/TCGA_PAAD/PAAD_clinicalMatrix"),header=T)
Bulk_TCGA_PAAD_Pheno_survival <- read.delim(paste0(PathName,"/Deconvolution/TCGA_PAAD/TCGA_PAAD_survival.txt"),header=T)

Test4_2_M1 <- cbind(row.names(Test4_2),Test4_2)
colnames(Test4_2_M1)[[1]] <- c("sample")
library(dplyr)

Test4_2_M1_Pheno1 <- left_join(Test4_2_M1,Bulk_TCGA_PAAD_Pheno_survival,by="sample")

Test4_2_Index1 <- as.data.frame(Test4_2_M1_Pheno1$CoreCD00)
Test4_2_M1_Pheno2 <- cbind(Test4_2_M1_Pheno1,Test4_2_Index1)
colnames(Test4_2_M1_Pheno2)[[length(Test4_2_M1_Pheno2[1,])]] <- c("Index1")
#Test4_2_Index1_Mean <- mean(Test4_2_M1_Pheno2$`NK_Neg/Pos`)
#Test4_2_Index1_SD <- sd(Test4_2_M1_Pheno2$`NK_Neg/Pos`)

Test4_2_Index1_Mean <- mean(Test4_2_M1_Pheno1$CoreCD00)
Test4_2_Index1_SD <- sd(Test4_2_M1_Pheno1$CoreCD00)


## Mean & SD
Test4_2_M1_Pheno2$CoreCD00_Index <- ""
for (i in c(1:length(Test4_2_M1_Pheno2[,1]))) {
  if(Test4_2_M1_Pheno2$CoreCD00[i] >= Test4_2_Index1_Mean+Test4_2_Index1_SD){
    Test4_2_M1_Pheno2$CoreCD00_Index[i] <- 1
  }else if(Test4_2_M1_Pheno2$CoreCD00[i] <= Test4_2_Index1_Mean-Test4_2_Index1_SD){
    Test4_2_M1_Pheno2$CoreCD00_Index[i] <- 2
  }else{
    Test4_2_M1_Pheno2$CoreCD00_Index[i] <- 3
  }
}
Test4_2_M1_Pheno2 <- Test4_2_M1_Pheno2[!Test4_2_M1_Pheno2$CoreCD00_Index==3,]


# https://blog.yjtseng.info/post/2020-05-13-survival-curve/
# http://www.sthda.com/english/wiki/survminer-r-package-survival-data-analysis-and-visualization
# https://www.rdocumentation.org/packages/survminer/versions/0.4.9/topics/ggsurvplot
library(survival)
library(survminer)
Test4_2_M1_Pheno3 <- as.data.frame(cbind(Test4_2_M1_Pheno2$sample,Test4_2_M1_Pheno2$CoreCD00_Index,Test4_2_M1_Pheno2$OS,Test4_2_M1_Pheno2$OS.time))
Test4_2_M1_Pheno3 <- Test4_2_M1_Pheno2[,c("sample","CoreCD00_Index","OS","OS.time","DSS","DSS.time","DFI","DFI.time","PFI","PFI.time")]

# OS
fit_OS <- survfit(Surv(OS.time, OS) ~ CoreCD00_Index, data = Test4_2_M1_Pheno3)
fit_OS_ggplot <-  ggsurvplot(fit_OS, pval = TRUE, 
                             risk.table = TRUE,
                             legend.title = "CoreCD00_Index",
                             legend.labs = c("CoreCD00_High", "CoreCD00_Low"),
                             title="OS",
                             log.rank.weights="1")
fit_OS_ggplot

# DSS
fit_DSS <- survfit(Surv(DSS.time, DSS) ~ CoreCD00_Index, data = Test4_2_M1_Pheno3)
fit_DSS_ggplot <- ggsurvplot(fit_DSS, pval = TRUE, 
                             risk.table = TRUE,
                             legend.title = "CoreCD00_Index",
                             legend.labs = c("CoreCD00_High", "CoreCD00_Low"),
                             title="DSS",
                             log.rank.weights="1")
fit_DSS_ggplot

# DFI
fit_DFI <- survfit(Surv(DFI.time, DFI) ~ CoreCD00_Index, data = Test4_2_M1_Pheno3)
fit_DFI_ggplot <- ggsurvplot(fit_DFI, pval = TRUE, 
                             risk.table = TRUE,
                             legend.title = "CoreCD00_Index",
                             legend.labs = c("CoreCD00_High", "CoreCD00_Low"),
                             title="DFI",
                             log.rank.weights="1")
fit_DFI_ggplot

# PFI
fit_PFI <- survfit(Surv(PFI.time, PFI) ~ CoreCD00_Index, data = Test4_2_M1_Pheno3)
fit_PFI_ggplot <- ggsurvplot(fit_PFI, pval = TRUE, 
                             risk.table = TRUE,
                             legend.title = "CoreCD00_Index",
                             legend.labs = c("CoreCD00_High", "CoreCD00_Low"),
                             title="PFI",
                             log.rank.weights="1")
fit_PFI_ggplot

arrange_ggsurvplots(list(fit_OS_ggplot, fit_DSS_ggplot,fit_DFI_ggplot,fit_PFI_ggplot),
                    ncol=2, nrow=2, 
                    risk.table.height = 0.3)
