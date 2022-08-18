## Ref: http://www.sthda.com/english/wiki/correlation-test-between-two-variables-in-r
##### Presetting ######
  # rm(list = ls()) # Clean variable
  memory.limit(300000)


##### Example #####
  my_data <- mtcars
  head(my_data, 6)
  
  library("ggpubr")
  ggscatter(my_data, x = "mpg", y = "wt", 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "pearson",
            xlab = "Miles/(US) gallon", ylab = "Weight (1000 lbs)")

##### Extract cell type #####
  library(tidyverse)
  scRNA.SeuObj@meta.data[["Cell_type"]] %>% unique()
  scRNA_Sub.SeuObj <- scRNA.SeuObj[,scRNA.SeuObj@meta.data[["Cell_type"]] %in% "Macrophage cell"]
  DimPlot(scRNA_Sub.SeuObj, reduction = "umap",group.by = "Cell_type")  %>% BeautifyggPlot(.,LegPos = c(0.05, 0.15))
  
##### Extract gene matrix #####
  scRNA_Sub.MT <- scRNA_Sub.SeuObj@assays[["RNA"]]@counts %>% 
                  as.data.frame() %>% t() %>% as.data.frame()
  
  library("ggpubr")
  ggscatter(scRNA_Sub.MT, x = "MARCO", y = "CD248", 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "pearson",
            xlab = "MARCO", ylab = "CD248" ,cor.coef.size = 7) %>% 
            BeautifyggPlot(.,LegPos = c(1, 0.5), OL_Thick = 1.5,
                           AxisTitleSize=1.5)
  
  ggscatter(scRNA_Sub.MT, x = "MARCO", y = "IL1B", 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "pearson",
            xlab = "MARCO", ylab = "IL1B" ,cor.coef.size = 7) %>% 
    BeautifyggPlot(.,LegPos = c(1, 0.5), OL_Thick = 1.5,
                   AxisTitleSize=1.5)

  ggscatter(scRNA_Sub.MT, x = "MARCO", y = "NSUN2", 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "pearson",
            xlab = "MARCO", ylab = "NSUN2" ,cor.coef.size = 7) %>% 
    BeautifyggPlot(.,LegPos = c(1, 0.5), OL_Thick = 1.5,
                   AxisTitleSize=1.5)

