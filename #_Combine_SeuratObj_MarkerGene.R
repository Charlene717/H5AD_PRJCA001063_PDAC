##### Plot #####
FeaturePlot(scRNA.SeuObj, features = c("TOP2A"))
# DimPlot(scRNA.SeuObj, reduction = "umap")
DimPlot(scRNA.SeuObj, reduction = "umap",label = T)
DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "Cell_type")
DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "DataSetID")


# gut endoderm
FeaturePlot(scRNA.SeuObj, features = c("HNF1B","HNF6","FOXA2","HNF4A","HEX","GATA4","GATA6"))

# Pre-pancreatic endoderm
FeaturePlot(scRNA.SeuObj, features = c("MNX1","PTF1A","PDX1","HNF1B","HNF6","FOXA2","HNF4",
                                       "HEX","GATA4","GATA6"))

# budding endoderm
FeaturePlot(scRNA.SeuObj, features = c("SOX9","NKX6-1","NKX2-2","MNX1","PTF1A","PDX1","HNF1B",
                                       "HNF6","FOXA2","HNF4","HEX","GATA4","GATA6"))

# prodo-differentiated epithelium
FeaturePlot(scRNA.SeuObj, features = c("HES1", "PROX1", "SOX9", "NKX6-1", "NKX2-2", "MNX1", 
                                       "PDX1", "HNF1B", "HNF6", "FOXA2", "HNF4A", "GATA4", "GATA6"))

# Mult-Potent Cells
FeaturePlot(scRNA.SeuObj, features = c("Ptfla-hi", "PTF1A", "CPA1", "c-Myc-hi", "MYC", "NR5A2",
                                       "HNF1B", "MNX1", "SOX9"))


# pro-acinar
FeaturePlot(scRNA.SeuObj, features = c("GATA4", "NR5A2", "Ptfla-hi", "PTF1A"))

# precursor cells in tubules
FeaturePlot(scRNA.SeuObj, features = c("SOX9", "HNF1B", "HNF6", "HES1", "NKX6-1", "NKX2-2", "PROX1"))

# Duct cell
FeaturePlot(scRNA.SeuObj, features = c("HNF1B", "HNF6", "FOXA2", "NR5A2", "GLIS3", "PROX1"))

# Acinar cell
FeaturePlot(scRNA.SeuObj, features = c("BHLHA15", "RBPJL", "GATA4", "GATA6", "NR5A2", "Ptf1a-hi"))

# islet precursor
FeaturePlot(scRNA.SeuObj, features = c("NEUROG3", "Sox9-lo", "NKX6-1", "NKX2-2", "PDX1", "HNF6", 
                                       "FOXA2","HNF4A", "GATA6", "INSM1", "SNAI2", "SNAI1", "ISL1", 
                                       "NEUROD1", "MAFB", "GLIS3", "RFX6"))


## cancer stem cell markers
## Ref: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5737740/#:~:text=Similar%20to%20stemness%2Drelated%20transcription,%2C%20CD117%2Fc%2Dkit%2C
FeaturePlot(scRNA.SeuObj, features = c("OCT4","SOX2","KLF4","SALL4"))
FeaturePlot(scRNA.SeuObj, features = c("ALDH1A1","CD24","CD44","CD133","CXCR4"))

## PHH
FeaturePlot(scRNA.SeuObj, features = c("ALDH1A1","ALDH1A2","ALDH1A3","NR5A2","KRT19","SOX9"))
