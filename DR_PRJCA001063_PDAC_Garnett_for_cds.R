############# Import files settings #############
## General setting
PathName = setwd(getwd())
RVersion = "Test"
dir.create(paste0(PathName,"/",RVersion))

## Marker genes file
# Marker_file_Name <- c("NAKAMURA_METASTASIS_MODEL_UP")
# Marker_file <- paste0(PathName,"/",Marker_file_Name,".txt")
# Marker_List <- read.delim(Marker_file,header=F,sep= c("\t"))
# Marker_List2 <- as.data.frame(Marker_List[-1:-2,])

############# Marker genes file (Old Version) #############
# Marker_file_Name <- c("NAKAMURA_METASTASIS_MODEL_M18483")
# Marker_file <- paste0(PathName,"/marker_file_",Marker_file_Name,".txt")
# Marker_List <- read.table(Marker_file,header=F,sep= c(","),stringsAsFactors = FALSE, fill = TRUE)
# library(stringr)
# Marker_List_1 <- Marker_List[2,1]
# Marker_List_2 <- str_replace_all(Marker_List_1,"expressed: ","")
# Marker_List <- str_trim(Marker_List[2,-1], side = c("both"))
# Marker_List <- c(Marker_List_2,Marker_List)
############# Marker genes file (Old Version) #############

Garnett_Marker_file_Name <- c("NAKAMURA_METASTASIS_MODEL_M18483")
Garnett_Marker_file <- paste0(PathName,"/marker_file_",Garnett_Marker_file_Name,".txt")


##############  Annotate your cells according to type (Custom Marker)  ##############

# ####################    Cell discrimination by AddModuleScore    ####################
# getFilePath("Monocle3_AddModuleScore.R")
# set.seed(1) # Fix the seed
# 
# Marker_PDAC_file_Name <- c("GRUETZMANN_PANCREATIC_CANCER_UP")
# Marker_PDAC_Name <- c("PDAC")
# cds <- Monocle3_AddModuleScore(Marker_PDAC_file_Name,Marker_PDAC_Name,marrow,cds)
# plot_cells(cds, color_cells_by= Marker_PDAC_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE)
# plot_cells(cds, color_cells_by= Marker_PDAC_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1.2) +
#   scale_colour_gradient2(low = "#440075", mid = "#ffd261", high = "#4aff8c", 
#                          guide = "colourbar",midpoint = 0.2, labs(fill = Marker_PDAC_Name))


####################   Cell discrimination by Garnett  ####################
Human_classifier_cds <- train_cell_classifier(cds = cds,
                                              marker_file = Garnett_Marker_file,   # Import the marker_file
                                              db=org.Hs.eg.db::org.Hs.eg.db, cds_gene_id_type = "SYMBOL",
                                              #num_unknown = 2215, max_training_samples = 10000,
                                              marker_file_gene_id_type = "SYMBOL",cores=8)


cds_Garnett <- classify_cells(cds, Human_classifier_cds, db = org.Hs.eg.db::org.Hs.eg.db,
                              cluster_extend = TRUE, cds_gene_id_type = "SYMBOL")

plot_cells(cds_Garnett,group_cells_by="cluster",cell_size=1.5,
           color_cells_by="cluster_ext_type", show_trajectory_graph = FALSE)

plot_cells(cds_Garnett,group_cells_by="cluster",cell_size=1.5,
           color_cells_by="cluster_ext_type",label_cell_groups=FALSE, show_trajectory_graph = FALSE)
