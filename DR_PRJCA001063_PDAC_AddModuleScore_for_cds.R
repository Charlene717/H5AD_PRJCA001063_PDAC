##############  Annotate your cells according to type (Custom Marker)  ##############

####################    Cell discrimination by AddModuleScore    ####################
getFilePath("Monocle3_AddModuleScore.R")
set.seed(1) # Fix the seed

Marker_PDAC_file_Name <- c("GRUETZMANN_PANCREATIC_CANCER_UP")
Marker_PDAC_Name <- c("PDAC")
cds <- Monocle3_AddModuleScore(Marker_PDAC_file_Name,Marker_PDAC_Name,marrow,cds)
plot_cells(cds, color_cells_by= Marker_PDAC_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE)
plot_cells(cds, color_cells_by= Marker_PDAC_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1.2) +
  scale_colour_gradient2(low = "#440075", mid = "#ffd261", high = "#4aff8c", 
                         guide = "colourbar",midpoint = 0.2, labs(fill = Marker_PDAC_Name))

Marker_EMT_file_Name <- c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
Marker_EMT_Name <- c("EMT")
cds <- Monocle3_AddModuleScore(Marker_EMT_file_Name,Marker_EMT_Name,marrow,cds)
plot_cells(cds, color_cells_by= Marker_EMT_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1.2) +
           scale_colour_gradient2(low = "#440075", mid = "#ffd261", high = "#4aff8c",
                                  guide = "colourbar",midpoint = 0.2, labs(fill = Marker_EMT_Name))

Marker_ACST_file_Name <- c("HP_ABNORMALITY_OF_CHROMOSOME_STABILITY")
Marker_ACST_Name <- c("ACST")
cds <- Monocle3_AddModuleScore(Marker_ACST_file_Name,Marker_ACST_Name,marrow,cds)
plot_cells(cds, color_cells_by= Marker_ACST_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1.2) +
           scale_colour_gradient2(low = "darkblue", mid = "#f7c211", high = "green",
                                  guide = "colourbar",midpoint = 0.15, labs(fill = Marker_ACST_Name))

Marker_Mig_file_Name <- c("GOBP_POSITIVE_REGULATION_OF_EPITHELIAL_CELL_MIGRATION")
Marker_Mig_Name <- c("Migration")
cds <- Monocle3_AddModuleScore(Marker_Mig_file_Name,Marker_Mig_Name,marrow,cds)
plot_cells(cds, color_cells_by= Marker_Mig_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1.2) +
           scale_colour_gradient2(low = "darkblue", mid = "#f7c211", high = "green",
                                  guide = "colourbar",midpoint = 0.15, labs(fill = Marker_Mig_Name))

Marker_Meta_file_Name <- c("NAKAMURA_METASTASIS_MODEL_UP")
Marker_Meta_Name <- c("Metastasis")
cds <- Monocle3_AddModuleScore(Marker_Meta_file_Name,Marker_Meta_Name,marrow,cds)
plot_cells(cds, color_cells_by= Marker_Meta_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1.2) +
           scale_colour_gradient2(low = "darkblue", mid = "#f7c211", high = "green",
                                  guide = "colourbar",midpoint = 0.15, labs(fill = Marker_Meta_Name))

Marker_NE_file_Name <- c("REACTOME_INITIATION_OF_NUCLEAR_ENVELOPE_NE_REFORMATION")
Marker_NE_Name <- c("NE")
cds <- Monocle3_AddModuleScore(Marker_NE_file_Name,Marker_NE_Name,marrow,cds)
plot_cells(cds, color_cells_by= Marker_NE_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1.2) +
  scale_colour_gradient2(low = "darkblue", mid = "#f7c211", high = "green",
                         guide = "colourbar",midpoint = 0.15, labs(fill = Marker_NE_Name))

Marker_NP_file_Name <- c("REACTOME_NUCLEAR_PORE_COMPLEX_NPC_DISASSEMBLY")
Marker_NP_Name <- c("NP")
cds <- Monocle3_AddModuleScore(Marker_NP_file_Name,Marker_NP_Name,marrow,cds)
plot_cells(cds, color_cells_by= Marker_NP_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1.2) +
  scale_colour_gradient2(low = "darkblue", mid = "#f7c211", high = "green",
                         guide = "colourbar",midpoint = 0.15, labs(fill = Marker_NP_Name))


Marker_ATR_file_Name <- c("REACTOME_ACTIVATION_OF_ATR_IN_RESPONSE_TO_REPLICATION_STRESS")
Marker_ATR_Name <- c("ATR")
cds <- Monocle3_AddModuleScore(Marker_ATR_file_Name,Marker_ATR_Name,marrow,cds)
plot_cells(cds, color_cells_by= Marker_ATR_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1.2) +
  scale_colour_gradient2(low = "darkblue", mid = "#f7c211", high = "green",
                         guide = "colourbar",midpoint = 0.15, labs(fill = Marker_ATR_Name))

Marker_G0G1_file_Name <- c("REACTOME_G0_AND_EARLY_G1")
Marker_G0G1_Name <- c("G0G1")
cds <- Monocle3_AddModuleScore(Marker_G0G1_file_Name,Marker_G0G1_Name,marrow,cds)
plot_cells(cds, color_cells_by= Marker_G0G1_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1.2) +
  scale_colour_gradient2(low = "darkblue", mid = "#f7c211", high = "green",
                         guide = "colourbar",midpoint = 0.15, labs(fill = Marker_G0G1_Name))

# ############    Plot Cell discrimination by AddModuleScore (AcinaDucT)   ############
# plot_cells(cds_sub_AcinaDucT_NewK_ReCluster, color_cells_by= Marker_PDAC_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1.2) +
#            scale_colour_gradient2(low = "#440075", mid = "#ffd261", high = "#4aff8c", 
#                                  guide = "colourbar",midpoint = 0.2, labs(fill = Marker_PDAC_Name))
# 
# plot_cells(cds_sub_AcinaDucT_NewK_ReCluster, color_cells_by= Marker_EMT_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1.2) +
#   scale_colour_gradient2(low = "#440075", mid = "#ffd261", high = "#4aff8c", 
#                          guide = "colourbar",midpoint = 0.2, labs(fill = Marker_EMT_Name))
# 
# plot_cells(cds_sub_AcinaDucT_NewK_ReCluster, color_cells_by= Marker_ChroSt_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1.2) +
#   scale_colour_gradient2(low = "darkblue", mid = "#f7c211", high = "green", 
#                          guide = "colourbar",midpoint = 0.15, labs(fill = Marker_ChroSt_Name))
# 
# plot_cells(cds_sub_AcinaDucT_NewK_ReCluster, color_cells_by= Marker_Mig_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1.2) +
#   scale_colour_gradient2(low = "darkblue", mid = "#f7c211", high = "green", 
#                          guide = "colourbar",midpoint = 0.15, labs(fill = Marker_Mig_Name))
# 
# plot_cells(cds_sub_AcinaDucT_NewK_ReCluster, color_cells_by= Marker_Meta_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1.2) +
#   scale_colour_gradient2(low = "darkblue", mid = "#f7c211", high = "green", 
#                          guide = "colourbar",midpoint = 0.15, labs(fill = Marker_Meta_Name))
# 
# plot_cells(cds_sub_AcinaDucT_NewK_ReCluster, color_cells_by= Marker_NE_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1.2) +
#   scale_colour_gradient2(low = "darkblue", mid = "#f7c211", high = "green", 
#                          guide = "colourbar",midpoint = 0.15, labs(fill = Marker_NE_Name))
# 
# plot_cells(cds_sub_AcinaDucT_NewK_ReCluster, color_cells_by= Marker_NP_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1.2) +
#   scale_colour_gradient2(low = "darkblue", mid = "#f7c211", high = "green", 
#                          guide = "colourbar",midpoint = 0.15, labs(fill = Marker_NP_Name))
# 
# plot_cells(cds_sub_AcinaDucT_NewK_ReCluster, color_cells_by= Marker_ATR_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1.2) +
#   scale_colour_gradient2(low = "darkblue", mid = "#f7c211", high = "green", 
#                          guide = "colourbar",midpoint = 0.15, labs(fill = Marker_ATR_Name))
# 
# plot_cells(cds_sub_AcinaDucT_NewK_ReCluster, color_cells_by= Marker_ACST_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1.2) +
#   scale_colour_gradient2(low = "darkblue", mid = "#f7c211", high = "green", 
#                          guide = "colourbar",midpoint = 0.15, labs(fill = Marker_ACST_Name))
