########################  DucT2 (choose_cells) ##########################
    cds_sub_DucT2 <- choose_cells(cds)
    #cds_subset <- reduce_dimension(cds_subset)

    plot_cells(cds_sub_DucT2, label_cell_groups=FALSE, show_trajectory_graph = FALSE)
    plot_cells(cds_sub_DucT2, color_cells_by = "cell_cycle",cell_size=2,
               label_cell_groups=FALSE, show_trajectory_graph = FALSE) + scale_color_manual(values = colors_cc)

    cds_sub_DucT2_NewK <- cluster_cells(cds_sub_DucT2,k = k_cds_sub_DucT2, resolution=1e-5)
    plot_cells(cds_sub_DucT2_NewK, color_cells_by = "cluster",cell_size=2, label_cell_groups=FALSE, show_trajectory_graph = FALSE)

        ################ Cell-Cycle Scoring and Regression (DucT2) ################
        ## Convert Monocle3 Object to Seurat Object    # getFilePath("Monocle3_To_Seurat.R")
        marrow_sub_DucT2_NewK <- Monocle3_To_Seurat(cds_sub_DucT2_NewK,"sub_DT2") #sub_DT2TOP2ACTR:sub_DucT2_TOP2ACenter

        ###### Insert the cell cycle results from Monocle3 cds_sub into the  Seurat  marrow_sub ######
        marrow_sub_DucT2_NewK@active.ident <- cds_sub_DucT2_NewK@colData@listData$cell_cycle
        RidgePlot(marrow_sub_DucT2_NewK,cols = colors_cc, features = c(Main), ncol = 1)

        plot_cells(cds_sub_DucT2_NewK, genes=c("TOP2A"),cell_size=2, label_cell_groups=FALSE, show_trajectory_graph = FALSE)
        plot_cells(cds_sub_DucT2_NewK, color_cells_by="cell_cycle",cell_size=2, label_cell_groups=FALSE, show_trajectory_graph = FALSE) + scale_color_manual(values = colors_cc)

        ############    Cell discrimination by AddModuleScore (DucT2)   ############
        ## getFilePath("Monocle3_AddModuleScore.R")
        PDAC_Marker_file_Name <- c("GRUETZMANN_PANCREATIC_CANCER_UP")
        PDAC_Marker_Name <- c("PDAC_Marker")

        cds_sub_DucT2_NewK <- Monocle3_AddModuleScore(PDAC_Marker_file_Name,PDAC_Marker_Name,marrow_sub_DucT2_NewK,cds_sub_DucT2_NewK)
        plot_cells(cds_sub_DucT2_NewK, color_cells_by= PDAC_Marker_Name,
                   cell_size=2, label_cell_groups=FALSE, show_trajectory_graph = FALSE)


        ############    Find marker genes expressed by each cluster (DucT2)   ############
        marker_test_res_DucT2 <- top_markers(cds_sub_DucT2_NewK, group_cells_by="cluster")

        top_specific_markers_DucT2 <- marker_test_res_DucT2 %>% filter(fraction_expressing >= 0.10) %>%
                                                                group_by(cell_group) %>% top_n(10, pseudo_R2)
        top_specific_marker_ids_DucT2 <- unique(top_specific_markers_DucT2  %>% pull(gene_id))

        plot_genes_by_group(cds_sub_DucT2_NewK, top_specific_marker_ids_DucT2, group_cells_by="cluster",
                            ordering_type="maximal_on_diag", max.size=3)
        plot_genes_by_group(cds_sub_DucT2_NewK, top_specific_marker_ids_DucT2,group_cells_by="cluster",
                            ordering_type="cluster_row_col",max.size=3)

        ## Get marker genes from each cluster
        top_specific_markers_DucT2_Sub1 <- top_specific_markers_DucT2[top_specific_markers_DucT2$cell_group =="1",]
        top_specific_markers_DucT2_Sub2 <- top_specific_markers_DucT2[top_specific_markers_DucT2$cell_group =="2",]
        #...
        marker_test_res_DucT2_Sub1 <- marker_test_res_DucT2[marker_test_res_DucT2$cell_group =="1",]
        marker_test_res_DucT2_Sub2 <- marker_test_res_DucT2[marker_test_res_DucT2$cell_group =="2",]
        #...

        ## Export a marker genes information file
        write.table(marker_test_res_DucT2, file=paste0(PathName,"/",RVersion,"/",RVersion,"_",
                                                       "DucT2_marker_test_res.txt"),  sep="\t", row.names=FALSE)

        ## Generate a Garnett file
        # Require that markers have at least JS specificty score > 0.1 and be significant in the logistic test for identifying their cell type:
        garnett_markers_DucT2 <- marker_test_res_DucT2 %>% filter(marker_test_q_value < 0.01 & specificity >= 0.1) %>%
                                                           group_by(cell_group) %>% top_n(100, marker_score)

        # # Exclude genes that are good markers for more than one cell type:
        # garnett_markers_DucT2 <- garnett_markers_DucT2 %>%
        #   group_by(gene_short_name) %>%
        #   filter(n() == 1)

        generate_garnett_marker_file(garnett_markers_DucT2,max_genes_per_group = 100,
                                   file=paste0(PathName,"/",RVersion,"/",RVersion,"_","DucT2_marker_Garnett.txt"))

        
        # #************************************************************************************************************************#    
        ########################  DucT2_TOP2ACenter trajectories ##########################
        for (i in c(1:8)) {

            cds_sub_DucT2_TOP2ACenter_Tn <- choose_graph_segments(cds_sub_DucT2 ,clear_cds = FALSE)
            plot_cells(cds_sub_DucT2_TOP2ACenter_Tn, color_cells_by="cell_cycle",cell_size=2,
                       label_cell_groups=FALSE) + scale_color_manual(values = colors_cc)

            ###### Convert Monocle3 Object to Seurat Object ######
            # getFilePath("Monocle3_To_Seurat.R")
            marrow_sub_DucT2_TOP2ACenter_Tn <- Monocle3_To_Seurat(cds_sub_DucT2_TOP2ACenter_Tn,paste0("sub_DT2TOP2ACTR_T", i)) #sub_DT2TOP2ACTR:sub_DucT2_TOP2ACenter

            ###### Insert the cell cycle results from Monocle3 cds_sub into the  Seurat  marrow_sub ######
            marrow_sub_DucT2_TOP2ACenter_Tn@active.ident <- cds_sub_DucT2_TOP2ACenter_Tn@colData@listData$cell_cycle

            assign(paste0("marrow_sub_DucT2_TOP2ACenter_T", i),marrow_sub_DucT2_TOP2ACenter_Tn)

            plot_cells(cds_sub_DucT2_TOP2ACenter_Tn, color_cells_by="cell_cycle", label_cell_groups=FALSE) + scale_color_manual(values = colors_cc)
            assign(paste0("cds_sub_DucT2_TOP2ACenter_T", i),cds_sub_DucT2_TOP2ACenter_Tn)

            ###### PCA Scores for finding significantly different genes at the endpoints ######
            getFilePath("PCA_Threshold.R")
            PCA_DT2TOP2ACTR_Tn <- assign(paste0("PCA_DT2TOP2ACTR_T", i),marrow_sub_DucT2_TOP2ACenter_Tn@reductions[["pca"]]@feature.loadings)
            assign(paste0("PCA_DT2TOP2ACTR_T", i,"_PC_Sum"),PCA_Threshold_Pos(PCA_DT2TOP2ACTR_Tn, i ,PCAThreshold_Pos))
            assign(paste0("PCA_DT2TOP2ACTR_T", i,"_NC_Sum"),PCA_Threshold_Neg(PCA_DT2TOP2ACTR_Tn, i ,PCAThreshold_Neg))

            rm(cds_sub_DucT2_TOP2ACenter_Tn,marrow_sub_DucT2_TOP2ACenter_Tn)
          }


        # #************************************************************************************************************************#
        ########################  Heterogeneity center and Ori Ductal2 ##########################
            cds_sub_HeteroCent_OriDucT2 <- choose_cells(cds)
            #cds_subset <- reduce_dimension(cds_subset)
            plot_cells(cds_sub_HeteroCent_OriDucT2, label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 2)

            cds_sub_HeteroCent_KNew <- cluster_cells(cds_sub_HeteroCent_OriDucT2,k = 50, resolution=1e-5)
            plot_cells(cds_sub_HeteroCent_KNew, label_cell_groups=FALSE, color_cells_by = "cluster", show_trajectory_graph = FALSE,cell_size = 2)

                ##  Find marker genes expressed by each cluster
                marker_test_HeteroCent_OriDucT2 <- top_markers(cds_sub_HeteroCent_KNew, group_cells_by="cluster")

                top_specific_markers_HeteroCent_OriDucT2 <- marker_test_HeteroCent_OriDucT2 %>%
                                                            filter(fraction_expressing >= 0.10) %>%
                                                            group_by(cell_group) %>%
                                                            top_n(25, pseudo_R2)

                top_specific_marker_ids_HeteroCent_OriDucT2  <- unique(top_specific_markers_HeteroCent_OriDucT2  %>% pull(gene_id))


                plot_genes_by_group(cds_sub_HeteroCent_KNew,
                                    top_specific_marker_ids_HeteroCent_OriDucT2,
                                    group_cells_by="cluster",
                                    ordering_type="maximal_on_diag",
                                    max.size=3)

                top_specific_markers_HeteroCent_OriDucT2 <- data.frame(top_specific_markers_HeteroCent_OriDucT2)

                top_marker_HeteroCent_OriDucT2_Sub1 <- top_specific_markers_HeteroCent_OriDucT2[top_specific_markers_HeteroCent_OriDucT2$cell_group =="1",]
                top_marker_HeteroCent_OriDucT2_Sub2 <- top_specific_markers_HeteroCent_OriDucT2[top_specific_markers_HeteroCent_OriDucT2$cell_group =="2",]

