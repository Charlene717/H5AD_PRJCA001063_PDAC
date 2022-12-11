########################  cds_sub_Try  ##########################
plot_cells(cds, genes=c("TOP2A"),cell_size=1, label_cell_groups=FALSE, show_trajectory_graph = FALSE)
plot_cells(cds, color_cells_by="cell_cycle",cell_size=1, label_cell_groups=FALSE, show_trajectory_graph = FALSE) + scale_color_manual(values = colors_cc)


cds_sub_Try <- choose_cells(cds_subset_NewK)

cds_sub_Try_Ori <- cds_sub_Try
plot_cells(cds_sub_Try_Ori, label_cell_groups=FALSE)
plot_cells(cds_sub_Try_Ori, label_cell_groups=FALSE, show_trajectory_graph = FALSE, cell_size = 2)
plot_cells(cds_sub_Try_Ori, label_cell_groups=FALSE, show_trajectory_graph = FALSE, cell_size = 2, color_cells_by="cell_cycle")+ scale_color_manual(values = colors_cc)
plot_cells(cds_sub_Try_Ori, genes=c(Main), label_cell_groups=FALSE, show_trajectory_graph = FALSE, cell_size = 2)


        ###### Convert Monocle3 Object to Seurat Object ######
        # getFilePath("Monocle3_To_Seurat.R")
        marrow_sub_sub_Try <- Monocle3_To_Seurat(cds_sub_Try,"sub_DT2TOP2ACTR") #sub_DT2TOP2ACTR:sub_DucT2_TOP2ACenter
        
        ###### Assign Cell-Cycle Scores ######
        # getFilePath("Cell-Cycle Scoring and Regression.R")
        marrow_sub_sub_Try <- CCScorReg(GeneNAFMT,marrow_sub_sub_Try) #這個function存在於Cell-Cycle Scoring and Regression.R裡面
        # view cell cycle scores and phase assignments
        head(marrow_sub_sub_Try[[]])
        
        ## Plot the RidgePlot
        RidgePlot(marrow_sub_sub_Try,cols = colors_cc, features = c(Main), ncol = 1)
        RidgePlot(marrow_sub_sub_Try,cols = colors_cc, features = c(Main_Group), ncol = 2,y.max = 100) 
        RidgePlot(marrow_sub_sub_Try,cols = colors_cc, features = c(Main_Group), ncol = 2,log=TRUE)
        
        ###### Insert the cell cycle results from Seurat into the  Monocle3 cds object ######
        cds_sub_Try@colData@listData$cell_cycle <- marrow_sub_sub_Try@active.ident
        
        # ENSG00000131747 TOP2A
        plot_cells(cds_sub_Try, genes=c("TOP2A"),cell_size=2, label_cell_groups=FALSE, show_trajectory_graph = FALSE)
        plot_cells(cds_sub_Try, color_cells_by="cell_cycle",cell_size=2, label_cell_groups=FALSE, show_trajectory_graph = FALSE) + scale_color_manual(values = colors_cc)
        plot_cells(cds_sub_Try_Ori, color_cells_by="cell_cycle",cell_size=2, label_cell_groups=FALSE, show_trajectory_graph = FALSE) + scale_color_manual(values = colors_cc)
        
        ## Plot the Violin Plot 
        cds_sub_Try_Maingroup <- cds_sub_Try[rowData(cds_sub_Try)$gene_short_name %in% Main_Group,]
        plot_genes_violin(cds_sub_Try_Maingroup, group_cells_by="cell_cycle", ncol=2, log_scale = FALSE) +
          scale_fill_manual(values = colors_cc) + 
          theme(axis.text.x=element_text(angle=45, hjust=1))
        plot_genes_violin(cds_sub_Try_Maingroup, group_cells_by="cell_cycle", ncol=2, log_scale = TRUE) +
          scale_fill_manual(values = colors_cc) + 
          theme(axis.text.x=element_text(angle=45, hjust=1))
        
        ## Plot the pseudotime
        MainGroup_lineage_sub_Try <- cds_sub_Try[rowData(cds_sub_Try)$gene_short_name %in% Main_Group]
        plot_genes_in_pseudotime(MainGroup_lineage_sub_Try ,
                                 color_cells_by="cell_cycle",cell_size=2,
                                 min_expr=0.5)+ scale_color_manual(values = colors_cc)
        
        
        ##############  Annotate your cells according to type (Custom Marker)  ##############
        
        ############    Cell discrimination by AddModuleScore    ############
        getFilePath("Monocle3_AddModuleScore.R")
        
        PDAC_Marker_file_Name <- c("ALONSO_METASTASIS_EMT_UP")
        PDAC_Marker_Name <- c("EMT_Marker")
        
        cds_sub_Try <- Monocle3_AddModuleScore(PDAC_Marker_file_Name,PDAC_Marker_Name,marrow_sub_sub_Try,cds_sub_Try)
        plot_cells(cds_sub_Try, color_cells_by= PDAC_Marker_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE)
        
        # cds_sub_DucT2 <- Monocle3_AddModuleScore(PDAC_Marker_file_Name,PDAC_Marker_Name,marrow_sub_DucT2,cds_sub_DucT2)
        # plot_cells(cds_sub_DucT2, color_cells_by= Marker_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE)
        
        
        