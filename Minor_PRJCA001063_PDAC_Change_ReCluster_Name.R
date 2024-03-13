# colData(cds_sub_AcinaDucT_NewK_ReCluster)$assigned_cell_type <- 
#   as.character(clusters(cds_sub_AcinaDucT_NewK_ReCluster)[colnames(cds_sub_AcinaDucT_NewK_ReCluster)])
colData(cds_sub_AcinaDucT_NewK_ReCluster)$assigned_cell_type <- 
  dplyr::recode(colData(cds_sub_AcinaDucT_NewK_ReCluster)$assigned_cell_type,
                "AC"="AC",
                
                "nAtD"="nAtD",
                "nAtD"="nAtD",
                
                "aAtD"="aAtD",
                
                "ND01"="ND01",
                "ND02"="ND02",
                "ND03"="ND03",
                "ND04"="ND04",
                
                "AD"="AD",
                
                "MD00"="CDOri",
                
                "MDC01"="CoreCD01",
                "MDC02"="CoreCD02",
                "MDC03" ="CoreCD03",
                "MDC00"="CoreCD00",
                "MDC04"="CoreCD04",
                "MDC05"="CoreCD05",
                "MDC06"="CoreCD06",
                "MDC07"="CoreCD07",
                "MDC08"="CoreCD08",
                
                "MDO01"="DistalCD01",
                "MDO02"="DistalCD02",
                "MDO02"="DistalCD02",
                "MDO03"="DistalCD03",
                "MDO04"="DistalCD04",
                "MDO05"="DistalCD05",
                "MDO06"="DistalCD06",
                "MDO07"="DistalCD07",
                "MDO07"="DistalCD07",
                "MDO08"="DistalCD08",
                "MDO09"="DistalCD09",
                "MDO10"="DistalCD10",
                "MDO11"="DistalCD11")


cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["ReCluster"]] <- cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["assigned_cell_type"]]
plot_cells(cds_sub_AcinaDucT_NewK_ReCluster, color_cells_by = "ReCluster",cell_size=2, label_cell_groups=FALSE, show_trajectory_graph = FALSE)
