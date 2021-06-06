## (scRNA-seq data analysis for H5AD files)

#############
rm(list = ls()) # Clean variable

memory.limit(150000)
set.seed(1) # Fix the seed
############# Library list #############
library(SummarizedExperiment)
library(Seurat)
library(SeuratDisk)
library(stringr)
library(SeuratWrappers)
library(monocle3)
library(AnnotationDbi)
library(org.Mm.eg.db)
library('org.Hs.eg.db')
library(Hmisc)
library(tidyverse)
library(garnett)
# library(cicero)
# detach("package:monocle", unload = TRUE)
# 錯誤: package 'monocle' is required by 'cicero' so will not be detached

############# Import files settings #############
  ## General setting
  PathName = setwd(getwd())
  RVersion = "20210606V1"
  dir.create(paste0(PathName,"/",RVersion))
  
  ## Marker genes file
  Marker_file_Name <- c("NAKAMURA_METASTASIS_MODEL_UP")
  Marker_file <- paste0(PathName,"/",Marker_file_Name,".txt")
  Marker_List <- read.delim(Marker_file,header=F,sep= c("\t"))
  Marker_List2 <- as.data.frame(Marker_List[-1:-2,])
  
  Garnett_Marker_file_Name <- c("NAKAMURA_METASTASIS_MODEL_M18483")
  Garnett_Marker_file <- paste0(PathName,"/marker_file_",Garnett_Marker_file_Name,".txt")
  
  ## Cell cycle genes file
  cc.genes_list <- read.csv(paste0(PathName,"/Cell cycle/regev_lab_cell_cycle_genesCh.csv")) # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat. 
  # cc: Cell-Cycle

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


############# Parameter setting #############
  ## Gene list of interest 
  Main = c("TOP2A")
  Main_Group = c("TOP2A","TP53","CGAS","PTK2")
  Main_Group2 = c("KRAS","EXO1","TP53","CGAS","H2AX","PTK2")
  
  Main_Group3 = c("KRAS","EXO1","NSUN2","MUC1","AMBP","FXYD2","TOP2B","CCNE1")
  EMT_Meta = c("ANLN","APLP2","CD63","CDH2","CLIC4","CTSB","CX3CR1","DSG2","EDNRB")
  candidates14 = c("BRIP1","KIF23","TOP2A","FOSL1","FAM25A","ANLN","NCAPH","KRT9","MCM4","CKAP2L","CENPE","RACGAP1","DTL","RAD51AP1")
  
  DREAM_complex= c("RBL2","E2F4","E2F5","TFDP1","TFDP2")
  Regulators= c("TP53","YBX1","E2F1")
  
  ## Color setting
  colors_cc <- c("#59c26b", "#2e6087", "#417034") ## Color for Cell-Cycle
  colors_cc2 <- c("#FF59c26b", "#FF2e6087", "#FF417034") ## Color for Cell-Cycle
 # colors_cc <- c("#FF9912B3", "#32CD3299", "#4169E1B3") ## Color for Cell-Cycle
  
  ## Format of data
  GeneNAFMT <- c("HuGSymbol") # Gene names format of data: HuGSymbol,MouGSymbol,HuENSEMBL,MouENSEMBL
  
  ## Cluster cells setting
  k_cds_sub_DucT2 <- c(5) # k for ductal cell type2
  k_cds_sub_AcinaDucT <- 7 # k for Acinar + ductal cell type

  k_cds_sub_DucT2_HG <- c(4)
  
  ## Threshold of PCA scores
  PCAThreshold_Pos <- 0.03
  PCAThreshold_Neg <- -0.03

#********************************************************************************************************************#
  
##################  Function setting ################## 
  
  ## Call function
  filePath <- ""
  #匯入 同一個資料夾中的R檔案
  getFilePath <- function(fileName) {
    # path <- setwd("~")  #專案資料夾絕對路徑
    path <- setwd(getwd()) 
    #字串合併無間隔
    # 「<<-」為全域變數給值的指派
    filePath <<- paste0(path ,"/" , fileName)  
    # 載入檔案
    sourceObj <- source(filePath)
    return(sourceObj)
  }

#********************************************************************************************************************#
  
############# Import raw data #############
library(SummarizedExperiment)
library(Seurat)
library(SeuratDisk)

## Convert h5ad to h5seurat
Convert(paste0(PathName,"/StdWf1_PRJCA001063_CRC_besca2.raw.h5ad"), "PRJCA001063.h5seurat", assay = "RNA",) # This creates a copy of this .h5ad object reformatted into .h5seurat inside the example_dir directory
seuratObject <- LoadH5Seurat(paste0(PathName,"/PRJCA001063.h5seurat")) # This .d5seurat object can then be read in manually
 
## Convert Seurat Object to Monocle3 Object
library(SeuratWrappers)
cds <- as.cell_data_set(seuratObject) # Convert objects to Monocle3 'cell_data_set' objects

#---------------------------------------------------------------------------------------------------------------------#

    ############# Run Monocle3 #############
    library(monocle3)
    cds <- estimate_size_factors(cds) # issues with cds object in monocle3 #54 # https://github.com/satijalab/seurat-wrappers/issues/54
    
    ###### Pre-process the data ######
    cds <- preprocess_cds(cds, num_dim = 100)
    plot_pc_variance_explained(cds)

    ###### Reduce dimensionality and visualize the cells ######
    
    ##### UMAP #####
    cds <- reduce_dimension(cds,preprocess_method = 'PCA')
    plot_cells(cds)
    
    ##### (UMAP) Plot genes #####
    plot_cells(cds, genes=c("TOP2A","TOP2B","TP53")) #error
    ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cds@rowRanges@elementMetadata@listData$gene_short_name <- cds@assays@data@listData[["counts"]]@Dimnames[[1]]
    plot_cells(cds, genes=c("TOP2A","TOP2B","TP53","CCNE1")) #ok
    
    plot_cells(cds, genes=c(Main),cell_size=1,label_cell_groups = FALSE, show_trajectory_graph = FALSE)
    plot_cells(cds, genes=c(Main_Group),cell_size=0.5,label_cell_groups = FALSE, show_trajectory_graph = FALSE)
    plot_cells(cds, genes=c(Main_Group2),cell_size=0.5,label_cell_groups = FALSE, show_trajectory_graph = FALSE)

    ##### (UMAP) Plot different phenotype #####
    plot_cells(cds, color_cells_by="Cell_type", label_cell_groups=FALSE, show_trajectory_graph = FALSE)
    plot_cells(cds, color_cells_by="Type", label_cell_groups=FALSE, show_trajectory_graph = FALSE)
    plot_cells(cds, color_cells_by="Patient", label_cell_groups=FALSE, show_trajectory_graph = FALSE)
    plot_cells(cds, color_cells_by="CONDITION", label_cell_groups=FALSE, show_trajectory_graph = FALSE)
    
    
    ##### Group cells into clusters ######
    set.seed(1) # Fix the seed
    cds <- cluster_cells(cds)
    plot_cells(cds, color_cells_by = "partition", label_cell_groups=FALSE, show_trajectory_graph = FALSE)
    plot_cells(cds, color_cells_by = "cluster", label_cell_groups=FALSE, show_trajectory_graph = FALSE)
    plot_cells(cds, label_cell_groups=FALSE, show_trajectory_graph = FALSE)
    
    ############ Cell-Cycle Scoring and Regression - Monocle3 & Seurat Mutual conversion #############
    
    ############# Run Seurat #############
        ## Load package
        library(Seurat)
        library(SummarizedExperiment) 
        library(AnnotationDbi)
        library(org.Mm.eg.db)
        library('org.Hs.eg.db')
        library(Hmisc)
    
        ###### Convert Monocle3 Object to Seurat Object ######
        getFilePath("Monocle3_To_Seurat.R")
        marrow <- Monocle3_To_Seurat(cds,"cds") #這個function存在於Monocle3_To_Seurat.R裡面
        
        ###### Assign Cell-Cycle Scores ######
        getFilePath("Cell-Cycle Scoring and Regression.R")
        marrow <- CCScorReg(GeneNAFMT,marrow) #這個function存在於Cell-Cycle Scoring and Regression.R裡面
        
        RidgePlot(marrow,cols = colors_cc, features = c(Main), ncol = 1)
        RidgePlot(marrow,cols = colors_cc, features = c(Main_Group), ncol = 2,y.max = 100) 
        RidgePlot(marrow,cols = colors_cc, features = c(Main_Group), ncol = 2,log=TRUE) 
        
        ###### Insert the cell cycle results from Seurat into the  Monocle3 cds object ######
        cds@colData@listData$cell_cycle <- marrow@active.ident
        # cds@colData@listData$cell_cycle <- marrow@meta.data[["Phase"]]
        
        plot_cells(cds, color_cells_by="cell_cycle", label_cell_groups=FALSE ,show_trajectory_graph = F) + scale_color_manual(values = colors_cc)
        
        ## Plot the violin diagram
        Maingroup_ciliated_genes <- c(Main_Group)
        cds_marrow_cc <- cds[rowData(cds)$gene_short_name %in% Maingroup_ciliated_genes,]
        
        plot_genes_violin(cds_marrow_cc, group_cells_by="cell_cycle", ncol=2, log_scale = FALSE)+ scale_fill_manual(values = colors_cc)
        plot_genes_violin(cds_marrow_cc, group_cells_by="cell_cycle", ncol=2, log_scale = T)+ scale_fill_manual(values = colors_cc)
        plot_genes_violin(cds_marrow_cc, group_cells_by="cell_cycle", ncol=2, log_scale = T)+ scale_fill_manual(values = colors_cc)+
                         geom_boxplot(width=0.1, fill="white") + theme(axis.text.x=element_text(angle=45, hjust=1))
        

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
        
        # Marker_G0G1_file_Name <- c("REACTOME_G0_AND_EARLY_G1")
        # Marker_G0G1_Name <- c("G0G1")
        # cds <- Monocle3_AddModuleScore(Marker_G0G1_file_Name,Marker_G0G1_Name,marrow,cds)
        # plot_cells(cds, color_cells_by= Marker_G0G1_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1.2) +
        #   scale_colour_gradient2(low = "darkblue", mid = "#f7c211", high = "green", 
        #                          guide = "colourbar",midpoint = 0.15, labs(fill = Marker_G0G1_Name))
         
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

    
    ####################### Constructing single-cell trajectories #######################
    cds2 <-cds
    cds <- learn_graph(cds)
    plot_cells(cds,
               color_cells_by = "cluster",
               label_cell_groups=FALSE, label_leaves=FALSE, label_branch_points=FALSE, graph_label_size=1.5)
    
    cds <- order_cells(cds)
    plot_cells(cds,
               color_cells_by = "pseudotime",
               label_cell_groups=FALSE, label_leaves=FALSE, label_branch_points=FALSE, graph_label_size=1.5)
    
    MainGroup_lineage_cds <- cds[rowData(cds)$gene_short_name %in% Main_Group]
    
    plot_genes_in_pseudotime(MainGroup_lineage_cds, color_cells_by="cell_cycle",cell_size=2,
                             min_expr=0.5)+ scale_color_manual(values = colors_cc)
    
#************************************************************************************************************************#    
#************************************************************************************************************************#    
    
        ######################################  cds_subset ########################################
            ##################  Grab specific terms ################## 
                ## grepl Ductal cell type 2
                cds_sub_DucT2 <- cds[,grepl("Ductal cell type 2", colData(cds)$Cell_type, ignore.case=TRUE)]
                plot_cells(cds_sub_DucT2, color_cells_by="partition")
                plot_cells(cds_sub_DucT2, color_cells_by="partition", show_trajectory_graph = F)
                plot_cells(cds_sub_DucT2, genes=c(Main),cell_size=1,label_cell_groups = FALSE, show_trajectory_graph = FALSE)
                plot_cells(cds_sub_DucT2, genes=c(Main_Group),cell_size=0.5,label_cell_groups = FALSE, show_trajectory_graph = FALSE)
                
                ## Ductal cell type & cds_sub_AcinaDucT
                cds_sub_AcinaDucT <- cds[,colData(cds)$Cell_type %in% c("Acinar cell","Ductal cell type 1","Ductal cell type 2")]
                plot_cells(cds_sub_AcinaDucT, color_cells_by="cluster", show_trajectory_graph = F)
        
                
        ########################  AcinaDucT (choose_cells) ##########################
                cds_sub_AcinaDucT <- choose_cells(cds)
                plot_cells(cds_sub_AcinaDucT, color_cells_by="cluster", show_trajectory_graph = F)
                
                plot_cells(cds_sub_AcinaDucT , genes=c(Main_Group), show_trajectory_graph = FALSE) #ok
                plot_cells(cds_sub_AcinaDucT , genes=c("H2AX"), cell_size = 0.8,show_trajectory_graph = FALSE,label_cell_groups=FALSE) #ok
                plot_cells(cds_sub_AcinaDucT , genes=c("EXO1"), cell_size = 0.8,show_trajectory_graph = FALSE,label_cell_groups=FALSE) #ok
                plot_cells(cds_sub_AcinaDucT , genes=c("MYC"), cell_size = 0.8,show_trajectory_graph = FALSE,label_cell_groups=FALSE) #ok
                
                
                set.seed(1) # Fix the seed
                cds_sub_AcinaDucT_NewK <- cluster_cells(cds_sub_AcinaDucT,k = k_cds_sub_AcinaDucT, resolution=1e-5)
                plot_cells(cds_sub_AcinaDucT_NewK, color_cells_by = "cluster",cell_size=2, label_cell_groups=FALSE, show_trajectory_graph = FALSE)
                plot_cells(cds_sub_AcinaDucT_NewK, color_cells_by = "cluster",cell_size=2, 
                           label_cell_groups=TRUE, show_trajectory_graph = FALSE, group_label_size = 5)
                
                ######## Reorganize the Cluster for AcinaDucT ####
                cds_sub_AcinaDucT_NewK_ReCluster <- cds_sub_AcinaDucT_NewK
                colData(cds_sub_AcinaDucT_NewK_ReCluster)$assigned_cell_type <- 
                        as.character(clusters(cds_sub_AcinaDucT_NewK_ReCluster)[colnames(cds_sub_AcinaDucT_NewK_ReCluster)])
                colData(cds_sub_AcinaDucT_NewK_ReCluster)$assigned_cell_type <- 
                        dplyr::recode(colData(cds_sub_AcinaDucT_NewK_ReCluster)$assigned_cell_type,
                                                                        "6"="AC",
                                      
                                                                        "24"="nAtD",
                                                                        "29"="nAtD",
                                      
                                                                        "14"="aAtD",
                                      
                                                                        "8"="ND01",
                                                                        "1"="ND02",
                                                                        "3"="ND03",
                                                                        "13"="ND04",
                                      
                                                                        "2"="AD",
                                      
                                                                        "28"="MD00",
                                      
                                                                        "18"="MDC01",
                                                                        "26"="MDC02",
                                                                     #  "?" ="MDC03",
                                                                        "4"="MDC00",
                                                                        "19"="MDC04",
                                                                        "9"="MDC05",
                                                                        "16"="MDC06",
                                                                        "15"="MDC07",
                                                                     #  "?" ="MDC08",
                                      
                                                                        "11"="MDO01",
                                                                        "25"="MDO02",
                                                                        "30"="MDO02",
                                                                        "22"="MDO03",
                                                                        "12"="MDO04",
                                                                        "21"="MDO05",
                                                                        "10"="MDO06",
                                                                        "17"="MDO07",
                                                                        "27"="MDO07",
                                                                        "7"="MDO08",
                                                                        "23"="MDO09",
                                                                        "5"="MDO10",
                                                                        "20"="MDO11")
                
                
                cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["ReCluster"]] <- cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["assigned_cell_type"]]
                
                ## MDC03
                cds_sub_AcinaDucT_NewK_ReCluster_MDC03 <- choose_cells(cds_sub_AcinaDucT_NewK_ReCluster)
                colData(cds_sub_AcinaDucT_NewK_ReCluster_MDC03)$assigned_cell_type <- "MDC03"
                cds_sub_AcinaDucT_NewK_ReCluster_MDC03@colData@listData[["ReCluster"]] <- cds_sub_AcinaDucT_NewK_ReCluster_MDC03@colData@listData[["assigned_cell_type"]]
                
                ## MDC08
                cds_sub_AcinaDucT_NewK_ReCluster_MDC08 <- choose_cells(cds_sub_AcinaDucT_NewK_ReCluster)
                colData(cds_sub_AcinaDucT_NewK_ReCluster_MDC08)$assigned_cell_type <- "MDC08"
                cds_sub_AcinaDucT_NewK_ReCluster_MDC08@colData@listData[["ReCluster"]] <- cds_sub_AcinaDucT_NewK_ReCluster_MDC08@colData@listData[["assigned_cell_type"]]
                
                
                colData(cds_sub_AcinaDucT_NewK_ReCluster)[colnames(cds_sub_AcinaDucT_NewK_ReCluster_MDC03),]$assigned_cell_type <- colData(cds_sub_AcinaDucT_NewK_ReCluster_MDC03)$assigned_cell_type
                colData(cds_sub_AcinaDucT_NewK_ReCluster)[colnames(cds_sub_AcinaDucT_NewK_ReCluster_MDC08),]$assigned_cell_type <- colData(cds_sub_AcinaDucT_NewK_ReCluster_MDC08)$assigned_cell_type
                
                cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["ReCluster"]] <- cds_sub_AcinaDucT_NewK_ReCluster@colData@listData[["assigned_cell_type"]]
                plot_cells(cds_sub_AcinaDucT_NewK_ReCluster, color_cells_by = "ReCluster",cell_size=2, label_cell_groups=FALSE, show_trajectory_graph = FALSE)
                plot_cells(cds_sub_AcinaDucT_NewK_ReCluster, color_cells_by = "ReCluster",cell_size=2, 
                           label_cell_groups=TRUE, show_trajectory_graph = FALSE, group_label_size =4)
                
                
                
                ############    Find marker genes expressed by each cluster (AcinaDucT)   ############
                set.seed(1) # Fix the seed
                marker_test_res_AcinaDucT <- top_markers(cds_sub_AcinaDucT_NewK_ReCluster,
                                                         genes_to_test_per_group = 50, group_cells_by="ReCluster")
                
                top_specific_markers_AcinaDucT <- marker_test_res_AcinaDucT %>% filter(fraction_expressing >= 0.10) %>%
                                                  group_by(cell_group) %>% top_n(10, pseudo_R2)
                top_specific_marker_ids_AcinaDucT <- unique(top_specific_markers_AcinaDucT  %>% pull(gene_id))
                
                plot_genes_by_group(cds_sub_AcinaDucT_NewK_ReCluster, top_specific_marker_ids_AcinaDucT, group_cells_by="cluster",
                                    ordering_type="maximal_on_diag", max.size=3)
                plot_genes_by_group(cds_sub_AcinaDucT_NewK_ReCluster, top_specific_marker_ids_AcinaDucT,group_cells_by="cluster",
                                    ordering_type="cluster_row_col",max.size=3)
                
                ## Get marker genes from each cluster
                top_specific_markers_AcinaDucT_SubAD <- top_specific_markers_AcinaDucT[top_specific_markers_AcinaDucT$cell_group =="AD",]
                #top_specific_markers_AcinaDucT_Sub2 <- top_specific_markers_AcinaDucT[top_specific_markers_AcinaDucT$cell_group =="2",]
                #...
                marker_test_res_AcinaDucT_SubAD <- marker_test_res_AcinaDucT[marker_test_res_AcinaDucT$cell_group =="AD",]
                #marker_test_res_AcinaDucT_Sub2 <- marker_test_res_AcinaDucT[marker_test_res_AcinaDucT$cell_group =="2",]
                #...
                
                ## Export a marker genes information file
                write.table(marker_test_res_AcinaDucT, file=paste0(PathName,"/",RVersion,"/",RVersion,"_", 
                                                               "AcinaDucT_marker_test_res_GPG50.txt"),  sep="\t", row.names=FALSE)
                
                ## Generate a Garnett file
                # Require that markers have at least JS specificty score > 0.1 and be significant in the logistic test for identifying their cell type:
                garnett_markers_AcinaDucT <- marker_test_res_AcinaDucT %>% filter(marker_test_q_value < 0.05 & specificity >= 0.1) %>%
                                             group_by(cell_group) %>% top_n(100, marker_score)
                
                # # Exclude genes that are good markers for more than one cell type:
                # garnett_markers_DucT2 <- garnett_markers_DucT2 %>% 
                #   group_by(gene_short_name) %>%
                #   filter(n() == 1)
                
                generate_garnett_marker_file(garnett_markers_AcinaDucT,max_genes_per_group = 100, 
                                             file=paste0(PathName,"/",RVersion,"/",RVersion,"_","AcinaDucT_marker_Garnett_GPG50_q005spe01.txt"))
                

                ################ Plot Cell-Cycle Scoring and Regression (AcinaDucT) ################ 
                ## Convert Monocle3 Object to Seurat Object    # getFilePath("Monocle3_To_Seurat.R")
                marrow_sub_AcinaDucT_NewK_ReCluster <- Monocle3_To_Seurat(cds_sub_AcinaDucT_NewK_ReCluster,"sub_AcinaDucT") #sub_DT2TOP2ACTR:sub_DucT2_TOP2ACenter
                
                  
                ###### Insert the cell cycle results from Monocle3 cds_sub into the  Seurat  marrow_sub ######
                marrow_sub_AcinaDucT_NewK_ReCluster@active.ident <- cds_sub_AcinaDucT_NewK_ReCluster@colData@listData$cell_cycle
                RidgePlot(marrow_sub_AcinaDucT_NewK_ReCluster,cols = colors_cc, features = c(Main_Group), ncol = 2)
                
                plot_cells(cds_sub_AcinaDucT_NewK_ReCluster, genes=c("TOP2A"),cell_size=2, label_cell_groups=FALSE, show_trajectory_graph = FALSE)
                plot_cells(cds_sub_AcinaDucT_NewK_ReCluster, color_cells_by="cell_cycle",cell_size=2, label_cell_groups=FALSE, show_trajectory_graph = FALSE) + scale_color_manual(values = colors_cc)
                
                ## Plot the violin diagram
                Maingroup_ciliated_genes <- c(Main_Group)
                cds_marrow_cc_AcinaDucT_NewK_ReCluster <- cds_sub_AcinaDucT_NewK_ReCluster[rowData(cds_sub_AcinaDucT_NewK_ReCluster)$gene_short_name %in% Maingroup_ciliated_genes,]
                
                plot_genes_violin(cds_marrow_cc_AcinaDucT_NewK_ReCluster, group_cells_by="cell_cycle", ncol=2, log_scale = F)+ scale_fill_manual(values = colors_cc)
                plot_genes_violin(cds_marrow_cc_AcinaDucT_NewK_ReCluster, group_cells_by="cell_cycle", ncol=2, log_scale = T)+ scale_fill_manual(values = colors_cc)
                plot_genes_violin(cds_marrow_cc, group_cells_by="cell_cycle", ncol=2, log_scale = T)+ scale_fill_manual(values = colors_cc)+
                  geom_boxplot(width=0.1, fill="white", alpha = 0.7) + theme(axis.text.x=element_text(angle=45, hjust=1))
                
                
                ############    Plot Cell discrimination by AddModuleScore (AcinaDucT)   ############
                plot_cells(cds_sub_AcinaDucT_NewK_ReCluster, color_cells_by= Marker_PDAC_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1.2) +
                           scale_colour_gradient2(low = "#440075", mid = "#ffd261", high = "#4aff8c", 
                                                 guide = "colourbar",midpoint = 0.2, labs(fill = Marker_PDAC_Name))
                
                plot_cells(cds_sub_AcinaDucT_NewK_ReCluster, color_cells_by= Marker_EMT_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1.2) +
                  scale_colour_gradient2(low = "#440075", mid = "#ffd261", high = "#4aff8c", 
                                         guide = "colourbar",midpoint = 0.2, labs(fill = Marker_EMT_Name))
                
                plot_cells(cds_sub_AcinaDucT_NewK_ReCluster, color_cells_by= Marker_ChroSt_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1.2) +
                  scale_colour_gradient2(low = "darkblue", mid = "#f7c211", high = "green", 
                                         guide = "colourbar",midpoint = 0.15, labs(fill = Marker_ChroSt_Name))
                
                plot_cells(cds_sub_AcinaDucT_NewK_ReCluster, color_cells_by= Marker_Mig_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1.2) +
                  scale_colour_gradient2(low = "darkblue", mid = "#f7c211", high = "green", 
                                         guide = "colourbar",midpoint = 0.15, labs(fill = Marker_Mig_Name))
                
                plot_cells(cds_sub_AcinaDucT_NewK_ReCluster, color_cells_by= Marker_Meta_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1.2) +
                  scale_colour_gradient2(low = "darkblue", mid = "#f7c211", high = "green", 
                                         guide = "colourbar",midpoint = 0.15, labs(fill = Marker_Meta_Name))
        
                plot_cells(cds_sub_AcinaDucT_NewK_ReCluster, color_cells_by= Marker_NE_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1.2) +
                  scale_colour_gradient2(low = "darkblue", mid = "#f7c211", high = "green", 
                                         guide = "colourbar",midpoint = 0.15, labs(fill = Marker_NE_Name))
                
                plot_cells(cds_sub_AcinaDucT_NewK_ReCluster, color_cells_by= Marker_NP_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1.2) +
                  scale_colour_gradient2(low = "darkblue", mid = "#f7c211", high = "green", 
                                         guide = "colourbar",midpoint = 0.15, labs(fill = Marker_NP_Name))
                
                plot_cells(cds_sub_AcinaDucT_NewK_ReCluster, color_cells_by= Marker_ATR_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1.2) +
                  scale_colour_gradient2(low = "darkblue", mid = "#f7c211", high = "green", 
                                         guide = "colourbar",midpoint = 0.15, labs(fill = Marker_ATR_Name))
                
                plot_cells(cds_sub_AcinaDucT_NewK_ReCluster, color_cells_by= Marker_ACST_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1.2) +
                  scale_colour_gradient2(low = "darkblue", mid = "#f7c211", high = "green", 
                                         guide = "colourbar",midpoint = 0.15, labs(fill = Marker_ACST_Name))
                
               ####################### Constructing single-cell trajectories (AcinaDucT) #######################
               cds3 <- as.cell_data_set(seuratObject) # Convert objects to Monocle3 'cell_data_set' objects
               cds3 <- estimate_size_factors(cds3) # issues with cds object in monocle3 #54 # https://github.com/satijalab/seurat-wrappers/issues/54
               
               ## Pre-process the data 
               cds3 <- preprocess_cds(cds3, num_dim = 100)
               plot_pc_variance_explained(cds3)
               
               ## Reduce dimensionality and visualize the cells
               
               ## UMAP
               cds3 <- reduce_dimension(cds3,preprocess_method = 'PCA')
               plot_cells(cds3)
               
               cds_sub_AcinaDucT2 <- choose_cells(cds3)
               cds_sub_AcinaDucT2 <- cluster_cells(cds_sub_AcinaDucT2)
               cds_sub_AcinaDucT2 <- learn_graph(cds_sub_AcinaDucT2, use_partition = F)
               plot_cells(cds_sub_AcinaDucT2,
                          color_cells_by = "cluster",
                          label_cell_groups=FALSE, label_leaves=FALSE, label_branch_points=FALSE, graph_label_size=1.5)
               
               cds_sub_AcinaDucT2 <- order_cells(cds_sub_AcinaDucT2)
               plot_cells(cds_sub_AcinaDucT2,
                          color_cells_by = "pseudotime",
                          label_cell_groups=FALSE, label_leaves=FALSE, label_branch_points=FALSE, graph_label_size=1.5)
               
               # MainGroup_lineage_cds_sub_AcinaDucT2 <- cds_sub_AcinaDucT2[rowData(cds_sub_AcinaDucT2)$gene_short_name %in% Main_Group]
               # 
               # plot_genes_in_pseudotime(MainGroup_lineage_cds_sub_AcinaDucT2, color_cells_by="cell_cycle",cell_size=2,
               #                          min_expr=0.5)+ scale_color_manual(values = colors_cc)
               
               ##
               cds_sub_AcinaToDucT2 <- choose_cells(cds3)
               cds_sub_AcinaToDucT2 <- cluster_cells(cds_sub_AcinaToDucT2)
               cds_sub_AcinaToDucT2 <- learn_graph(cds_sub_AcinaToDucT2, use_partition = F)
               plot_cells(cds_sub_AcinaToDucT2,
                          color_cells_by = "cluster",
                          label_cell_groups=FALSE, label_leaves=FALSE, label_branch_points=FALSE, graph_label_size=1.5)
               
               ######   PCA for trajectories
               for (i in c(8:8)) {
                 
                 cds_sub_DucT2_TOP2ACenter_Tn <- choose_graph_segments(cds_sub_AcinaToDucT2 ,clear_cds = FALSE)
                 plot_cells(cds_sub_DucT2_TOP2ACenter_Tn, color_cells_by="cluster",cell_size=2, 
                            label_cell_groups=FALSE) + scale_color_manual(values = colors_cc)
                 
                 ###### Convert Monocle3 Object to Seurat Object ######
                 # getFilePath("Monocle3_To_Seurat.R")
                 marrow_sub_DucT2_TOP2ACenter_Tn <- Monocle3_To_Seurat(cds_sub_DucT2_TOP2ACenter_Tn,paste0("sub_DT2TOP2ACTR_T", i)) #sub_DT2TOP2ACTR:sub_DucT2_TOP2ACenter
                 
                 # ###### Insert the cell cycle results from Monocle3 cds_sub into the  Seurat  marrow_sub ######
                 # marrow_sub_DucT2_TOP2ACenter_Tn@active.ident <- cds_sub_DucT2_TOP2ACenter_Tn@colData@listData$cell_cycle
                 
                 assign(paste0("marrow_sub_DucT2_TOP2ACenter_T", i),marrow_sub_DucT2_TOP2ACenter_Tn)
                 
                 plot_cells(cds_sub_DucT2_TOP2ACenter_Tn, color_cells_by="cluster", label_cell_groups=FALSE) + scale_color_manual(values = colors_cc)
                 assign(paste0("cds_sub_DucT2_TOP2ACenter_T", i),cds_sub_DucT2_TOP2ACenter_Tn)
                 
                 ###### PCA Scores for finding significantly different genes at the endpoints ######
                 getFilePath("PCA_Threshold.R")
                 PCA_DT2TOP2ACTR_Tn <- assign(paste0("PCA_DT2TOP2ACTR_T", i),marrow_sub_DucT2_TOP2ACenter_Tn@reductions[["pca"]]@feature.loadings)
                 assign(paste0("PCA_DT2TOP2ACTR_T", i,"_PC_Sum"),PCA_Threshold_Pos(PCA_DT2TOP2ACTR_Tn, i ,PCAThreshold_Pos))
                 assign(paste0("PCA_DT2TOP2ACTR_T", i,"_NC_Sum"),PCA_Threshold_Neg(PCA_DT2TOP2ACTR_Tn, i ,PCAThreshold_Neg))
                 
                 rm(cds_sub_DucT2_TOP2ACenter_Tn,marrow_sub_DucT2_TOP2ACenter_Tn)
               }
               
               cds_sub_DucT2_TOP2ACenter_T4@rowRanges@elementMetadata@listData$gene_short_name <- cds_sub_DucT2_TOP2ACenter_T4@assays@data@listData[["counts"]]@Dimnames[[1]]
               
               plot_cells(cds_sub_DucT2_TOP2ACenter_T4, label_cell_groups=FALSE, show_trajectory_graph = T)
               
              
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
        
        
        
        ########################  DucT2_TOP2A_Center ##########################
            cds_sub_DucT2_TOP2ACenter <- choose_cells(cds_sub_AcinaDucT_NewK_ReCluster)
            plot_cells(cds_sub_DucT2_TOP2ACenter, genes=c("TOP2A"),cell_size=2, label_cell_groups=FALSE, show_trajectory_graph = FALSE)
            plot_cells(cds_sub_DucT2_TOP2ACenter, color_cells_by="cell_cycle",cell_size=2, label_cell_groups=FALSE, show_trajectory_graph = FALSE) + scale_color_manual(values = colors_cc)
  
            
                ###### Convert Monocle3 Object to Seurat Object ######
                # getFilePath("Monocle3_To_Seurat.R")
                marrow_sub_DucT2_TOP2ACenter <- Monocle3_To_Seurat(cds_sub_DucT2_TOP2ACenter,"sub_DT2TOP2ACTR") #sub_DT2TOP2ACTR:sub_DucT2_TOP2ACenter
                
                ###### Insert the cell cycle results from Monocle3 cds_sub into the  Seurat  marrow_sub ######
                marrow_sub_DucT2_TOP2ACenter@active.ident <- cds_sub_DucT2_TOP2ACenter@colData@listData$cell_cycle
        
                ## Plot the RidgePlot
                RidgePlot(marrow_sub_DucT2_TOP2ACenter,cols = colors_cc, features = c(Main), ncol = 1)
                RidgePlot(marrow_sub_DucT2_TOP2ACenter,cols = colors_cc, features = c(Main_Group), ncol = 2,log=TRUE) 
                RidgePlot(marrow_sub_DucT2_TOP2ACenter,cols = colors_cc, features = c(Main_Group), ncol = 2,y.max = 100) 
                
                
                ## Plot the Violin Plot 
                cds_sub_DT2TOP2ACTR_Maingroup <- cds_sub_DucT2_TOP2ACenter[rowData(cds_sub_DucT2_TOP2ACenter)$gene_short_name %in% Main_Group,]
                plot_genes_violin(cds_sub_DT2TOP2ACTR_Maingroup, group_cells_by="cell_cycle", ncol=2, log_scale = FALSE) +
                                  scale_fill_manual(values = colors_cc) + 
                                  theme(axis.text.x=element_text(angle=45, hjust=1))
      
                plot_genes_violin(cds_sub_DT2TOP2ACTR_Maingroup, group_cells_by="cell_cycle", ncol=2, log_scale = T) +
                  scale_fill_manual(values = colors_cc) + theme(axis.text.x=element_text(angle=45, hjust=1))
                plot_genes_violin(cds_sub_DT2TOP2ACTR_Maingroup, group_cells_by="cell_cycle", ncol=2, log_scale = T)+ scale_fill_manual(values = colors_cc)+
                  geom_boxplot(width=0.1, fill="white",alpha=(0.7)) + theme(axis.text.x=element_text(angle=45, hjust=1))
                
                ## Plot pseudotime
                MainGroup_lineage_sub_DT2TOP2ACTR <- cds_sub_DucT2_TOP2ACenter[rowData(cds_sub_DucT2_TOP2ACenter)$gene_short_name %in% Main_Group]
                plot_genes_in_pseudotime(MainGroup_lineage_sub_DT2TOP2ACTR,
                                         color_cells_by="cell_cycle",cell_size=2,
                                         min_expr=0.5)+ scale_color_manual(values = colors_cc)
            
                
            
        #************************************************************************************************************************#    
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
        
            
        #************************************************************************************************************************#    
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
        
        
        
        #************************************************************************************************************************#    
        ############(ERROR) Finding modules of co-regulated genes ############
        ciliated_cds_pr_test_res <- graph_test(cds_subTra2, neighbor_graph="principal_graph", cores=4)
        pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))
        gene_module_df <- find_gene_modules(cds_subTra2[pr_deg_ids,], resolution=1e-2)
        
        cell_group_df <- tibble::tibble(cell=row.names(colData(cds_subTra2)),
                                        cell_group=partitions(cds)[colnames(cds_subTra2)])
        agg_mat <- aggregate_gene_expression(cds_subTra2, gene_module_df, cell_group_df)
        row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
        colnames(agg_mat) <- stringr::str_c("Partition ", colnames(agg_mat))
        
        pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                           scale="column", clustering_method="ward.D2",
                           fontsize=6)


#************************************************************************************************************************#    
############ Working with 3D trajectories ############
    cds_3d <- reduce_dimension(cds, max_components = 3,preprocess_method = 'PCA')
    cds_3d <- cluster_cells(cds_3d)
    cds_3d <- learn_graph(cds_3d)
    # cds_3d <- order_cells(cds_3d, root_pr_nodes=get_earliest_principal_node(cds))
    # # Error in get_earliest_principal_node(cds) : 
    # #   沒有這個函數 "get_earliest_principal_node"
    #
    # cds_3d <- order_cells(cds_3d)
    
    cds_3d_plot_obj <- plot_cells_3d(cds_3d, color_cells_by="partition")
    
        plot_cells_3d(cds_3d)
        plot_cells_3d(cds_3d, color_cells_by="cluster", show_trajectory_graph = FALSE)
        plot_cells_3d(cds_3d, color_cells_by="Cell_type", show_trajectory_graph = FALSE)
        plot_cells_3d(cds_3d, color_cells_by="Type", show_trajectory_graph = FALSE)
        plot_cells_3d(cds_3d, color_cells_by="Patient", show_trajectory_graph = FALSE)
    
        plot_cells_3d(cds_3d, color_cells_by="PDAC_Marker", show_trajectory_graph = FALSE)
        
        plot_cells_3d(cds_3d, genes = Main, show_trajectory_graph = FALSE)
        plot_cells_3d(cds_3d, color_cells_by="cell_cycle", show_trajectory_graph = FALSE)
        
    cds_3d_sub_DucT2_TOP2ACenter <- reduce_dimension(cds_sub_DucT2_TOP2ACenter, max_components = 3,preprocess_method = 'PCA')
        plot_cells_3d(cds_3d_sub_DucT2_TOP2ACenter, genes = Main, show_trajectory_graph = FALSE)
        plot_cells_3d(cds_3d_sub_DucT2_TOP2ACenter, color_cells_by="cell_cycle", show_trajectory_graph = FALSE)
        
        cds_3d_sub_DucT2 <- reduce_dimension(cds_sub_DucT2, max_components = 3,preprocess_method = 'PCA')
        plot_cells_3d(cds_3d_sub_DucT2, genes = Main, show_trajectory_graph = FALSE)
        plot_cells_3d(cds_3d_sub_DucT2, color_cells_by="cell_cycle", show_trajectory_graph = FALSE)

 
