## (scRNA-seq data analysis for H5AD files)

#############
rm(list = ls()) # Clean variable

memory.limit(150000)

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
  RVersion = "20210525V1"
  dir.create(paste0(PathName,"/",RVersion))
  
############# Parameter setting #############
  ## Gene list of interest 
  Main = c("TOP2A")
  Main_Group = c("TOP2A","TOP2B","TP53","CCNE1")
  Main_Group2 = c("KRAS","EXO1","NSUN2","MUC1","AMBP","FXYD2")
  EMT_Meta = c("ANLN","APLP2","CD63","CDH2","CLIC4","CTSB","CX3CR1","DSG2","EDNRB")
  candidates14 = c("BRIP1","KIF23","TOP2A","FOSL1","FAM25A","ANLN","NCAPH","KRT9","MCM4","CKAP2L","CENPE","RACGAP1","DTL","RAD51AP1")
  
  DREAM_complex= c("RBL2","E2F4","E2F5","TFDP1","TFDP2")
  Regulators= c("TP53","YBX1","E2F1")
  
  ## Color setting
  colors_cc <- c("#FF9912B3", "#32CD3299", "#4169E1B3") ## Color for Cell-Cycle
  colors_cc2 <- c("#FF59c26b", "#FF2e6087", "#FF417034") ## Color for Cell-Cycle
  
    
########## Save the Figures ##########
    ##### Save as png #####
    png(paste0(PathName,"/",RVersion,"/",RVersion,"_","CC_Violin_Main.png")) # 設定輸出圖檔
    plot_genes_violin(cds_marrow_cc, group_cells_by="cell_cycle", ncol=2, log_scale = FALSE)+ 
      scale_fill_manual(values = colors_cc)+
      theme(axis.text.x=element_text(angle=45, hjust=1))
    dev.off() # 關閉輸出圖檔
    
    ##### Save as pdf #####
    pdf(paste0(PathName,"/",RVersion,"/",RVersion,"_","CC_Violin_Main.pdf")) # 設定輸出圖檔
    plot_genes_violin(cds_marrow_cc, group_cells_by="cell_cycle", ncol=2, log_scale = FALSE)+ 
      scale_fill_manual(values = colors_cc)+
      theme(axis.text.x=element_text(angle=45, hjust=1))
    dev.off() # 關閉輸出圖檔
    
    #************************************************************************************************************************#
    #*  ################ generating new theme  ################
    theme_uwv <- 
      function(base_size = 22,                                                                 #general font size
               base_family = ""){                                                              #general font type
        theme_hc(base_size = base_size,                                                      #basic theme to start from
                 base_family = base_family)   %+replace%  
          theme(axis.text.x = element_text(face="bold",  size=14),
                axis.text.y = element_text(face="bold",size=14),
                axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
                axis.title = element_text(size = rel(1.5),face="bold"),
                plot.title = element_text(color="black", size=20, face="bold.italic",hjust = 0.5,vjust =-1),
                legend.title = element_text(size=12, color = "black", face="bold"),
                legend.text = element_text(colour="black", size=12,face="bold"),
                aspect.ratio=1, #square plot
                complete = TRUE) 
      }
    
    ##!!!!!! https://github.com/tidyverse/ggplot2/issues/3474   
    
    #************************************************************************************************************************#    
    ##################  Phenotype ##################
    plot_cells(cds, color_cells_by= "Type", 
               label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1)  + 
#      scale_color_manual(values = colors_cc)+
      ggtitle("Type")+
      theme(axis.text.x = element_text(face="bold",  size=14),
            axis.text.y = element_text(face="bold",size=14),
            axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
            axis.title = element_text(size = rel(1.5),face="bold"),
            plot.title = element_text(color="black", size=20, 
                                      face="bold.italic",hjust = 0.1,vjust =-8), # margin = margin(t = 0.5, b = -7),
            legend.title = element_text(size=12, color = "black", face="bold"),
            legend.text = element_text(colour="black", size=12,face="bold"),
            legend.position = c(0.1, 0.18),
            aspect.ratio=1) + #square plot
      theme(axis.line.x = element_line(colour = "black", size = 0.8),
            axis.line.y = element_line(colour = "black", size = 0.8))

    plot_cells(cds, color_cells_by= "Cell_type", 
               label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1)  + 
      #      scale_color_manual(values = colors_cc)+
      ggtitle("Cell type")+
      theme(axis.text.x = element_text(face="bold",  size=14),
            axis.text.y = element_text(face="bold",size=14),
            axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
            axis.title = element_text(size = rel(1.5),face="bold"),
            plot.title = element_text(color="black", size=20, 
                                      face="bold.italic",hjust = 0.1,vjust =-8), # margin = margin(t = 0.5, b = -7),
            legend.title = element_text(size=12, color = "black", face="bold"),
            legend.text = element_text(colour="black", size=12,face="bold"),
          #  legend.position = c(0.1, 0.18),
            aspect.ratio=1) + #square plot
      theme(axis.line.x = element_line(colour = "black", size = 0.8),
            axis.line.y = element_line(colour = "black", size = 0.8))
    
        
    
    plot_cells(cds, color_cells_by= "Patient", 
               label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1)  + 
      #      scale_color_manual(values = colors_cc)+
      ggtitle("Patient")+
      theme(axis.text.x = element_text(face="bold",  size=14),
            axis.text.y = element_text(face="bold",size=14),
            axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
            axis.title = element_text(size = rel(1.5),face="bold"),
            plot.title = element_text(color="black", size=20, 
                                      face="bold.italic",hjust = 0.1,vjust =-8), # margin = margin(t = 0.5, b = -7),
            legend.title = element_text(size=12, color = "black", face="bold"),
            legend.text = element_text(colour="black", size=12,face="bold"),
          #  legend.position = c(0.1, 0.18),
            aspect.ratio=1) + #square plot
      theme(axis.line.x = element_line(colour = "black", size = 0.8),
            axis.line.y = element_line(colour = "black", size = 0.8))
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
    cds_sub_AcinaDucT <- cds[,colData(cds)$Cell_type %in% c("cds_sub_AcinaDucT","Ductal cell type 1","Ductal cell type 2")]
    plot_cells(cds_sub_AcinaDucT, color_cells_by="cluster", show_trajectory_graph = F)
    
    
         ########################  AcinaDucT (choose_cells) ##########################
         #######   Genes  ####### 
          cds_sub_AcinaDucT <- choose_cells(cds)
          plot_cells(cds_sub_AcinaDucT, color_cells_by="cluster", show_trajectory_graph = F)
          
          
          plot_cells(cds_sub_AcinaDucT, genes = Main_Group, show_trajectory_graph = F,label_cell_groups = F
                     ,cell_size =1 )+
            scale_colour_gradient2(low = "black", mid = "#3528c7", high = "#ff85e9", 
                                   guide = "colourbar",midpoint = 0, labs(fill = "Exp")) +
            ggtitle("Test Title")+
            theme(axis.text.x = element_text(face="bold",  size=14),
                  axis.text.y = element_text(face="bold",size=14),
                  axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
                  axis.title = element_text(size = rel(1.5),face="bold"),
                  plot.subtitle = element_text(color="black", size=20, face="bold.italic",hjust = 0.5,vjust =-1),
                  legend.title = element_text(size=12, color = "black", face="bold"),
                  legend.text = element_text(colour="black", size=12,face="bold"),
                  legend.background = element_rect(fill = alpha("white", 0.5)),
                  legend.position = c(0.5, 0.5),
                  legend.direction = "horizontal",
                  aspect.ratio=1) + #square plot
            theme(axis.line.x = element_line(colour = "black", size = 1),
                  axis.line.y = element_line(colour = "black", size = 1))
          
          ##!!!!!! 
          Main = "TOP2A"
          plot_cells(cds_Patient1, genes = Main, show_trajectory_graph = F,label_cell_groups = F
                     ,cell_size =1 )+
            scale_colour_gradient2(low = "black", mid = "#3528c7", high = "#ff85e9", 
                                   guide = "colourbar",midpoint = 0, labs(fill ="Exp")) +
            ggtitle(Main)+
            theme(axis.text.x = element_text(face="bold",  size=14),
                  axis.text.y = element_text(face="bold",size=14),
                  axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
                  axis.title = element_text(size = rel(1.5),face="bold"),
                  plot.title = element_text(color="black", size=17, 
                                            face="bold.italic",hjust = 0.1,vjust =-10), # margin = margin(t = 0.5, b = -7),
             #     plot.background = element_rect(fill = 'chartreuse'),
                  legend.title = element_text(size=12, color = "black", face="bold"),
                  legend.text = element_text(colour="black", size=12,face="bold"),
                  legend.background = element_rect(fill = alpha("white", 0.5)),
                  legend.position = c(0.1, 0.18),
             #     plot.text = element_text(size = 20),
                  aspect.ratio=1) + #square plot
            theme(axis.line.x = element_line(colour = "black", size = 0.8),
                  axis.line.y = element_line(colour = "black", size = 0.8))
 
          
          ##!!!!!! Sub
          Main = "TOP2A"
          Sub_Name = "Patient1"
          Sub = cds_Patient1
          plot_cells(Sub, genes = Main, show_trajectory_graph = F,label_cell_groups = F
                     ,cell_size =1 )+
            scale_colour_gradient2(low = "black", mid = "#3528c7", high = "#ff85e9", 
                                   guide = "colourbar",midpoint = 0, labs(fill ="Exp")) +
            #            ggtitle(Main)+
            ggtitle(paste0(Main,"(",Sub_Name,")"))+            
            theme(axis.text.x = element_text(face="bold",  size=14),
                  axis.text.y = element_text(face="bold",size=14),
                  axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
                  axis.title = element_text(size = rel(1.5),face="bold"),
                  plot.title = element_text(color="black", size=17, 
                                            face="bold.italic",hjust = 0.1,vjust =-10), # margin = margin(t = 0.5, b = -7),
                  #     plot.background = element_rect(fill = 'chartreuse'),
                  legend.title = element_text(size=12, color = "black", face="bold"),
                  legend.text = element_text(colour="black", size=12,face="bold"),
                  legend.background = element_rect(fill = alpha("white", 0.5)),
                  legend.position = c(0.1, 0.18),
                  #     plot.text = element_text(size = 20),
                  aspect.ratio=1) + #square plot
            theme(axis.line.x = element_line(colour = "black", size = 0.8),
                  axis.line.y = element_line(colour = "black", size = 0.8))
          
          
          #######   Cell cycle  #######   
          colors_cc <-c("#FF9912B3", "#2e6087", "#417034")  ## Color for Cell-Cycle
          ## AcinaDucT
          plot_cells(cds_sub_AcinaDucT_NewK_ReCluster, color_cells_by= "cell_cycle", 
                     label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1)  + 
            scale_color_manual(values = colors_cc)+
            ggtitle("Cell cycle")+
            theme(axis.text.x = element_text(face="bold",  size=14),
                  axis.text.y = element_text(face="bold",size=14),
                  axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
                  axis.title = element_text(size = rel(1.5),face="bold"),
                  plot.title = element_text(color="black", size=20, 
                                            face="bold.italic",hjust = 0.1,vjust =-10), # margin = margin(t = 0.5, b = -7),
                  legend.title = element_text(size=12, color = "black", face="bold"),
                  legend.text = element_text(colour="black", size=12,face="bold"),
                  legend.position = c(0.1, 0.18),
                  aspect.ratio=1) + #square plot
            theme(axis.line.x = element_line(colour = "black", size = 0.8),
                  axis.line.y = element_line(colour = "black", size = 0.8))
          
          ## DucT2_TOP2ACente
          plot_cells(cds_sub_DucT2_TOP2ACenter, color_cells_by= "cell_cycle", 
                     label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 2)  + 
            scale_color_manual(values = colors_cc)+
            ggtitle("Cell cycle")+
            theme(axis.text.x = element_text(face="bold",  size=14),
                  axis.text.y = element_text(face="bold",size=14),
                  axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
                  axis.title = element_text(size = rel(1.5),face="bold"),
                  plot.title = element_text(color="black", size=20, 
                                            face="bold.italic",hjust = 0.1,vjust =-10), # margin = margin(t = 0.5, b = -7),
                  legend.title = element_text(size=12, color = "black", face="bold"),
                  legend.text = element_text(colour="black", size=12,face="bold"),
                  legend.position = c(0.1, 0.18),
                  aspect.ratio=1) + #square plot
            theme(axis.line.x = element_line(colour = "black", size = 0.8),
                  axis.line.y = element_line(colour = "black", size = 0.8))
          
          ##!!! candidates14
          Main = candidates14
          plot_cells(cds_sub_AcinaDucT, genes = Main, show_trajectory_graph = F,label_cell_groups = F
                     ,cell_size =0.5)+
            scale_colour_gradient2(low = "black", mid = "#3528c7", high = "#ff85e9", 
                                   guide = "colourbar",midpoint = 0, labs(fill ="Exp")) +
            #      ggtitle(Main)+
            
            theme(axis.text.x = element_text(face="bold",  size=10),
                  axis.text.y = element_text(face="bold",size=10),
                  axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
                  axis.title = element_text(size = rel(1.5),face="bold"),
                  plot.title = element_text(color="black", size=20, 
                                            face="bold.italic",hjust = 0.1,vjust =-10), # margin = margin(t = 0.5, b = -7),
                  #     plot.background = element_rect(fill = 'chartreuse'),
                  legend.title = element_text(size=12, color = "black", face="bold"),
                  legend.text = element_text(colour="black", size=12,face="bold"),
                  legend.background = element_rect(fill = alpha("white", 0.5)),
                  legend.direction = "horizontal",
                  legend.position = c(0.75, 0.08),
                  #     plot.text = element_text(size = 20),
                  aspect.ratio=1) + #square plot
            theme(axis.line.x = element_line(colour = "black", size = 0.6),
                  axis.line.y = element_line(colour = "black", size = 0.6))
          
          
          ##!!!!!! Sub (Patients)
          Main = "TOP2A"
          Sub_Name = "Patient1"
          Sub = cds_Patient1
          #######   Cell cycle  #######           
          plot_cells(Sub, color_cells_by= "cell_cycle", 
                     label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1)  + 
            scale_color_manual(values = colors_cc)+
            ggtitle(paste0("Cell cycle","(",Sub_Name,")"))+
            theme(axis.text.x = element_text(face="bold",  size=14),
                  axis.text.y = element_text(face="bold",size=14),
                  axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
                  axis.title = element_text(size = rel(1.5),face="bold"),
                  plot.title = element_text(color="black", size=17, 
                                            face="bold.italic",hjust = 0.1,vjust =-7), # margin = margin(t = 0.5, b = -7),
                  legend.title = element_text(size=12, color = "black", face="bold"),
                  legend.text = element_text(colour="black", size=12,face="bold"),
                  legend.position = c(0.1, 0.18),
                  aspect.ratio=1) + #square plot
            theme(axis.line.x = element_line(colour = "black", size = 0.8),
                  axis.line.y = element_line(colour = "black", size = 0.8))
          
          #######   Recluster  #######           
          plot_cells(cds_sub_AcinaDucT_NewK_ReCluster, color_cells_by= "ReCluster", label_cell_groups=T,group_label_size = 4.5, 
                     show_trajectory_graph = FALSE,cell_size = 1)  + 
#            scale_color_manual(values = colors_cc)+
#            ggtitle("cluster")+
            theme(axis.text.x = element_text(face="bold",  size=14),
                  axis.text.y = element_text(face="bold",size=14),
                  axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
                  axis.title = element_text(size = rel(1.5),face="bold"),
                  plot.title = element_text(color="black", size=20, 
                                            face="bold.italic",hjust = 0.1,vjust =-10), # margin = margin(t = 0.5, b = -7),
                  legend.title = element_text(size=12, color = "black", face="bold"),
                  legend.text = element_text(colour="black", size=12,face="bold"),
#                  legend.position = c(0.1, 0.18),
                  aspect.ratio=1) + #square plot
            theme(axis.line.x = element_line(colour = "black", size = 0.8),
                  axis.line.y = element_line(colour = "black", size = 0.8))
         
          #######   Recluster V2  #######           
          plot_cells(cds_sub_AcinaDucT_NewK_ReCluster, color_cells_by= "ReCluster", label_cell_groups=F,group_label_size = 4.5, 
                     show_trajectory_graph = FALSE,cell_size = 1)  + 
            #            scale_color_manual(values = colors_cc)+
            #            ggtitle("cluster")+
            theme(axis.text.x = element_text(face="bold",  size=14),
                  axis.text.y = element_text(face="bold",size=14),
                  axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
                  axis.title = element_text(size = rel(1.5),face="bold"),
                  plot.title = element_text(color="black", size=20, 
                                            face="bold.italic",hjust = 0.1,vjust =-10), # margin = margin(t = 0.5, b = -7),
                  legend.title = element_text(size=12, color = "black", face="bold"),
                  legend.text = element_text(colour="black", size=9),
                  #                  legend.position = c(0.1, 0.18),
                  aspect.ratio=1) + #square plot
            theme(axis.line.x = element_line(colour = "black", size = 0.8),
                  axis.line.y = element_line(colour = "black", size = 0.8))
          
          
          #######   cluster  #######           
          plot_cells(cds_sub_AcinaDucT2, color_cells_by= "cluster", label_cell_groups=FALSE, 
                     show_trajectory_graph = T,cell_size = 1)  + 
            #            scale_color_manual(values = colors_cc)+
            #            ggtitle("cluster")+
            theme(axis.text.x = element_text(face="bold",  size=14),
                  axis.text.y = element_text(face="bold",size=14),
                  axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
                  axis.title = element_text(size = rel(1.5),face="bold"),
                  plot.title = element_text(color="black", size=20, 
                                            face="bold.italic",hjust = 0.1,vjust =-10), # margin = margin(t = 0.5, b = -7),
                  legend.title = element_text(size=12, color = "black", face="bold"),
                  legend.text = element_text(colour="black", size=12,face="bold"),
                  #                  legend.position = c(0.1, 0.18),
                  aspect.ratio=1) + #square plot
            theme(axis.line.x = element_line(colour = "black", size = 0.8),
                  axis.line.y = element_line(colour = "black", size = 0.8))

          #######   pseudotime  #######           
          plot_cells(cds_sub_AcinaDucT2, color_cells_by= "pseudotime", label_cell_groups=FALSE, 
                     show_trajectory_graph = T,cell_size = 1)  + 
            #            scale_color_manual(values = colors_cc)+
            #            ggtitle("pseudotime")+
            theme(axis.text.x = element_text(face="bold",  size=14),
                  axis.text.y = element_text(face="bold",size=14),
                  axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
                  axis.title = element_text(size = rel(1.5),face="bold"),
                  plot.title = element_text(color="black", size=20, 
                                            face="bold.italic",hjust = 0.1,vjust =-10), # margin = margin(t = 0.5, b = -7),
                  legend.title = element_text(size=12, color = "black", face="bold"),
                  legend.text = element_text(colour="black", size=12,face="bold"),
                  #                  legend.position = c(0.1, 0.18),
                  aspect.ratio=1) + #square plot
            theme(axis.line.x = element_line(colour = "black", size = 0.8),
                  axis.line.y = element_line(colour = "black", size = 0.8))
          
          
          
          #######   Marker  ####### 
          ### PDAC 
          plot_cells(cds_sub_AcinaDucT, color_cells_by= Marker_PDAC_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1.2) +
            scale_colour_gradient2(low = "#440075", mid = "#ffd261", high = "#4aff8c", 
                                   guide = "colourbar",midpoint = 0.2, labs(fill =  "Exp"))+
            ggtitle(Marker_PDAC_Name)+
            theme(axis.text.x = element_text(face="bold",  size=14),
                  axis.text.y = element_text(face="bold",size=14),
                  axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
                  axis.title = element_text(size = rel(1.5),face="bold"),
                  plot.title = element_text(color="black", size=20, 
                                            face="bold.italic",hjust = 0.1,vjust =-10), # margin = margin(t = 0.5, b = -7),
                  legend.title = element_text(size=12, color = "black", face="bold"),
                  legend.text = element_text(colour="black", size=12,face="bold"),
                  legend.position = c(0.1, 0.18),
                  aspect.ratio=1) + #square plot
           theme(axis.line.x = element_line(colour = "black", size = 0.8),
                  axis.line.y = element_line(colour = "black", size = 0.8))
          
          ### EMT 
          plot_cells(cds_sub_AcinaDucT, color_cells_by= Marker_EMT_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1.2) +
            scale_colour_gradient2(low = "#440075", mid = "#ffd261", high = "#4aff8c", 
                                   guide = "colourbar",midpoint = 0.2, labs(fill =  "Exp"))+
            ggtitle(Marker_EMT_Name)+
            theme(axis.text.x = element_text(face="bold",  size=14),
                  axis.text.y = element_text(face="bold",size=14),
                  axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
                  axis.title = element_text(size = rel(1.5),face="bold"),
                  plot.title = element_text(color="black", size=20, 
                                            face="bold.italic",hjust = 0.1,vjust =-10), # margin = margin(t = 0.5, b = -7),
                  legend.title = element_text(size=12, color = "black", face="bold"),
                  legend.text = element_text(colour="black", size=12,face="bold"),
                  legend.position = c(0.1, 0.18),
                  aspect.ratio=1) + #square plot
            theme(axis.line.x = element_line(colour = "black", size = 0.8),
                  axis.line.y = element_line(colour = "black", size = 0.8))
          
         
          ### ChroSt 
          plot_cells(cds_sub_AcinaDucT, color_cells_by= Marker_ChroSt_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1.2) +
            scale_colour_gradient2(low = "#440075", mid = "#ffd261", high = "#4aff8c", 
                                   guide = "colourbar",midpoint = 0.15, labs(fill =  "Exp"))+
            ggtitle(Marker_ChroSt_Name)+
            theme(axis.text.x = element_text(face="bold",  size=14),
                  axis.text.y = element_text(face="bold",size=14),
                  axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
                  axis.title = element_text(size = rel(1.5),face="bold"),
                  plot.title = element_text(color="black", size=20, 
                                            face="bold.italic",hjust = 0.1,vjust =-10), # margin = margin(t = 0.5, b = -7),
                  legend.title = element_text(size=12, color = "black", face="bold"),
                  legend.text = element_text(colour="black", size=12,face="bold"),
                  legend.position = c(0.1, 0.18),
                  aspect.ratio=1) + #square plot
            theme(axis.line.x = element_line(colour = "black", size = 0.8),
                  axis.line.y = element_line(colour = "black", size = 0.8))
          
          ### Migration 
          plot_cells(cds_sub_AcinaDucT, color_cells_by= Marker_Mig_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1.2) +
            scale_colour_gradient2(low = "#440075", mid = "#ffd261", high = "#4aff8c", 
                                   guide = "colourbar",midpoint = 0.15, labs(fill =  "Exp"))+
            ggtitle(Marker_Mig_Name)+
            theme(axis.text.x = element_text(face="bold",  size=14),
                  axis.text.y = element_text(face="bold",size=14),
                  axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
                  axis.title = element_text(size = rel(1.5),face="bold"),
                  plot.title = element_text(color="black", size=20, 
                                            face="bold.italic",hjust = 0.1,vjust =-10), # margin = margin(t = 0.5, b = -7),
                  legend.title = element_text(size=12, color = "black", face="bold"),
                  legend.text = element_text(colour="black", size=12,face="bold"),
                  legend.position = c(0.1, 0.18),
                  aspect.ratio=1) + #square plot
            theme(axis.line.x = element_line(colour = "black", size = 0.8),
                  axis.line.y = element_line(colour = "black", size = 0.8))
          
          ### Metastasis 
          plot_cells(cds_sub_AcinaDucT, color_cells_by= Marker_Meta_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1.2) +
            scale_colour_gradient2(low = "#440075", mid = "#ffd261", high = "#4aff8c", 
                                   guide = "colourbar",midpoint = 0.15, labs(fill =  "Exp"))+
            ggtitle(Marker_Meta_Name)+
            theme(axis.text.x = element_text(face="bold",  size=14),
                  axis.text.y = element_text(face="bold",size=14),
                  axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
                  axis.title = element_text(size = rel(1.5),face="bold"),
                  plot.title = element_text(color="black", size=20, 
                                            face="bold.italic",hjust = 0.1,vjust =-10), # margin = margin(t = 0.5, b = -7),
                  legend.title = element_text(size=12, color = "black", face="bold"),
                  legend.text = element_text(colour="black", size=12,face="bold"),
                  legend.position = c(0.1, 0.18),
                  aspect.ratio=1) + #square plot
            theme(axis.line.x = element_line(colour = "black", size = 0.8),
                  axis.line.y = element_line(colour = "black", size = 0.8))

          
          
          ### ATR 
          plot_cells(cds_sub_AcinaDucT_NewK_ReCluster, color_cells_by= Marker_ATR_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1.2) +
            scale_colour_gradient2(low = "#440075", mid = "#ffd261", high = "#4aff8c", 
                                   guide = "colourbar",midpoint = 0.12, labs(fill =  "Exp"))+
            ggtitle(Marker_ATR_Name)+
            theme(axis.text.x = element_text(face="bold",  size=14),
                  axis.text.y = element_text(face="bold",size=14),
                  axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
                  axis.title = element_text(size = rel(1.5),face="bold"),
                  plot.title = element_text(color="black", size=20, 
                                            face="bold.italic",hjust = 0.1,vjust =-10), # margin = margin(t = 0.5, b = -7),
                  legend.title = element_text(size=12, color = "black", face="bold"),
                  legend.text = element_text(colour="black", size=12,face="bold"),
                  legend.position = c(0.1, 0.18),
                  aspect.ratio=1) + #square plot
            theme(axis.line.x = element_line(colour = "black", size = 0.8),
                  axis.line.y = element_line(colour = "black", size = 0.8))
          
          ### ACST 
          plot_cells(cds_sub_AcinaDucT_NewK_ReCluster, color_cells_by= Marker_ACST_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1.2) +
            scale_colour_gradient2(low = "#440075", mid = "#ffd261", high = "#4aff8c", 
                                   guide = "colourbar",midpoint = 0.12, labs(fill =  "Exp"))+
            ggtitle(Marker_ACST_Name)+
            theme(axis.text.x = element_text(face="bold",  size=14),
                  axis.text.y = element_text(face="bold",size=14),
                  axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
                  axis.title = element_text(size = rel(1.5),face="bold"),
                  plot.title = element_text(color="black", size=20, 
                                            face="bold.italic",hjust = 0.1,vjust =-10), # margin = margin(t = 0.5, b = -7),
                  legend.title = element_text(size=12, color = "black", face="bold"),
                  legend.text = element_text(colour="black", size=12,face="bold"),
                  legend.position = c(0.1, 0.18),
                  aspect.ratio=1) + #square plot
            theme(axis.line.x = element_line(colour = "black", size = 0.8),
                  axis.line.y = element_line(colour = "black", size = 0.8))
          
         ### Marker          
          plot_cells(cds_sub_AcinaDucT, color_cells_by= Marker_Meta_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1.2) +
            scale_colour_gradient2(low = "darkblue", mid = "#f7c211", high = "green", 
                                   guide = "colourbar",midpoint = 0.15, labs(fill = Marker_Meta_Name))+
            theme(axis.text.x = element_text(face="bold",  size=14),
                  axis.text.y = element_text(face="bold",size=14),
                  axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
                  axis.title = element_text(size = rel(1.5),face="bold"),
                  plot.title = element_text(color="black", size=20, face="bold.italic",hjust = 0.5,vjust =-1),
                  legend.title = element_text(size=12, color = "black", face="bold"),
                  legend.text = element_text(colour="black", size=12,face="bold"),
                  legend.background = element_rect(fill = alpha("white", 0.5)),
                  legend.position = c(0.12, 0.18),
                  aspect.ratio=1) + #square plot
            theme(axis.line.x = element_line(colour = "black", size = 1),
                  axis.line.y = element_line(colour = "black", size = 1))
    
    ##########################################
          ##!!!!!! 
          Main = "CDT1"
          Main = "CDC6"
          Main = "RAD17"
          Main = "GMNN"
          
          
          plot_cells(cds_sub_AcinaDucT_NewK_ReCluster, genes = Main, show_trajectory_graph = F,label_cell_groups = F
                     ,cell_size =1 )+
            scale_colour_gradient2(low = "black", mid = "#3528c7", high = "#ff85e9", 
                                   guide = "colourbar",midpoint = 0, labs(fill ="Exp")) +
            ggtitle(Main)+
            theme(axis.text.x = element_text(face="bold",  size=14),
                  axis.text.y = element_text(face="bold",size=14),
                  axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
                  axis.title = element_text(size = rel(1.5),face="bold"),
                  plot.title = element_text(color="black", size=17, 
                                            face="bold.italic",hjust = 0.1,vjust =-10), # margin = margin(t = 0.5, b = -7),
                  #     plot.background = element_rect(fill = 'chartreuse'),
                  legend.title = element_text(size=12, color = "black", face="bold"),
                  legend.text = element_text(colour="black", size=12,face="bold"),
                  legend.background = element_rect(fill = alpha("white", 0.5)),
                  legend.position = c(0.1, 0.18),
                  #     plot.text = element_text(size = 20),
                  aspect.ratio=1) + #square plot
            theme(axis.line.x = element_line(colour = "black", size = 0.8),
                  axis.line.y = element_line(colour = "black", size = 0.8))
          

