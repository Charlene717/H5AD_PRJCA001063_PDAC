##### Current path and new folder setting* #####
ProjectName = "TOP2A"
Sampletype = "PDAC"
#ProjSamp.Path = paste0(Sampletype,"_",ProjectName)

Version = paste0(Sys.Date(),"_",ProjectName,"_",Sampletype)
Save.Path = paste0(getwd(),"/",Version)
## Create new folder
if (!dir.exists(Save.Path)){
  dir.create(Save.Path)
}

IGene = "TOP2A"

#####
cds_subset <- cds[,colData(cds)$Cell_type %in% c("Acinar cell","Ductal cell type 1","Ductal cell type 2")]
cds_subset <- choose_cells(cds)
plot_cells(cds_subset, genes=c("TOP2A"),show_trajectory_graph = FALSE)
plot_cells(cds_subset, genes=c("TOP2A"),show_trajectory_graph = FALSE,norm_method = c("log"))

library(monocle3)
library(ggplot2)
library(tidyverse)
Main = c("TOP2A","H2AX","PTK2","NSUN2","TP53","MYC")
pdf(
  file = paste0(getwd(), "/",Version,"/", Sys.Date(), "_GeneUMAP_log_STR.pdf"),
  width = 7, height = 7
)
for (i in 1:length(Main)) {
  plot_cells(cds_subset, genes = Main[i], show_trajectory_graph = F,label_cell_groups = F
             ,cell_size =1,norm_method = c("log"),scale_to_range=T)+ # norm_method = c("log", "size_only"),
    scale_colour_gradient2(low = "#0c3999", mid = "#7dbd91", high = "#edde87", 
                           guide = "colourbar",midpoint = 0, labs(fill ="Exp")) +
    ggtitle(Main[i])+
    #ggtitle(paste0(Main,"(",Sub_Name,")"))+            
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
          axis.line.y = element_line(colour = "black", size = 0.8)) -> P
  print(P)
  
}
graphics.off()

##### Pesudotime test #####
cds <- learn_graph(cds)
plot_cells(cds_subset,
           color_cells_by = "Cell_type",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

#cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
cds <- order_cells(cds)
cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "Cell_type",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)


cds_subTTT <- choose_graph_segments(cds)
cds_subTTT <- cds[,colData(cds)$ReCluster %in% c("AD","AC","aAtD")]
cds_subTTT <- order_cells(cds_subTTT)
#cds_subTTT <- learn_graph(cds_subTTT)



AFD_genes <- c("TOP2A", "PTK2", "NSUN2")
AFD_lineage_cds <- cds_subTTT[rowData(cds_subTTT)$gene_short_name %in% AFD_genes,]
plot_genes_in_pseudotime(AFD_lineage_cds,
                         color_cells_by="ReCluster",
                         min_expr=0.5)
plot_genes_in_pseudotime(AFD_lineage_cds,
                         color_cells_by="EMT",
                         min_expr=0.5)
plot_cells(cds_subTTT,
           color_cells_by = "Cell_type",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

##### Correlation #####

## Extract dataframe
Anno.df <- data.frame(PDAC = cds@colData@listData[["PDAC"]],
                      EMT = cds@colData@listData[["EMT"]],
                      CST =  cds@colData@listData[["CST"]],
                      Migration = cds@colData@listData[["Migration"]],
                      Metastasis = cds@colData@listData[["Metastasis"]],
                      NE = cds@colData@listData[["NE"]],
                      NP = cds@colData@listData[["NP"]],
                      Cell_type = cds@colData@listData[["Cell_type"]],
                      ReCluster = cds@colData@listData[["ReCluster"]],
                      CELL = cds@colData@listData[["CELL"]])

matrix.df <- cds@assays@data@listData[["counts"]] %>% as.data.frame
matrix.df  <- data.frame(gene = row.names(matrix.df),matrix.df)
matrixlog.df <- cds@assays@data@listData[["logcounts"]] %>% as.data.frame()
matrixlog.df  <- data.frame(gene = row.names(matrixlog.df),matrixlog.df)

## Export tsv files 
write.table(Anno.df,"Anno.tsv",row.names = F,sep = "\t")
write.table(matrix.df ,"matrix.tsv",row.names = F,sep = "\t")
write.table(matrixlog.df ,"matrixlog.tsv",row.names = F,sep = "\t")

##### GSEA test #####
