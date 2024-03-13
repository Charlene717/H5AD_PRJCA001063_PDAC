## Ref: https://github.com/sqjin/CellChat
## Ref: https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html

##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

##### Parameter setting* #####
  SignalingType = "All" # Secreted Signaling, ECM-Receptor, Cell-Cell Contact, All
  Species = "Human" # Human, Mouse
  nPatternsOut = 4 # Set patterns number for identify global communication in outgoing signaling 
  nPatternsIn = 4
  
##### Current path and new folder setting*  #####
  ProjectName = "All" # Secret, ECM, CC, All
  Version = paste0(Sys.Date(),"_",ProjectName,"_PADC")
  Save.Path = paste0(getwd(),"/",Version)
  ## Create new folder
  if (!dir.exists(Save.Path)){
    dir.create(Save.Path)
  }
  
#### Load the required libraries ####
  #### Basic installation ####
  ## Package.set
  Package.set <- c("tidyverse","patchwork","NMF","ggalluvial")
  ## Check whether the installation of those packages is required
  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      install.packages(Package.set[i])
    }
  }
  ## Load Packages
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)
  
  #### BiocManager installation ####
  ## Package.set
  Package.set <- c("ComplexHeatmap")
  ## Check whether the installation of those packages is required from BiocManager
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      BiocManager::install(Package.set[i])
    }
  }
  ## Load Packages
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)

  #### GitHub installation ####
  if (!require("remotes", quietly = TRUE))
    install.packages("remotes")
  devtools::install_github("sqjin/CellChat")
  remotes::install_github("mojaveazure/seurat-disk")
  library(SeuratDisk)
  
  # if (!require("devtools", quietly = TRUE))
  #   install.packages("devtools")
  # devtools::install_github("satijalab/seurat-data")
  # library(SeuratData)

  
  
##### Part I: Data input & processing and initialization of CellChat object #####  
  ####------ Load data ------####
  #### Converse h5ad to Seurat ####
  library(SeuratDisk)
  
  # This creates a copy of this .h5ad object reformatted into .h5seurat inside the example_dir directory
  Convert("StdWf1_PRJCA001063_CRC_besca2.annotated.h5ad", "PRJCA001063.h5seurat")
  
  # This .d5seurat object can then be read in manually
  seuratObject <- LoadH5Seurat("PRJCA001063.h5seurat")
  
  #### Extract the CellChat input files from a Seurat V3 object ####
    # Ref: https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/Interface_with_other_single-cell_analysis_toolkits.html
    library(Seurat)
    data.input <- GetAssayData(seuratObject, assay = "RNA", slot = "data") # normalized data matrix
    labels <- Idents(seuratObject)
    # meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
    meta <- seuratObject@meta.data # create a dataframe of the cell labels

    # ## Prepare input data for CelChat analysis 
    # cell.use = rownames(meta)[meta$condition == "LS"] # extract the cell names from disease data
    # data.input = data.input[, cell.use]
    # meta = meta[cell.use, ]
    # meta = data.frame(labels = meta$labels[cell.use], row.names = colnames(data.input)) # manually create a dataframe consisting of the cell labels
    unique(meta$Cell_type) # check the cell labels
  
  #### Create a CellChat object ####
    cellchat <- createCellChat(object = data.input, meta = meta, group.by = "Cell_type")
    
  #### Set the ligand-receptor interaction database ####
    if(Species == "Human"){
      CellChatDB <- CellChatDB.human 
    }else if(Species == "Mouse"){
      CellChatDB <- CellChatDB.mouse 
    }else{
      print("Error in Species setting: Please set the Species as Human or Mouse.")
    }
    
    showDatabaseCategory(CellChatDB)
    
    # Show the structure of the database
    dplyr::glimpse(CellChatDB$interaction)
  
    # use a subset of CellChatDB for cell-cell communication analysis
    if(SignalingType == "All"){
      # use all CellChatDB for cell-cell communication analysis
      CellChatDB.use <- CellChatDB # simply use the default CellChatDB
    }else{
      CellChatDB.use <- subsetDB(CellChatDB, search = SignalingType) 
      # CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
    }
    
    # set the used database in the object
    cellchat@DB <- CellChatDB.use
  
  #### Preprocessing the expression data for cell-cell communication analysis ####
    # subset the expression data of signaling genes for saving computation cost
    cellchat <- subsetData(cellchat) # Subset the expression data of signaling genes for saving computation cost. This step is necessary even if using the whole database
    future::plan("multiprocess", workers = 4) # do parallel
  
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
  
    # project gene expression data onto PPI network (optional)
    cellchat <- projectData(cellchat, PPI.human)

##### Part II: Inference of cell-cell communication network #####    
  #### Compute the communication probability and infer cellular communication network ####
    cellchat <- computeCommunProb(cellchat)
    # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    
  #### Infer the cell-cell communication at a signaling pathway level ####
    cellchat <- computeCommunProbPathway(cellchat)

  #### Calculate the aggregated cell-cell communication network ####  
    cellchat <- aggregateNet(cellchat)  
    
    groupSize <- as.numeric(table(cellchat@idents))
    
    ## Circle plot
    par(mfrow = c(1,2), xpd=TRUE)
    netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
    netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
    

##### Part III: Visualization of cell-cell communication network #####
  ##### Summary #####  
    ## Create new folder
    PathSum <- paste0(Save.Path,"/Summary")
    if (!dir.exists(PathSum)){
      dir.create(PathSum)
    }
    
    #### Visualize summarize signaling pathway using Hierarchy plot, Circle plot or Chord diagram ####
    groupSize <- as.numeric(table(cellchat@idents))
    pdf(file = paste0(PathSum,"/",ProjectName,"_Sum_Communication_Network_Main.pdf"),
        width = 7,  height = 7)
      
      ## Circle plot
      netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
      netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
      
      ## Heatmap
      par(mfrow=c(1,1))
      netVisual_heatmap(cellchat, color.heatmap = "Reds")
      
      ## Chord diagram  
      par(mfrow=c(1,1))
      netVisual_aggregate(cellchat, signaling = cellchat@netP[["pathways"]], layout = "chord")
      
      ## Barplot: Contribution of each ligand-receptor
      netAnalysis_contribution(cellchat, signaling = cellchat@netP[["pathways"]])
      
    dev.off()
    
    ## CirclePlot Sup
    pdf(file = paste0(PathSum,"/",ProjectName,"_Sum_Communication_Network_CirclePlot_Sup.pdf"),
        width = 15,  height = 12)
    par(mfrow = c(1,2), xpd=TRUE)
    gg1 <-netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
    gg2 <-netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
    
    mat <- cellchat@net$weight
    par(mfrow = c(3,4), xpd=TRUE)
    for (i in 1:nrow(mat)) {
      mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
      mat2[i, ] <- mat[i, ]
      netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
    }
    dev.off()
    rm(i)
    
    #### Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways ####  
    ## Summary Bubble
    pdf(file = paste0(PathSum,"/",ProjectName,"_Sum_Communication_Network_Bubble.pdf"),
        width = 15,  height = 20)
      ## Bubble plot
      # show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
      netVisual_bubble(cellchat,remove.isolate = FALSE)
    dev.off()
    
    ## Summary Violin
    pdf(file = paste0(PathSum,"/",ProjectName,"_Sum_Communication_Network_Violin.pdf"),
        width = 10,  height = 20)
    ## Violin Plot
    plotGeneExpression(cellchat, signaling = cellchat@netP[["pathways"]])
    dev.off()
    
    

    
  #### All pathways #### 
    ## Create new folder
    PathDetail <- paste0(Save.Path,"/Detail")
    if (!dir.exists(PathDetail)){
      dir.create(PathDetail)
    }
    
    pathway.set <- cellchat@netP[["pathways"]]
    
    #### Visualize each signaling pathway using Hierarchy plot, Circle plot or Chord diagram ####
      #### Main: Plot all pathway ####
      pdf(file = paste0(PathDetail,"/",ProjectName,"_AllPT_Communication_Network_Main_01_AllPT.pdf"),
          width = 7,  height = 7
      )
      
        for (i in 1:length(pathway.set)) {
        pathways.show <- pathway.set[i] # pathways.show <- c("CXCL") 
        
        # # Hierarchy plot
        # # Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
        # vertex.receiver = seq(1,4) # a numeric vector. 
        # netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
        
        # Circle plot
        par(mfrow=c(1,1))
        netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
        
        # Heatmap
        par(mfrow=c(1,1))
        Heatmap <- netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
        print(Heatmap)
        #> Do heatmap based on a single object
        
        
        # Chord diagram
        par(mfrow=c(1,1))
        netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
        
        ### Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
        # Barchart
        p <- netAnalysis_contribution(cellchat, signaling = pathways.show)
        print(p)
        
      }
    dev.off()
    rm(i,p,Heatmap)
      
      #### Main: Plot all pathway and LR ####
        for (i in 1:length(pathway.set)) {
          pathways.show <- pathway.set[i] # pathways.show <- c("CXCL") 
          
          pdf(file = paste0(PathDetail,"/",ProjectName,"_AllPT_Communication_Network_Main_",pathways.show,"_LR.pdf"),
              width = 7,  height = 7
          )
  
            # # Hierarchy plot
            # # Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
            # vertex.receiver = seq(1,4) # a numeric vector. 
            # netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
            
            # Circle plot
            par(mfrow=c(1,1))
            netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
  
            # Heatmap
            par(mfrow=c(1,1))
            Heatmap <- netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
            print(Heatmap)
            #> Do heatmap based on a single object
            
            
            # Chord diagram
            par(mfrow=c(1,1))
            netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
            
            # # Chord diagram
            # group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
            # names(group.cellType) <- levels(cellchat@idents)
            # netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))
            # 
            # #> Plot the aggregated cell-cell communication network at the signaling pathway level
            # #> Note: The first link end is drawn out of sector 'Inflam. FIB'.
            
            
            ### Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
            # Barchart
            p <- netAnalysis_contribution(cellchat, signaling = pathways.show)
            print(p)
            
            pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
            for (j in 1:nrow(pairLR.CXCL)) {
               LR.show <- pairLR.CXCL[j,] # show one ligand-receptor pair
               # # Hierarchy plot
               # vertex.receiver = seq(1,4) # a numeric vector
               # netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
               
               # Circle plot
               netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
               # Chord diagram
               netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
               
            }
            rm(j)
          
          dev.off()
        }
      rm(i,p,Heatmap)
        
    #### Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways ####  
      #### Bubble plot ####
      ## Plot the signaling gene expression distribution using dot plot
      # show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
      netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)
  
      ## Sum
      pdf(file = paste0(PathDetail,"/",ProjectName,"_AllLRPair_Bubble_Sum.pdf"),
          width = 15,  height = 20
      )
        netVisual_bubble(cellchat,  remove.isolate = FALSE)
      dev.off()
      
  
      ## All
      pdf(file = paste0(PathDetail,"/",ProjectName,"_AllLRPair_Bubble_All.pdf"),
          width = 5,  height = 10
      )
        for (i in 1:ncol(mat)) {
        try({
    
        # show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
        P <- netVisual_bubble(cellchat, sources.use = i, remove.isolate = FALSE)
        print(P)
  
      })
      }
      dev.off()
      rm(i,p)
      
      #### Chord diagram ####
      # show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
      pdf(file = paste0(PathDetail,"/",ProjectName,"_AllLRPair_ChordDiagram.pdf"),
          width = 10,  height = 10
      )
      for (i in 1:ncol(mat)) {
      try({
        
        # show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
        P1 <- netVisual_chord_gene(cellchat, sources.use = i,  lab.cex = 0.5,legend.pos.y = 30, title.name = paste0("Signaling from ",levels(cellchat@idents)[i]))
        print(P1)
        P2 <- netVisual_chord_gene(cellchat, targets.use =i , lab.cex = 0.5,legend.pos.y = 30, title.name = paste0("Signaling received by ",levels(cellchat@idents)[i]))
        print(P2)
      })
      }
      dev.off()
      rm(i,P1,P2)
      
      #### Violin ####
      ## Plot the signaling gene expression distribution using violin plot
      pdf(file = paste0(PathDetail,"/",ProjectName,"_AllLRPair_Violin.pdf"),
          width = 10,  height = 10
      )
        for (i in 1:length(pathway.set)) {
        pathways.show <- pathway.set[i] # pathways.show <- c("CXCL") 
          
        P <- plotGeneExpression(cellchat, signaling = pathways.show)
        print(P)
        }
      dev.off()
      rm(i,P)
  
##### Part IV: Systems analysis of cell-cell communication network #####
  ## Create new folder
  PathSys <- paste0(Save.Path,"/SysAna")
  if (!dir.exists(PathSys)){
    dir.create(PathSys)
  }
  ##### Identify signaling roles (e.g., dominant senders, receivers) of cell groups as well as the major contributing signaling ##### 
    #### Compute and visualize the network centrality scores ####
      # Compute the network centrality scores
      cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
      # Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
      netAnalysis_signalingRole_network(cellchat, signaling = pathway.set[1], width = 8, height = 2.5, font.size = 10)
      
      pdf(file = paste0(PathSys,"/",ProjectName,"_SystemsAnalysis_Heatmap_NWCentralityScores.pdf"),
          width = 7,  height = 7
      )
        for (i in 1:length(pathway.set)) {
          pathways.show <- pathway.set[i] # pathways.show <- c("CXCL") 
          
          P <-  netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
          print(P)
        }
      dev.off()
      rm(i,P)

    #### Visualize the dominant senders (sources) and receivers (targets) in a 2D space #### 
      # Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
      gg1 <- netAnalysis_signalingRole_scatter(cellchat) + ggtitle("All signaling pathways")+ 
        theme(plot.title = element_text(color="black", size=14, face="bold"))
      # Signaling role analysis on the cell-cell communication networks of interest
      gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = pathway.set[1]) + 
        ggtitle(paste0(pathway.set[1]," signaling pathway network"))+ 
        theme(plot.title = element_text(color="black", size=14, face="bold"))
      gg1 + gg2
      
      
      pdf(file = paste0(PathSys,"/",ProjectName,"_SystemsAnalysis_SourcesTargets_2Dspace.pdf"),
          width = 7,  height = 7
      )
        # Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
        gg1 <- netAnalysis_signalingRole_scatter(cellchat) + ggtitle("All signaling pathways")+ 
          theme(plot.title = element_text(color="black", size=14, face="bold"))
        gg1
        
        for (i in 1:length(pathway.set)) {
          pathways.show <- pathway.set[i] # pathways.show <- c("CXCL") 
          #> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
          # Signaling role analysis on the cell-cell communication networks of interest
          gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = pathways.show) + 
            ggtitle(paste0(pathways.show," signaling pathway network"))+ 
            theme(plot.title = element_text(color="black", size=14, face="bold"))
          print(gg2)
          rm(gg2)
      }
      dev.off()
      rm(i,gg1,gg2)
      
      
    #### Identify signals contributing most to outgoing or incoming signaling of certain cell groups ####
      # Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
      ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
      ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
      ht1 + ht2
      # Signaling role analysis on the cell-cell communication networks of interest
      ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = pathways.show)
      ht   
      
      pdf(file = paste0(PathSys,"/",ProjectName,"_SystemsAnalysis_Heatmap_mostOutIn.pdf"),
          width = 12,  height = 8
      )
        # Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
        ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
        ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
        ht1 + ht2
        
        for (i in 1:length(pathway.set)) {
    
          pathways.show <- pathway.set[i] # pathways.show <- c("CXCL") 
     
          # Signaling role analysis on the cell-cell communication networks of interest
          try({
          ht3 <- netAnalysis_signalingRole_heatmap(cellchat, signaling = pathways.show, pattern = "outgoing")
          print(ht3)
          })
          try({
          ht4 <- netAnalysis_signalingRole_heatmap(cellchat, signaling = pathways.show, pattern = "incoming")
          print(ht4)
          })
          rm(ht3,ht4)
        }
      dev.off()
      rm(i,ht1,ht2)
  
  ##### Identify global communication patterns to explore how multiple cell types and signaling pathways coordinate together #####
    #### Identify and visualize outgoing communication pattern of secreting cells ####    
      library(NMF)
      library(ggalluvial)
      P.outgoing <- selectK(cellchat, pattern = "outgoing")
      P.incoming <- selectK(cellchat, pattern = "incoming")
      
      pdf(file = paste0(PathSys,"/",ProjectName,"_SystemsAnalysis_GlobalPatterns_outgoing.pdf"),
          width = 12,  height = 8
      )
  
        cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatternsOut)
        # P.nHeatmap.outgoing <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
        # P.nHeatmap.outgoing
        # river plot
        netAnalysis_river(cellchat, pattern = "outgoing")
        #> Please make sure you have load `library(ggalluvial)` when running this function 
           
        # dot plot
        netAnalysis_dot(cellchat, pattern = "outgoing")
        
        P.outgoing
        #graphics.off()
      dev.off()
      
    #### Identify and visualize incoming communication pattern of target cells ####
      
      pdf(file = paste0(PathSys,"/",ProjectName,"_SystemsAnalysis_GlobalPatterns_incoming.pdf"),
          width = 12,  height = 8
      )
        cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatternsIn)
        
        # river plot
        netAnalysis_river(cellchat, pattern = "incoming")
        #> Please make sure you have load `library(ggalluvial)` when running this function
        
        # dot plot
        netAnalysis_dot(cellchat, pattern = "incoming")
        
        P.incoming
        #graphics.off()
      dev.off()
   
  ##### Manifold and classification learning analysis of signaling networks #####

    #### Identify signaling groups based on their functional similarity ####
      cellchat <- computeNetSimilarity(cellchat, type = "functional")
      cellchat <- netEmbedding(cellchat, type = "functional")
      #> Manifold learning of the signaling networks for a single dataset
      cellchat <- netClustering(cellchat, type = "functional")
      #> Classification learning of the signaling networks for a single dataset
      # Visualization in 2D-space
      netVisual_embedding(cellchat, type = "functional", label.size = 3.5)

    #### Identify signaling groups based on structure similarity  ####
      cellchat <- computeNetSimilarity(cellchat, type = "structural")
      cellchat <- netEmbedding(cellchat, type = "structural")
      #> Manifold learning of the signaling networks for a single dataset
      cellchat <- netClustering(cellchat, type = "structural")
      #> Classification learning of the signaling networks for a single dataset
      # Visualization in 2D-space
      netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
      
    pdf(file = paste0(PathSys,"/",ProjectName,"_SystemsAnalysis_Classification.pdf"),
        width = 7,  height = 7
    )  
      netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
      netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)
      netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
      netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)
    dev.off()
  
##### Part V: Save the CellChat object #####    
    saveRDS(cellchat, file = paste0(Save.Path,"/",Version,".rds")) 
    #### Automatically save the plots of the all inferred network for quick exploration ####
    # # Access all the signaling pathways showing significant communications
    # pathways.show.all <- cellchat@netP$pathways
    # # check the order of cell identity to set suitable vertex.receiver
    # levels(cellchat@idents)
    # vertex.receiver = seq(1,4)
    # for (i in 1:length(pathways.show.all)) {
    #   # Visualize communication network associated with both signaling pathway and individual L-R pairs
    #   netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
    #   # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
    #   gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
    #   ggsave(filename=paste0(Version,"/",pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
    # }
      
      
##### Save CellChatDataBase #####
    PathDB <- paste0(Save.Path,"/DataBase")
    ## Create new folder
    if (!dir.exists(PathDB)){
      dir.create(PathDB)
    }
    
    #### Export Database Category ####
    pdf(file = paste0(PathDB,"/",ProjectName,"_CellChatDB.pdf"),
        width = 7,  height = 7
    )
      showDatabaseCategory(CellChatDB)
    dev.off()
    
    #### Export all database ####
    DB_Interact_All.df <- data.frame(Term = row.names(CellChatDB[["interaction"]]), CellChatDB[["interaction"]])
    write.table(DB_Interact_All.df, 
                file=paste0(PathDB,"/",ProjectName,"_DBAll_Interact.tsv"),sep="\t",
                row.names=F, quote = FALSE)
    DB_Complex_All.df <- data.frame(Term = row.names(CellChatDB[["complex"]]), CellChatDB[["complex"]])
    write.table(DB_Complex_All.df, 
                file=paste0(PathDB,"/",ProjectName,"_DBAll_Complex.tsv"),sep="\t",
                row.names=F, quote = FALSE)
    DB_Cofactor_All.df <- data.frame(Term = row.names(CellChatDB[["cofactor"]]), CellChatDB[["cofactor"]])
    write.table(DB_Cofactor_All.df, 
                file=paste0(PathDB,"/",ProjectName,"_DBAll_Cofactor.tsv"),sep="\t",
                row.names=F, quote = FALSE)
    DB_GeneInfo_All.df <- data.frame(Term = row.names(CellChatDB[["geneInfo"]]), CellChatDB[["geneInfo"]])
    write.table(DB_GeneInfo_All.df, 
                file=paste0(PathDB,"/",ProjectName,"_DBAll_GeneInfo.tsv"),sep="\t",
                row.names=F, quote = FALSE)
  
    #### Export used database ####
    DB_Interact.df <- data.frame(Term = row.names(CellChatDB.use[["interaction"]]), CellChatDB.use[["interaction"]])
    write.table(DB_Interact.df, 
                file=paste0(PathDB,"/",ProjectName,"_DB_Interact.tsv"),sep="\t",
                row.names=F, quote = FALSE)
    DB_Complex.df <- data.frame(Term = row.names(CellChatDB.use[["complex"]]), CellChatDB.use[["complex"]])
    write.table(DB_Complex.df, 
                file=paste0(PathDB,"/",ProjectName,"_DB_Complex.tsv"),sep="\t",
                row.names=F, quote = FALSE)
    DB_Cofactor.df <- data.frame(Term = row.names(CellChatDB.use[["cofactor"]]),CellChatDB.use[["cofactor"]])
    write.table(DB_Cofactor.df, 
                file=paste0(PathDB,"/",ProjectName,"_DB_Cofactor.tsv"),sep="\t",
                row.names=F, quote = FALSE)
    DB_GeneInfo.df <- data.frame(Term = row.names(CellChatDB.use[["geneInfo"]]),CellChatDB.use[["geneInfo"]])
    write.table(DB_GeneInfo.df, 
                file=paste0(PathDB,"/",ProjectName,"_DB_GeneInfo.tsv"),sep="\t",
                row.names=F, quote = FALSE)

    #### Catch the significant path ####
      DB_Interact_Sig.df <- DB_Interact.df[DB_Interact.df$pathway_name %in% pathway.set,]
      DB_Interact_Sig.df$pathway_name %>% unique()
      
      write.table(DB_Interact.df, 
                  file=paste0(PathDB,"/",ProjectName,"_DBSig_Interact.tsv"),sep="\t",
                  row.names=F, quote = FALSE)
      
  rm(list = str_subset(objects(), pattern = "DB_"))
#### Save the RData ####
  rm(list=setdiff(ls(), c("cellchat","CellChatDB","P.incoming","P.outgoing",
                          "pathway.set","seuratObject","Save.Path","Version","mat")))
  save.image(paste0(Save.Path,"/",Version,".RData")) 
  
  cellchatDB <- CellChatDB
  rm(list=setdiff(ls(), str_subset(objects(), pattern = "cellchat")))
 
      