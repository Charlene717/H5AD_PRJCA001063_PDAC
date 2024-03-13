library(Seurat)

# Function to read marker genes from file and add module scores
AddModuleScoreFromFile <- function(markerFileName, markerName, seuratObject) {
  markerFilePath <- paste0(getwd(), "/", markerFileName, ".txt")
  markerGenes <- read.delim(markerFilePath, header = TRUE, stringsAsFactors = FALSE)
  
  # Assuming the first column contains gene names
  markerGenesList <- as.list(markerGenes[[1]])
  
  # Adding module score
  seuratObject <- AddModuleScore(
    object = seuratObject,
    features = list(markerGenesList),
    name = markerName
  )
  
  return(seuratObject)
}

# Function to plot cells colored by a given module score
PlotCellsByModule <- function(cds, moduleName, lowColor = "darkblue", midColor = "#f7c211", highColor = "green", midpoint = 0.15) {
  plot_cells(cds, color_cells_by = moduleName, label_cell_groups = FALSE, show_trajectory_graph = FALSE, cell_size = 1.2) +
    scale_colour_gradient2(low = lowColor, mid = midColor, high = highColor, guide = "colourbar", midpoint = midpoint, labs(fill = moduleName))
}

# Marker file names and their corresponding module names
markers <- list(
  PDAC = "GRUETZMANN_PANCREATIC_CANCER_UP",
  EMT = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  ACST = "HP_ABNORMALITY_OF_CHROMOSOME_STABILITY",
  Migration = "GOBP_POSITIVE_REGULATION_OF_EPITHELIAL_CELL_MIGRATION",
  Metastasis = "NAKAMURA_METASTASIS_MODEL_UP",
  NE = "REACTOME_INITIATION_OF_NUCLEAR_ENVELOPE_NE_REFORMATION",
  NP = "REACTOME_NUCLEAR_PORE_COMPLEX_NPC_DISASSEMBLY",
  ATR = "REACTOME_ACTIVATION_OF_ATR_IN_RESPONSE_TO_REPLICATION_STRESS",
  G0G1 = "REACTOME_G0_AND_EARLY_G1"
)

# Loop through each marker and add module scores
for (markerName in names(markers)) {
  markerFileName <- markers[[markerName]]
  cds <- AddModuleScoreFromFile(markerFileName, markerName, cds)
  PlotCellsByModule(cds, markerName)
}



