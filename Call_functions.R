# https://dotblogs.com.tw/TingI/2018/12/02/230722

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

## Convert Monocle3 Object to Seurat Object
getFilePath("Monocle3_To_Seurat.R")
print(paste("filePath:",filePath))
marrow <- Monocle3_To_Seurat(cds) #這個function存在於Monocle3_To_Seurat.R裡面

## Assign Cell-Cycle Scores
getFilePath("Cell-Cycle Scoring and Regression.R")
marrow <- CCScorReg(GeneNAFMT,marrow,colors_cc,Main) #這個function存在於Cell-Cycle Scoring and Regression.R裡面

## PlotViolin 
getFilePath("PlotViolin.R")
TTT <- PlotViolin (cds,colors_cc,Main)
