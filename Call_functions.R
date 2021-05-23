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

getFilePath("Call_functionsTTT.R")
print(paste("filePath:",filePath))
Convert.stringToInt("1")   #這個function存在於methods.R裡面
