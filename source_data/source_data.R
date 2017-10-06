
source("source_files/additional_functions1.R")
source("source_files/additional_objects1.R")
options(warn=-1)

#loading simulation results
##Load Zhang simulation data and source count matrix
result.Zhang <- readRDS("source_files/results.All.Zhang.RData")
Zhang.source <- readRDS("source_files/Zhang_Counts.RData")
Zhang.DESeq2.results <- readRDS("source_files/Zhang_full_DESeq2_analysis_result.RData")


##Load NGP Nutlin simulation data
result.NGP <- readRDS("source_files/results.All.NGP.RData")
NGP.source <- readRDS("source_files/celine_neuroblastoma_data.RData")
NGP.DESeq2.results <- readRDS("source_files/NGP_full_DESeq2_analysis_result.RData")


##Load GTEx simulation data
result.GTEx <- readRDS("source_files/results.All.GTEx.RData")
GTEx.source <- readRDS("source_files/GTEx_data_full.RData")
GTEx.DESeq2.results <- readRDS("source_files/GTEx_full_DESeq2_analysis_result.RData")



#Loading table of methods
methods.table <- read.table("source_files/methods summary text.txt", header = T, sep="\t")
methods.table2 <- as.data.frame(apply(methods.table, 2, function(x) return(as.character(x))))
colnames(methods.table2) <- c("DE.tool.full", "DE.tool", "refTool", "countAssum", "normMethod", "description",
                              "numFactors", "outAdjust", "package", "packageVersion", "refPackage")




  
