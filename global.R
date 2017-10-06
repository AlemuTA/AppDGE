#Load required packages
library(shiny)        #CRAN
library(shinyjs)      #CRAN
library(shinyBS)      #CRAN
library(shinythemes)  #CRAN 

library(reshape2)      #CRAN 
library(ggplot2)       #CRAN 
library(data.table)    #CRAN 
library(plyr)          #CRAN 
library(gridExtra)     #CRAN 
library(pheatmap)      #CRAN

library(edgeR)         #Bioconductor
library(DESeq)        #Bioconductor
library(DESeq2)        #Bioconductor
library(limma)         #Bioconductor
library(PoissonSeq)        #CRAN 
library(samr)              #CRAN 
library(preprocessCore)    #Bioconductor


#Load source files
source("source_files/source_data.R")
source("ui.R")
source("server.R")

shinyApp(ui = ui, server = server)
