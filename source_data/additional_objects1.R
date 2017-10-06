methods.names <- c("edgeR.exact", "edgeR.GLM", "edgeR.robust", "edgeR.QL", "DESeq", "DESeq2", "limmaQN",
                   "limmaVoom", "limmaVoom_QW", "limmaVst", "PoissonSeq", "SAMSeq", "QuasiSeq")
methods.names2 <- c("edgeR exact", "edgeR GLM", "edgeR robust", "edgeR QL", "DESeq", "DESeq2", "limmaQN",
                    "limmaVoom", "limmaVoom+QW", "limmaVst", "PoissonSeq", "SAMSeq", "QuasiSeq")
cols <-c("slateblue4", "palegreen1", "palegreen3", "palegreen4", 
         "skyblue2","royalblue4",
         "goldenrod1", "salmon1", "salmon4", "lightpink3", 
         "purple1", 
         "gray41", 
         "mediumvioletred")

cols.biotype <- list(mRNA="deepskyblue3", lncRNA="orange")

ggplot2.theme <- theme(axis.text.x=element_text(size = 15, colour = "gray30"), 
                       axis.text.y=element_text(size = 15, colour = "gray30"),
                       axis.title = element_text(size = 17, colour = "black"),
                       plot.title =element_text(size=18, face="bold"),
                       panel.background = element_rect(fill = "gray90"),
                       panel.grid.major = element_line(colour = "white"),
                       strip.background = element_rect(colour = "white", fill = "gray60"),
                       strip.text = element_text(size = 17, colour = "black"),
                       panel.grid.major.y = element_blank(),
                       legend.text= element_text(size=15),
                       legend.key.size  = unit(2,"line"),
                       legend.key = element_rect(colour = "transparent", fill = "white"),
                       legend.title = element_blank(),
                       legend.position = c(0.5, 0.8),
                       legend.box = "horizontal",
                       legend.justification = "left")

