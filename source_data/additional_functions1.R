#Filtrarttion 
filter.count    <- function(counts, group, method="minimum read counts", min.read=1, 
                            normalized.lib.sizes=TRUE, min.cpm = 0.001, ...){
  if(is.null(rownames(counts))) stop("The count matrix has no unique row names.")
  group = as.factor(group)
  if(method == "minimum read counts"){
    counts.f <- counts[which(rowSums(counts[, group == levels(group)[1]])>=min.read &
                               rowSums(counts[, group == levels(group)[2]])>=min.read), ]
    return(list(counts=counts.f, group=group, method=method, min.read=min.read,  
                number.of.filtered.tags = nrow(counts)- nrow(counts.f)))
  }
  else if(method == "CPM"){
    #require(edgeR)
    counts.cpm <- cpm(counts, normalized.lib.size=snormalized.lib.sizes)
    counts.f   <- counts.cpm[which(rowMeans(counts.cpm[, group == levels(group)[1]])>=min.cpm &
                                     rowMeans(counts.cpm[, group == levels(group)[2]])>=min.cpm), ]
    counts.f   <- counts[rownames(counts.f), ]
    return(list(counts=counts.f, group=group, method=method, min.cpm=min.cpm, 
                number.of.filtered.tags = nrow(counts)- nrow(counts.f)))
  }
} 

#Normalization function
normalize.count <- function(counts, group, norm.method = "DESeq", ...){
  group <- as.factor(group)
  if(norm.method=="QN"){
    #Quantile Normalization
    #require(preprocessCore)
    norm.counts <- normalize.quantiles(as.matrix(counts))
    
    colnames(norm.counts) <- colnames(counts)
    rownames(norm.counts) <- rownames(counts)
    return(norm.counts)
  }
  else if(norm.method=="limmaQN"){
    #Quantile Normalization
    #require(limma)
    counts.log.dat=log2(counts+0.5)
    norm.counts=normalizeBetweenArrays(counts.log.dat,method='quantile')
    
    colnames(norm.counts) <- colnames(counts)
    return(norm.counts)
  }
  
  else if(norm.method=="TMM"){
    
    #TMM Normalization
    #require(edgeR)
    y <- DGEList(counts=counts, group=group)
    y <- calcNormFactors(y)
    y <- estimateCommonDisp(y) 
    
    norm.counts <- y$pseudo.counts
    
    #LS <- apply(counts, 2, sum)
    #v = calcNormFactors(y, method="TMM")$samples[,3]*LS
    
    #norm.counts <- matrix(NA, ncol=ncol(counts), nrow=nrow(counts))
    #for(i in 1:ncol(counts)){
    #  norm.counts[, i] <- counts[,i]/v[i]*1e6
    #}
    
    colnames(norm.counts) <- colnames(counts)    
    return(norm.counts)
  }
  
  #DESeq Normalization
  else if(norm.method=="DESeq"){
    #require(DESeq)
    cds <- newCountDataSet(counts, group)
    ## estimate size factor
    cds <- estimateSizeFactors(cds)
    norm.counts <- counts(cds, normalized=TRUE)
    
    colnames(norm.counts) <- colnames(counts)
    return(norm.counts)
  }
  
  else if(norm.method=="PoissonSeq"){
    #require(PoissonSeq)
    
    seq.depth <- PS.Est.Depth(counts, ct.sum=1, ct.mean=0)
    norm.counts <- counts
    for(i in 1:ncol(counts)){
      norm.counts[,i] <- counts[,i]/seq.depth[i]
    }
    
    colnames(norm.counts) <- colnames(counts)
    return(norm.counts)
  }
  
  
  else if(norm.method=="SAMSeq"){
    #require(samr)
    
    x <- counts
    y <- ifelse(group==levels(group)[1],1,2)  
    samfit <- SAMseq(x, y, geneid=rownames(counts), resp.type="Two class unpaired", fdr.output=1.0)
    
    ls <-apply(counts, 2, sum)
    depth <- samfit$samr.obj$depth
    require(psych)
    gm <- geometric.mean(depth)
    factor <- diag(gm/depth)
    
    norm.counts <- as.matrix(x)%*%factor
    colnames(norm.counts) <- colnames(counts)
    return(norm.counts)
  }
}

#Correlation summary
cor.summary <- function(counts, group, prob=c(0,0.25, 0.5, 0.75, 1)){
  group <- as.factor(group)
  counts <- counts[, c(which(group == levels(group)[1]), 
                       which(group == levels(group)[2]))]
  
  WC1 <- cor(counts[, which(group == levels(group)[1])])
  WC1 <- WC1[lower.tri(WC1, diag = FALSE)]
  WC2 <- cor(counts[, which(group == levels(group)[2])])
  WC2 <- WC2[lower.tri(WC2, diag = FALSE)]
  
  WC.summary <- quantile(c(WC1, WC2), prob=prob)
  
  BC <- cor(counts)[1:length(which(group == levels(group)[1])),
                    (length(which(group == levels(group)[1]))+1):length(group)]
  BC <- BC[lower.tri(BC, diag = FALSE)]
  BC.summary <- quantile(BC, prob=prob)
  
  return(list(within.cor=WC.summary, between.cor =BC.summary ))
}


attributes.names <- function(st){
  if(st=="DE.tool")        return("DE tool [short name]")
  else if(st=="DE.tool.full")   return("DE tool")
  else if(st=="package")        return("R Package name")
  else if(st=="description")    return("Description")
  else if(st=="normMethod")     return("Normalization method")
  else if(st=="countAssum")     return("Count distribution assumption")
  else if(st=="refTool")        return("DE tool reference")
  else if(st=="numFactors")     return("Possible number of factors supported")
  else if(st=="outAdjust")      return("Integrated outlier treatment method")
  else if(st=="refPackage")     return("R package reference")
  else if(st=="packageVersion") return("R package version")
}



convert.name <- function(x){
  x <- as.character(x)
  y=x
  if(x=="edgeR.exact") {y="edgeR exact"}
  else if(x=="edgeR.GLM") {y="edgeR GLM"}
  else if(x=="edgeR.robust") {y="edgeR robust"}
  else if(x=="edgeR.QL") {y="edgeR QL"}
  else if(x=="limmaVoom_QW") {y="limmaVoom+QW"}
  else if(x=="limmaVoom") {y="limmaVoom"}
  else if(x=="DESeq") {y="DESeq"}
  else if(x=="DESeq2") {y="DESeq2"}
  else if(x=="QuasiSeq") {y="QuasiSeq"}
  else if(x=="limmaVst") {y="limmaVst"}
  else if(x=="limmaQN") {y="limmaQN"}
  else if(x=="PoissonSeq") {y="PoissonSeq"}
  else if(x=="SAMSeq") {y="SAMSeq"}
  return(as.character(y))
}

convert.name2 <- function(x){
  x <- as.character(x)
  y=x
  if(x=="edgeR exact") {y="edgeR.exact"}
  else if(x=="edgeR GLM") {y="edgeR.GLM"}
  else if(x=="edgeR robust") {y="edgeR.robust"}
  else if(x=="edgeR QL") {y="edgeR.QL"}
  else if(x=="limmaVoom+QW") {y="limmaVoom_QW"}
  else if(x=="limmaVoom") {y="limmaVoom"}
  else if(x=="DESeq") {y="DESeq"}
  else if(x=="DESeq2") {y="DESeq2"}
  else if(x=="QuasiSeq") {y="QuasiSeq"}
  else if(x=="limmaVst") {y="limmaVst"}
  else if(x=="limmaQN") {y="limmaQN"}
  else if(x=="PoissonSeq") {y="PoissonSeq"}
  else if(x=="SAMSeq") {y="SAMSeq"}
  return(as.character(y))
}

#class(convert.name("edgeR ex"))

DE.tools.codes <- function(x){
  if(x=="edgeR exact")        {fun=run_edgeR_exact}
  else if(x=="edgeR GLM")     {fun=run_edgeR_glm}
  else if(x=="edgeR robust")  {fun=run_edgeR_robust}
  else if(x=="edgeR QL")      {fun=run_edgeR_ql}
  else if(x=="limmaVoom+QW")  {fun=run_limmaVoom_QW}
  else if(x=="limmaVoom")     {fun=run_limmaVoom}
  else if(x=="DESeq")         {fun=run_DESeq}
  else if(x=="DESeq2")        {fun=run_DESeq2}
  else if(x=="QuasiSeq")      {fun=run_QuasiSeq}
  else if(x=="limmaVst")      {fun=run_limmaVst}
  else if(x=="limmaQN")       {fun=run_limmaQN}
  else if(x=="PoissonSeq")    {fun=run_PoissonSeq}
  else if(x=="SAMSeq")        {fun=run_SAMSeq}
  format(fun)
}

