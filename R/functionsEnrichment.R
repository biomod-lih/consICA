#' @title Enrichment analysis of GO-terms based on Ensembl IDs 
#' @description Enrichment analysis of GO-terms for independent components with 
#' Ensembl IDs based on topGO package
#' @param genes character vector with list of ENSEBML IDs 
#' @param fdr numeric vector of FDR for each gene
#' @param fc numeric vector of logFC for each gene
#' @param ntop number of first taken genes
#' @param thr.fdr significance threshold for FDR 
#' @param thr.fc significance threshold for absolute logFC
#' @param db name of GO database: "BP","MF","CC"
#' @param genome R-package for genome annotation used. For human -
#' 'org.Hs.eg.db'
#' @param do.sort if TRUE - resulted functions sorted by p-value 
#' @param randomFraction for testing only, the fraction of the genes to be 
#' randomized
#' @param return.genes If TRUE include genes in output. Default value is FALSE
#' @return list with terms and stats
#' @author Petr V. Nazarov
enrichGOens = function(genes,
                       fdr=NULL,
                       fc=NULL,
                       ntop=NA,
                       thr.fdr=0.05,
                       thr.fc=NA,
                       db="BP",
                       genome="org.Hs.eg.db",
                       do.sort=TRUE,
                       randomFraction=0,
                       return.genes=FALSE){
  if (is.na(ntop)) cat("ntop is NA => using FDR (FC) as limits\n")

  ## create score depending on threshold and paradigm
  if (is.null(fc) | is.na(thr.fc)){
    score = (-log10(fdr))
    if (is.na(ntop)){
      score[fdr>=thr.fdr]=0
    }else{
      score[sort(score,index.return=TRUE,decreasing=TRUE)$ix[-(1:ntop)]]=0
    }
  }else{
    score = (-log10(fdr)*abs(fc))
    if (is.na(ntop)){
      score[fdr>=thr.fdr | abs(fc)<=thr.fc]=0
    }else{
      score[sort(score,index.return=TRUE,decreasing=TRUE)$ix[-(1:ntop)]]=0
    }
  }
  names(score)=genes
  
  ## add randomness if required, to test stability
  if (randomFraction>0){ 
    ## define remove probability: low score have more chances
    prob1 = 1/(1+score)
    prob1[is.na(prob1)]=0
    prob1[score == 0] = 0
    ## define add probability: high score has more chances
    prob2 = -log10(fdr)*abs(fc)
    prob2[is.na(prob2)]=0
    prob2[score > 0] = 0
    ## add and remove
    n=round(sum(score>0)*randomFraction)
    score[sample(1:length(genes),n,prob=prob2)]=1+rexp(n,1/mean(score[score>0]))
    score[sample(1:length(genes),n,prob=prob1)]=0
  }
  
  ## prepare gene categories and score
  myGO2genes <- topGO::annFUN.org(whichOnto=db, mapping = genome, ID="ensembl")

  ## create topGOdata object
  SelectScore = function(sc){return(sc>0)} ## simple function for significance
  GOdata = new("topGOdata",  ##constructor
               ontology = db,
               allGenes = score,
               geneSelectionFun = SelectScore,
               annot = annFUN.GO2genes,
               GO2genes = myGO2genes)
  #annot = annFUN.org,
  #mapping = genome) 
  #ID = "ensembl")
  ## run testing
  resFisher = runTest(GOdata, algorithm = "classic", statistic = "fisher")
  ## transform results into a table
  enrichRes = GenTable(GOdata, classicFisher = resFisher, 
                       ranksOf = "classicFisher",
                       topNodes = length(resFisher@score))
  enrichRes$classicFisher[grep("<",enrichRes$classicFisher)] = "1e-31"
  enrichRes$classicFisher = as.double(enrichRes$classicFisher)
  enrichRes$FDR = p.adjust(enrichRes$classicFisher,"fdr")
  enrichRes$Score = -log10(enrichRes$FDR)
  ## by default sorted by p-value. If needed - sort by ID
  if (!do.sort) enrichRes = sortDataFrame(enrichRes,"GO.ID") ## remove sorting
  rownames(enrichRes) = enrichRes$GO.ID
  
  if (return.genes){
    sig.genes = genes[score>0]
    g2g = matrix("",nrow = length(names(GOdata@graph@nodeData@data)),ncol=4)
    rownames(g2g) =  names(GOdata@graph@nodeData@data)
    colnames(g2g) = c("n.sig","n.all","genes.sig","genes.all")
    for (i in 1:nrow(g2g)) {
      x = ls(GOdata@graph@nodeData@data[[i]][[1]])
      x.sig = x[x%in%sig.genes]
      g2g[i,"genes.all"] = paste(x,collapse=",")
      g2g[i,"genes.sig"] = paste(x.sig,collapse=",")
      g2g[i,"n.all"] = length(x)
      g2g[i,"n.sig"] = length(x.sig)
      if(i%%1000 == 0) print(i)
    }
    enrichRes$n.sig = as.integer(g2g[enrichRes$GO.ID,"n.sig"])
    enrichRes$n.all = as.integer(g2g[enrichRes$GO.ID,"n.all"])
    enrichRes$genes.sig = g2g[enrichRes$GO.ID,"genes.sig"]
    enrichRes$genes.all = g2g[enrichRes$GO.ID,"genes.all"]
  }
  return(enrichRes)
}

## Gene Symbol
enrichGO = function(genes,
                    fdr=NULL,fc=NULL,ntop=NA,thr.fdr=0.05,thr.fc=NA,
                    db="BP",
                    genome="org.Hs.eg.db", 
                    id = c("entrez", "alias", "ensembl", "symbol", 
                           "genename", "unigene"),
                    algorithm = "weight",
                    do.sort=TRUE,randomFraction=0,return.genes=FALSE){
  #requireNamespace("topGO")
  groupGOTerms() 
  if (is.na(ntop)) cat("ntop is NA => using FDR (FC) as limits\n")
  ## create score depending on threshold and paradigm
  if (is.null(fc) | is.na(thr.fc)){
    score = (-log10(fdr))
    if (is.na(ntop)){
      score[fdr>=thr.fdr]=0
    }else{
      score[sort(score,index.return=TRUE,decreasing=TRUE)$ix[-(1:ntop)]]=0
    }
  }else{
    score = (-log10(fdr)*abs(fc))
    if (is.na(ntop)){
      score[fdr>=thr.fdr | abs(fc)<=thr.fc]=0
    }else{
      score[sort(score,index.return=TRUE,decreasing=TRUE)$ix[-(1:ntop)]]=0
    }
  }
  names(score)=genes
  
  cat("Significant:",sum(score>0),"\n")
  ## prepare gene categories and score
  
  if (length(id)>1){ 
    ## if so - find the best, i.e. the one which is overrepresented
    common = integer(length(id))
    names(common) = id
    for ( i in 1:length(id)){
      myGO2genes <- topGO::annFUN.org(whichOnto=db,mapping = genome,ID = id[i])
      common[i] = sum(genes %in% unique(unlist(myGO2genes)))
    }
    print(common)
    id = names(which.max(common))
  cat("Checking DB. The best overlap with [",id,"] gene IDs:",max(common),"\n")
  }
  myGO2genes <- topGO::annFUN.org(whichOnto=db, mapping = genome, ID = id)
  
  ## create topGOdata object
  SelectScore = function(sc){return(sc>0)} ## simple function for significance
  GOdata = new("topGOdata",  ##constructor
               ontology = db,
               allGenes = score,
               geneSelectionFun = SelectScore,
               annot = annFUN.GO2genes,
               GO2genes = myGO2genes)
  ## run testing
  resFisher = runTest(GOdata, algorithm = algorithm, statistic = "fisher")
  ## transform results into a table
  enrichRes = GenTable(GOdata, classicFisher = resFisher, 
                       ranksOf = "classicFisher",
                       topNodes = length(resFisher@score))
  enrichRes$classicFisher[grep("<",enrichRes$classicFisher)] = "1e-31"
  enrichRes$classicFisher = as.double(enrichRes$classicFisher)
  enrichRes$FDR = p.adjust(enrichRes$classicFisher,"fdr")
  enrichRes$Score = -log10(enrichRes$FDR)
  ## by default sorted by p-value. If needed - sort by ID
  if (!do.sort) enrichRes = sortDataFrame(enrichRes,"GO.ID") ## remove sorting
  rownames(enrichRes) = enrichRes$GO.ID
  
  if (return.genes){
    sig.genes = genes[score>0]
    g2g = matrix("",nrow = length(names(GOdata@graph@nodeData@data)),ncol=2)
    rownames(g2g) =  names(GOdata@graph@nodeData@data)
    colnames(g2g) = c("genes.all","genes.sig")
    for (i in 1:nrow(g2g)) {
      x = ls(GOdata@graph@nodeData@data[[i]][[1]])
      x.sig = x[x%in%sig.genes]
      g2g[i,"genes.all"] = paste(x,collapse=",")
      g2g[i,"genes.sig"] = paste(x.sig,collapse=",")
      if(i%%1000 == 0) print(i)
    }
    enrichRes$genes.sig = g2g[enrichRes$GO.ID,"genes.sig"]
    enrichRes$genes.all = g2g[enrichRes$GO.ID,"genes.all"]
  }
  return(enrichRes)
}
