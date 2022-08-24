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
#' @param id id
#' @param algorithm algorithm for `runTest()`
#' @param do.sort if TRUE - resulted functions sorted by p-value 
#' @param randomFraction for testing only, the fraction of the genes to be 
#' randomized
#' @param return.genes If TRUE include genes in output. Default value is FALSE
#' @return list with terms and stats
#' @author Petr V. Nazarov
enrichGO <- function(genes,
                    fdr=NULL,fc=NULL,ntop=NA,thr.fdr=0.05,thr.fc=NA,
                    db="BP",
                    genome="org.Hs.eg.db", 
                    id = c("entrez", "alias", "ensembl", "symbol", 
                           "genename"),
                    algorithm = "weight",
                    do.sort=TRUE,randomFraction=0,return.genes=FALSE){

  groupGOTerms() 
  if (is.na(ntop)) message("ntop is NA => using FDR (FC) as limits\n")
  
  ## create score depending on threshold and paradigm
  score <- get_score(genes, fc, thr.fc, fdr, thr.fdr, ntop)
  
  message("Significant:",sum(score>0),"\n")
  
  ## prepare gene categories and score
  if (length(id)>1){ 
    ## if so - find the best, i.e. the one which is overrepresented
    common <- vapply(id, 
                      function(x_id){
                        myGO2genes <- annFUN.org(whichOnto=db,mapping = genome,
                                                 ID = x_id)
                        sum(genes %in% unique(unlist(myGO2genes)))
                      },
                      FUN.VALUE = numeric(1))
    
    id <- names(which.max(common))
    message("Checking DB. The best overlap with [",id,"] gene IDs:",
            max(common),"\n")
  }
  
  myGO2genes <- annFUN.org(whichOnto=db, mapping = genome, ID = id)

  ## create topGOdata object
  GOdata = new("topGOdata",  ##constructor
               ontology = db,
               allGenes = score,
               geneSelectionFun = SelectScore,
               annot = annFUN.GO2genes,
               GO2genes = myGO2genes)
  ## run testing
  resFisher <- runTest(GOdata, algorithm = algorithm, statistic = "fisher")
  ## transform results into a table
  enrichRes <- GenTable(GOdata, classicFisher = resFisher, 
                       ranksOf = "classicFisher",
                       topNodes = length(score(resFisher)))
  enrichRes$classicFisher[grep("<",enrichRes$classicFisher)] <- "1e-31"
  enrichRes$classicFisher <- as.double(enrichRes$classicFisher)
  enrichRes$FDR <- p.adjust(enrichRes$classicFisher,"fdr")
  enrichRes$Score <- -log10(enrichRes$FDR)
  ## by default sorted by p-value. If needed - sort by ID
  if (!do.sort) enrichRes <- sortDataFrame(enrichRes,"GO.ID") ## remove sorting
  rownames(enrichRes) <- enrichRes$GO.ID
  
  graph_data <- nodeData(graph(GOdata))
  
  if (return.genes){
    sig.genes <- genes[score>0]
    g2g <- matrix("",nrow = length(names(graph_data)),ncol=2)
    rownames(g2g) <-  names(graph_data)
    colnames(g2g) <- c("genes.all","genes.sig")
    for (i in seq.int(1,nrow(g2g))) {
      x <- ls(graph_data[[i]][[1]])
      x.sig <- x[x%in%sig.genes]
      g2g[i,"genes.all"] <- paste(x,collapse=",")
      g2g[i,"genes.sig"] <- paste(x.sig,collapse=",")
    }
    enrichRes$genes.sig <- g2g[enrichRes$GO.ID,"genes.sig"]
    enrichRes$genes.all <- g2g[enrichRes$GO.ID,"genes.all"]
  }
  return(enrichRes)
}

#' @title Create score depending on threshold and paradigm
#' @param genes character vector with list of ENSEBML IDs 
#' @param fdr numeric vector of FDR for each gene
#' @param fc numeric vector of logFC for each gene
#' @param ntop number of first taken genes
#' @param thr.fdr significance threshold for FDR 
#' @param thr.fc significance threshold for absolute logFC
#' @return numeric score vector
get_score <- function(genes, fc, thr.fc, fdr, thr.fdr, ntop){
  
  ## create score depending on threshold and paradigm
  if (is.null(fc) | is.na(thr.fc)){
    score <- (-log10(fdr))
    if (is.na(ntop)){
      score[fdr>=thr.fdr]<-0
    }else{
      score[sort(score,index.return=TRUE,decreasing=TRUE)$ix[
        -(seq_along(ntop))]] <- 0
    }
  }else{
    score <- (-log10(fdr)*abs(fc))
    if (is.na(ntop)){
      score[fdr>=thr.fdr | abs(fc)<=thr.fc] <- 0
    }else{
      score[sort(score,index.return=TRUE,decreasing=TRUE)$ix[
        -seq_along(ntop)]] <- 0
    }
  }
  names(score) <- genes
  return(score)
}

# Simple function for significance
SelectScore <- function(sc){return(sc>0)} 
