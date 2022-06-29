#' @title Assigns IC signatures to Gene Ontologies
#' @description Assigns extracted independent components to Gene Ontologies
#' @param IC  list compliant to `consICA()` result
#' @param alpha value in [0,1] interval. Used to filter features with 
#' FDR < `alpha`. Default value is 0.05
#' @param genenames alternative names of genes
#' @param genome R-package for genome annotation used. For human -
#' 'org.Hs.eg.db'
#' @param db name of GO database: "BP","MF","CC"
#' @return list for each db chosen (BP, CC, MM), with dataframes `pos` for 
#' positive and `neg` for negative affecting features for each component:
#'     \item{GO.ID}{id of gene ontology term}
#'     \item{Term}{name of term}
#'     \item{Annotated}{number of annotated genes}
#'     \item{Significant}{number of significant genes}
#'     \item{Expected}{estimate of the number of annotated genes if the 
#'     significant} 
#'     \item{genes would be randomly selected from the gene universe
#'     classisFisher}{F-test}
#'     \item{FDR}{false discovery rate value}
#'     \item{Score}{genes score}
#' @author Petr V. Nazarov
#' @examples
#' data("samples_data")
#' # cica <- consICA(samples_data, ncomp=40, ntry=1, show.every=0)
#' cica <- consICA(samples_data, ncomp=2, ntry=1, show.every=0) #exp timesave
#' GOs <- getGO(cica, db = "BP")
#' @export
getGO <- function(IC,
                 alpha=0.05, 
                 genenames=NULL, 
                 genome="org.Hs.eg.db",
                 db=c("BP", "CC", "MF")){
  
  is_bp <- "BP" %in% db
  is_cc <- "CC" %in% db
  is_mf <- "MF" %in% db
  
  Genes <- list()
  if(is_bp)GOBP <- list()
  if(is_cc)GOCC <- list()
  if(is_mf)GOMF <- list()

  if(!(is_bp | is_cc | is_mf)){
    message("Check db parameter!\n")
    return(NULL)
  }
  
  icomp <- 1
  for (icomp in seq.int(1,ncol(IC$S))){
    pv <- pnorm(IC$S[,icomp],mean=median(IC$S[,icomp]),sd=mad(IC$S[,icomp]))
    pv[pv>0.5] <- 1-pv[pv>0.5]
    fdr <- p.adjust(pv,method="BH")
    fdr.pos <- fdr
    fdr.neg <- fdr
    fdr.pos[fdr>alpha | scale(IC$S[,icomp])<0] <- 1
    fdr.neg[fdr>alpha | scale(IC$S[,icomp])>0] <- 1
    
    if (is.null(genenames)) {
      genes <- rownames(IC$S) 
    }else{
      if(!is.null(names(genenames))){
        genes <- genenames[rownames(IC$S)]
      }else{
        genes <- genenames
      }
    }
    
    if(is_bp)GOBP[[sprintf("ic%02d",icomp)]] <- list()
    if(is_cc)GOCC[[sprintf("ic%02d",icomp)]] <- list()
    if(is_mf)GOMF[[sprintf("ic%02d",icomp)]] <- list()
    message("---- Component ",icomp," ----\n\n")
    if(is_bp)GOBP[[icomp]]$pos <- enrichGO(genes = genes,
                                          fdr = fdr.pos,thr.fdr=alpha,db="BP",
                                          id= c("entrez", "ensembl", "symbol",
                                                "genename"),
                                          genome=genome)
    if(is_cc)GOCC[[icomp]]$pos <- enrichGO(genes = genes,
                                          fdr = fdr.pos,thr.fdr=alpha,db="CC",
                                          id= c("entrez", "ensembl", "symbol", 
                                                "genename"),
                                          genome=genome)
    if(is_mf)GOMF[[icomp]]$pos <- enrichGO(genes = genes,
                                          fdr = fdr.pos,thr.fdr=alpha,db="MF",
                                          id= c("entrez", "ensembl", "symbol", 
                                                 "genename"),
                                          genome=genome)

    if(is_bp)GOBP[[icomp]]$neg <- enrichGO(genes = genes,
                                          fdr = fdr.neg,thr.fdr=alpha,db="BP",
                                          id= c("entrez", "ensembl", "symbol", 
                                                "genename"),
                                          genome=genome)
    if(is_cc)GOCC[[icomp]]$neg <- enrichGO(genes = genes,
                                          fdr = fdr.neg,thr.fdr=alpha,db="CC",
                                          id= c("entrez", "ensembl", "symbol", 
                                                "genename"),
                                          genome=genome)
    if(is_mf)GOMF[[icomp]]$neg <- enrichGO(genes = genes,
                                          fdr = fdr.neg,thr.fdr=alpha,db="MF",
                                          id= c("entrez", "ensembl", "symbol", 
                                                "genename"),
                                          genome=genome)

  }
  GO <- list()
  if(is_bp) {
    GO <- append(GO, list(GOBP))
    names(GO) <- c(names(GO)[-length(names(GO))], "GOBP")
  }
  if(is_cc) {
    GO <- append(GO, list(GOCC))
    names(GO) <- c(names(GO)[-length(names(GO))], "GOCC")
  }
  if(is_mf) {
    GO <- append(GO, list(GOMF))
    names(GO) <- c(names(GO)[-length(names(GO))], "GOMF")
  }
  rm(GOBP, GOCC, GOMF)
  gc()
  return(GO)
}
