#' @title Assigns IC signatures to Gene Ontologies
#' @description Assigns extracted independent components to Gene Ontologies and
#' rotate independent components (`S` matrix) to set most significant 
#' Gene Ontologies as positive affecting features. Set `ncores` param for 
#' paralleled calculations.
#' @param cica  list compliant to `consICA()` result
#' @param alpha value in [0,1] interval. Used to filter features with 
#' FDR < `alpha`. Default value is 0.05
#' @param genenames alternative names of genes. If NULL we use rownames of `S`
#' matrix. We automatically identify type of gene identifier, you can use  
#' Ensembl, Symbol, Entrez, Alias, Genename IDs.
#' @param genome R-package for genome annotation used. For human -
#' 'org.Hs.eg.db'
#' @param db name of GO database: "BP","MF","CC"
#' @param ncores number of cores for parallel calculation. Default 
#' value is 4
#' @param rotate rotate components in `S` and `M` matricies in `cica` object 
#' to set most significant Gene Ontologies as positive effective features. 
#' Default is TRUE
#' @return rotated (if need) `cica` object with added `GO` - 
#' list for each db chosen (BP, CC, MM), 
#' with dataframes `pos` for 
#' positive and `neg` for negative affecting features for each component:
#'     \item{GO.ID}{id of Gene Ontology term}
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
#' # Calculate ICA (run with ntry=1 for quick test, use more in real analysis)
#' cica <- consICA(samples_data, ncomp=2, ntry=1, ncores=1, show.every=0) 
#' # cica <- consICA(samples_data, ncomp=40, ntry=20, show.every=0)
#' 
#' # Annotate independent components with gene ontoligies
#' cica <- getGO(cica, db = "BP", ncores=1)
#' # Positively affected GOs for 2nd independent component
#' head(cica$GO$GOBP$ic02$pos)
#' @export

getGO <- function(cica,
                  alpha=0.05, 
                  genenames=NULL, 
                  genome="org.Hs.eg.db",
                  db=c("BP", "CC", "MF"),
                  ncores = 4,
                  rotate=TRUE){
  
  if(!is.consICA(cica)) {
    message("First parameter should be compliant to `consICA()` result\n")
    return (NULL)
  }
  IC <- cica
  
  is_bp <- "BP" %in% db
  is_cc <- "CC" %in% db
  is_mf <- "MF" %in% db
  
  genes_id_type <- c("entrez", "ensembl", "symbol", "genename")

  if(!(is_bp | is_cc | is_mf)){
    message("Check db parameter!\n")
    return(NULL)
  }
  
  par_getGO<-function(X,IC=IC,genenames,is_bp,is_cc,is_mf,genes_id_type,genome,
                        enrichGO=enrichGO,
                        get_score=get_score, SelectScore=SelectScore){
    suppressPackageStartupMessages({
      requireNamespace("topGO")
      requireNamespace("org.Hs.eg.db")
      requireNamespace("GO.db")
    })
    
    GOBP <- list()
    GOCC <- list()
    GOMF <- list()
    
    icomp <- X
    
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

    GOBP[[sprintf("ic%02d",icomp)]] <- list()
    GOCC[[sprintf("ic%02d",icomp)]] <- list()
    GOMF[[sprintf("ic%02d",icomp)]] <- list()
    message("---- Component ",icomp," ----\n")
    
    message(str(GOBP))
    if(is_bp)GOBP[[1]]$pos <- enrichGO(genes = genes,
                                       fdr = fdr.pos,thr.fdr=alpha,db="BP",
                                       id = genes_id_type,
                                       genome=genome)
    if(is_cc)GOCC[[1]]$pos <- enrichGO(genes = genes,
                                       fdr = fdr.pos,thr.fdr=alpha,db="CC",
                                       id= genes_id_type,
                                       genome=genome)
    if(is_mf)GOMF[[1]]$pos <- enrichGO(genes = genes,
                                       fdr = fdr.pos,thr.fdr=alpha,db="MF",
                                       id= genes_id_type,
                                       genome=genome)

    if(is_bp)GOBP[[1]]$neg <- enrichGO(genes = genes,
                                       fdr = fdr.neg,thr.fdr=alpha,db="BP",
                                       id= genes_id_type,
                                       genome=genome)
    if(is_cc)GOCC[[1]]$neg <- enrichGO(genes = genes,
                                       fdr = fdr.neg,thr.fdr=alpha,db="CC",
                                       id= genes_id_type,
                                       genome=genome)
    if(is_mf)GOMF[[1]]$neg <- enrichGO(genes = genes,
                                       fdr = fdr.neg,thr.fdr=alpha,db="MF",
                                       id= genes_id_type,
                                       genome=genome)
    
    return(list(GOBP=GOBP, GOCC=GOCC,GOMF=GOMF))
  }
  
  bpparam <-NULL
  bp_param <- set_bpparam(ncores, BPPARAM = bpparam)
  bp_param$progressbar <- TRUE
  seqicomp <- seq.int(1,ncol(IC$S))
  RESGO <- bplapply(X=seqicomp, FUN = par_getGO, BPPARAM = bp_param,
                    IC=IC,genenames=genenames,is_bp=is_bp, is_cc=is_cc, 
                    is_mf=is_mf,genes_id_type=genes_id_type,genome=genome, 
                    enrichGO=enrichGO, get_score=get_score, 
                    SelectScore=SelectScore)
  
  GO <- list()
  if(is_bp) {
    la <- lapply(seqicomp, function(i) RESGO[[i]]$GOBP[[1]])
    names_ics <- sprintf("ic%02d",seqicomp)
    names(la) <- names_ics
    GO <- append(GO, list("GOBP"=la))
    names(GO) <- c(names(GO)[-length(names(GO))], "GOBP")
  }
  if(is_cc) {
    la <- lapply(seqicomp, function(i) RESGO[[i]]$GOCC[[1]])
    names_ics <- sprintf("ic%02d",seqicomp)
    names(la) <- names_ics
    GO <- append(GO, list("GOCC"=la))
    names(GO) <- c(names(GO)[-length(names(GO))], "GOCC")
  }
  if(is_mf) {
    la <- lapply(seqicomp, function(i) RESGO[[i]]$GOMF[[1]])
    names_ics <- sprintf("ic%02d",seqicomp)
    names(la) <- names_ics
    GO <- append(GO, list("GOMF"=la))
    names(GO) <- c(names(GO)[-length(names(GO))], "GOMF")
  }
  cica$GO <- GO
  if(rotate){
    cica <- setOrientation(cica)
  }
  rm(GO, RESGO)
  gc()
  return(cica)
}
