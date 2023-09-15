#' @title Calculate consensus Independent Component Analysis 
#' @description calculate consensus independent component analysis (ICA) 
#' Implements efficient ICA calculations.
#' @param X input data with features in rows and samples in columns. 
#' Could be a `SummarizedExperiment` object, matrix or `Seurat` object. 
#' For `SummarizedExperiment` with multiple assays or `Seurat` pass the name 
#' with `assay_string` parameter, otherwise the first will be taken.
#' See \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
#' @param ncomp number of components
#' @param ntry number of consensus runs. Default value is 1
#' @param show.every numeric logging period in iterations (disabled for 
#' `ncore`s > 1). Default value is 1
#' @param filter.thr Filter out genes (rows) with max value lower than this 
#' value from `X` 
#' @param ncores number of cores for parallel calculation. Default 
#' value is 4
#' @param bpparam parameters from the `BiocParallel`
#' @param reduced If TRUE returns reduced result (no `X`, `i.best`, 
#' see 'return')
#' @param fun the functional form of the G function used in the approximation 
#' to neg-entropy in fastICA. Default value is "logcosh"
#' @param alg.typ parameter for fastICA(). If alg.typ == "deflation" the 
#' components are extracted one at a time. If alg.typ == "parallel" the 
#' components are extracted simultaneously. Default value is "deflation"
#' @param verbose logic TRUE or FALSE. Use TRUE for print process steps. 
#'     Default value is FALSE
#' @param assay_string name of assay for `SummarizedExperiment` or `Seurat` 
#' input object `X`. Default value is NULL
#' @return a list with
#'         \item{X}{input object}
#'         \item{nsamples, nfeatures}{dimension of X}
#'         \item{S, M}{consensus metagene and weight matrix}
#'         \item{ncomp}{number of components}
#'         \item{X_num}{input data in matrix format}
#'         \item{mr2}{mean R2 between rows of M}
#'         \item{stab}{stability, mean R2 between consistent columns of S in 
#'         multiple tries. Applicable only for `ntry` > 1}
#'         \item{i.best}{number of best iteration}
#' @author Petr V. Nazarov
#' @examples 
#' data("samples_data")
#' # Deconvolve into independent components
#' cica <- consICA(samples_data, ncomp=15, ntry=10, show.every=0)
#' # X = S * M, where S - independent signals matrix, M - weights matrix
#' dim(samples_data)
#' dim(cica$S)
#' dim(cica$M)
#' @seealso \code{\link[fastICA]{fastICA}}
#' @export

#' @importFrom fastICA fastICA ica.R.def ica.R.par
#' @importFrom BiocParallel bplapply SnowParam MulticoreParam bpparam
#' @import topGO org.Hs.eg.db 
#' @import GO.db 
#' @importFrom graph nodeData
#' @importFrom sm sm.density
#' @importFrom graphics axis barplot box boxplot legend lines par plot.new 
#' plot.window points polygon rect text title
#' @importFrom grDevices dev.off pdf
#' @importFrom stats aov cor density mad median p.adjust pnorm quantile rexp 
#' rnorm
#' @importFrom SummarizedExperiment assay colData
#' @importFrom methods new
#' @importFrom Rfast mat.mult

consICA <- function(X, 
                    ncomp=10,
                    ntry=1,
                    show.every=1,
                    filter.thr=NULL, 
                    ncores=1,
                    bpparam=NULL,
                    reduced=FALSE,
                    fun="logcosh",
                    alg.typ="deflation",
                    verbose=FALSE,
                    assay_string = NULL){
  
    if(! (inherits (X, "SummarizedExperiment") | inherits (X, "matrix") | 
        inherits (X, "Seurat")) ){
      message("X must be a matrix, SummarizedExperiment or Seurat object")
      return(NULL)
    }
  
    if(inherits (X, "SummarizedExperiment")){
      
      if(is.null(assay_string)){
        if (!is.null(filter.thr)) {
          ind <- apply(assay(X),1,max)>filter.thr
          X <-  X[ind,]
        }
        Xse <- X
        X <- as.matrix(assay(X))
      } else {
        if(length(which(names(SummarizedExperiment::assays(X)) == 
                                                        assay_string)) == 0){
          message(
            "Given name of assay was not found in X. Check `assay_string`")
          return(NULL)
        }
        if (!is.null(filter.thr)) {
          ind <- apply(SummarizedExperiment::assays(X)[[assay_string]],1,max) >
                                                                    filter.thr
          X <-  X[ind,]
        }
        Xse <- X
        X <- as.matrix(SummarizedExperiment::assays(X)[[assay_string]])
      }
    }
  
    if(inherits (X, "matrix")){
      if (!is.null(filter.thr)) {
        ind <- apply(X,1,max)>filter.thr
        X <-  X[ind,]
      }
      Xse <- X
    }
  
    if(inherits (X, "Seurat")){
      suppressPackageStartupMessages({
        requireNamespace("Seurat")
      })
      
      if(is.null(assay_string)) assay_string <- "data"
      
      if(length(which(c("data","counts","scale.data") == assay_string)) == 0){
        message(
          "Given name of assay was not found in X. Check `assay_string`")
        #return(NULL)
      }
      
      if (!is.null(filter.thr)) {
        ind <- apply(as.matrix(
                          GetAssayData(object = X, slot = assay_string)
                      ),1,max)>filter.thr
        X <-  subset(X, features = rownames(X)[ind])
      }
      Xse <- X
      X <- as.matrix(GetAssayData(object = X, slot = assay_string))
    }
  
    if(length(X) == 0){
      message("After filtering input matrix is empty. Check `filter.thr`\n")
      return(NULL)
    }
  
    Res <- list() # output
    Res$X <- Xse
    S <- list()
    M <- list()
    Res$X_num <- X
    
    ## metagenes  - matrix of features
    S[[1]] <- matrix(nrow=nrow(X),ncol=ncomp)
    rownames(S[[1]]) <- rownames(X)
    colnames(S[[1]]) <- sprintf("ic.%d",seq.int(1,ncomp))
    
    ## mixing matrix - matrix of samples
    M[[1]] <- matrix(nrow=ncomp,ncol=ncol(X))
    colnames(M[[1]]) <- colnames(X)
    rownames(M[[1]]) <- sprintf("ic.%d",seq.int(1,ncomp))
    
    itry <- 1
    Res$S <- S[[1]]
    Res$M <- M[[1]]
    Res$mr2 <- NA  ## mean correlation bw mixing profiles
   
     ## do multiple tries
    itry <- 1
    
    ###############################
    ## Parallel section starts
    ##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    preICA <- outICA(X, n.comp=ncomp, verbose = verbose)
    
    if(verbose){
      message("*** Starting ",ifelse(ncores>1,"parallel ",""),"calculation on ",
          ncores,"core(s)...\n")
      message("*** System: ",.Platform$OS.type,"\n")
      message("*** ",ncomp," components, ",ntry," runs, ",nrow(X)," features, ",
          ncol(X),"samples.\n")
      message("*** Start time: ",as.character(Sys.time()),"\n")
    }
    t0 <- Sys.time()
    ## multi run ICA
    if (ncores > 1) {
      
      par_fica <- function(x=x, Res=Res,Z,ncomp=ncomp,alg.typ=alg.typ,fun=fun){
        suppressPackageStartupMessages({
          requireNamespace("fastICA")
          requireNamespace("Rfast")
        })

        SP <- Res$S + NA
        MP <- Res$M + NA

        #ic <- fastICA(Z,n.comp=ncomp,alg.typ=alg.typ,fun=fun)
        ic <- coreICA(Z, n.comp=ncomp,preICA=preICA, alg.typ=alg.typ,fun=fun)
        SP[,] <- ic$S
        MP[,] <- ic$A
        
        return(list(S=SP,M=MP))
      }
      
      bp_param <- set_bpparam(ncores, BPPARAM = bpparam)
      bp_param$progressbar <- TRUE
      seqntry <- seq.int(ntry)
      MRICA <- bplapply(X=seqntry, FUN = par_fica, BPPARAM = bp_param,
                        Res=Res,Z=X,ncomp=ncomp,alg.typ=alg.typ,fun=fun)

    }else{
      if(verbose){
        m <- ifelse(show.every > 0, 
                     paste("showing progress every",show.every,"run(s)\n"),
                     "\n")
        message("Execute one-core analysis ", m)
      }
        
        MRICA <- list()
        for(itry in seq.int(1,ntry)){
            MRICA[[itry]] <- list()
            MRICA[[itry]]$S <- Res$S + NA
            MRICA[[itry]]$M <- Res$M + NA
            
            #ic <- fastICA(X, n.comp=ncomp, alg.typ=alg.typ,fun=fun)
            ic <- coreICA(X, n.comp=ncomp,preICA=preICA,alg.typ=alg.typ,fun=fun)
            MRICA[[itry]]$S[,] <- ic$S
            MRICA[[itry]]$M[,] <- ic$A
            
            if (itry%%show.every == 0 & show.every > 0) {
              if(verbose) message("try #",itry," of ",ntry,"\n")
            }
        }
    }
    
    if(verbose){
      message("*** Done!", "\n")
      message("*** End time:",as.character(Sys.time()),"\n")
      message(Sys.time()-t0, "\n")
    }
   
    ## end of parallel section
    ###############################
    
    S <- lapply(MRICA,function(x)x$S)
    M <- lapply(MRICA,function(x)x$M)
    rm(MRICA)
    
    if(verbose) {
      message("Calculate ||X-SxM|| and r2 between component weights\n")
    }
    
    Res$mr2 <- vapply(seq.int(1,ntry), 
                       function(itry){
                         mean((cor(t(M[[itry]]))^2)
                              [upper.tri(matrix(0,nrow=ncomp,ncol=ncomp))])
                       }, numeric(1))
    
    ## who is the best: min correlated (hoping that we capture majority of 
    ## independent signals
    Res$i.best <- which.min(Res$mr2)
    ## simply smallest error?
    if (length(Res$i.best) == 0 ) Res$i.best <- 1
    ## if only one try - return
    if (ntry == 1) {
        Res$S <- S[[1]]
        Res$M <- M[[1]]
        Res$ncomp  <- ncomp
        Res$nsamples <- ncol(X)
        Res$nfeatures <- nrow(X)
        return(Res)
    }
    
    ## correlate results
    ## s.cor - to which ic of the BEST decomposition we should address?
    if(verbose) message("Correlate rows of S between tries\n")
    s.cor <- matrix(nrow=ntry,ncol=ncomp)
    s.cor[Res$i.best,] <- seq.int(1,ncomp)
    itry <- 1
    for (itry in seq.int(1,ntry)[-Res$i.best]) {
        r <- cor(S[[itry]],S[[Res$i.best]])
        s.cor[itry,] <- apply((r)^2,2,which.max)
        for (ic in seq.int(1,ncomp))
            s.cor[itry,ic] <- s.cor[itry,ic] * sign(r[s.cor[itry,ic],ic])
    }
    ## build consensus S, M
    
    if(verbose) message("Build consensus ICA\n")
    Res$S[,] <- S[[1]]
    Res$M[,] <- M[[1]]
    itry<-2 # delete
    for (itry in seq.int(2,ntry)){ ## itry=2, because 1 is already there
        for (ic in seq.int(1,ncomp)) {
            Res$S[,ic] <- Res$S[,ic] + S[[itry]][,abs(s.cor[itry,ic])]* 
                sign(s.cor[itry,ic])
            Res$M[ic,] <- Res$M[ic,] + M[[itry]][abs(s.cor[itry,ic]),]* 
                sign(s.cor[itry,ic])
        }
    }
    Res$S <- Res$S / ntry
    Res$M <- Res$M / ntry
    
    ## use consensus S, M to analyze stability
    if(verbose) message("Analyse stability\n")
    Res$stab <- t(vapply(seq.int(1,ntry), 
                      function(itry){
                        diag(cor(Res$S,S[[itry]][,abs(s.cor[itry,])])^2)
                      }, numeric(ncomp)))
    colnames(Res$stab) <- colnames(Res$S)
    
    if (reduced) {
        Res$X <- NULL
        Res$X_num <- NULL
        Res$i.best <- NULL
    }
    
    Res$ncomp  <- ncomp
    Res$nsamples <- ncol(X)
    Res$nfeatures <- nrow(X)
    if(!is.null(assay_string)) Res$assay_name <- assay_string
    
    if(verbose) message("consensus ICA done\n")
    return(Res)
}

#' @title Runs \code{\link[fastICA]{fastICA}}
#' @description Runs \code{\link[fastICA]{fastICA}} once and store in a 
#' consICA manner
#' @param X input data with features in rows and samples in columns. 
#' Could be a `SummarizedExperiment` object, matrix or `Seurat` object. 
#' For `SummarizedExperiment` with multiple assays or `Seurat` pass the name 
#' with `assay_string` parameter, otherwise the first will be taken.
#' See \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
#' @param ncomp number of components. Default value is 10
#' @param filter.thr filter rows in input matrix with max value > `filter.thr`.
#' Default value is NULL
#' @param reduced If TRUE returns reduced result (no X, see 'return')
#' @param fun the functional form of the G function used in the approximation 
#' to neg-entropy in fastICA. Default value is "logcosh"
#' @param alg.typ parameter for fastICA(). if alg.typ == "deflation" the 
#' components are extracted one 
#' at a time. if alg.typ == "parallel" the components are extracted 
#' simultaneously. Default value is "deflation"
#' @param assay_string name of assay for `SummarizedExperiment` or `Seurat` 
#' input object `X`. Default value is NULL
#' @return a list with
#'         \item{X}{input `SummarizedExperiment` object}
#'         \item{nsamples, nfeatures}{dimension of X assay}
#'         \item{S, M}{consensus metagene and weight matrix}
#'         \item{ncomp}{number of components}
#' @author Petr V. Nazarov
#' @examples 
#' data("samples_data")
#' res <- oneICA(samples_data)
#' @seealso \code{\link[fastICA]{fastICA}}
#' @export
oneICA <- function(X,
                  ncomp=10, 
                  filter.thr=NULL,
                  reduced=FALSE, 
                  fun="logcosh", 
                  alg.typ="deflation",
                  assay_string=NULL){ 
  
    if(! (inherits (X, "SummarizedExperiment") | inherits (X, "matrix") | 
          inherits (X, "Seurat")) ){
      message("X must be a matrix, SummarizedExperiment or Seurat object")
      return(NULL)
    }
  
    ## create containers
    # if (!is.null(filter.thr)) {
    #   ind <- apply(assay(X),1,max)>filter.thr
    #   X <-  X[ind,]
    # }
    # Xse <- X
    # X <- as.matrix(assay(X))    
    
    if(inherits (X, "SummarizedExperiment")){
      
      if(is.null(assay_string)){
        if (!is.null(filter.thr)) {
          ind <- apply(assay(X),1,max)>filter.thr
          X <-  X[ind,]
        }
        Xse <- X
        X <- as.matrix(assay(X))
      } else {
        if(length(which(names(SummarizedExperiment::assays(X)) == assay_string)) 
                                                                          == 0){
          message(
            "Given name of assay was not found in X. Check `assay_string`")
          return(NULL)
        }
        if (!is.null(filter.thr)) {
          ind <- apply(SummarizedExperiment::assays(X)[[assay_string]],1,max)>
                                                                      filter.thr
          X <-  X[ind,]
        }
        Xse <- X
        X <- as.matrix(SummarizedExperiment::assays(X)[[assay_string]])
      }
    }
    
    if(inherits (X, "matrix")){
      if (!is.null(filter.thr)) {
        ind <- apply(X,1,max)>filter.thr
        X <-  X[ind,]
      }
      Xse <- X
    }
    
    if(inherits (X, "Seurat")){
      suppressPackageStartupMessages({
        requireNamespace("Seurat")
      })
      
      if(is.null(assay_string)) assay_string <- "data"
      
      if(length(which(c("data","counts","scale.data") == assay_string)) == 0){
        message(
          "Given name of assay was not found in X. Check `assay_string`")
        #return(NULL)
      }
      
      if (!is.null(filter.thr)) {
        ind <- apply(as.matrix(
          GetAssayData(object = X, slot = assay_string)
        ),1,max)>filter.thr
        X <-  subset(X, features = rownames(X)[ind])
      }
      Xse <- X
      X <- as.matrix(GetAssayData(object = X, slot = assay_string))
    }
    
    Res <- list()
    if (!reduced) Res$X <- Xse
    
    ## matrix of signal
    Res$S <- matrix(nrow=nrow(X),ncol=ncomp)
    rownames(Res$S) <- rownames(X)
    colnames(Res$S) <- sprintf("ic.%d", seq.int(1,ncomp))
    
    ## mixing matrix
    Res$M <- matrix(nrow=ncomp,ncol=ncol(X))
    colnames(Res$M) <- colnames(Res$X)
    rownames(Res$M) <- sprintf("ic.%d", seq.int(1,ncomp))
    
    ## run fastICA
    ic <- fastICA(X, n.comp=ncomp, alg.typ=alg.typ, fun=fun)
    
    ## prepare data to be returned
    Res$S[,] <- ic$S
    Res$M[,] <- ic$A
    colnames(Res$M) <- colnames(X)
    
    Res$ncomp <- ncomp
    Res$nsamples <- ncol(X)
    Res$nfeatures <- nrow(X)
    if(!is.null(assay_string)) Res$assay_name <- assay_string
    
    return(Res)
}
