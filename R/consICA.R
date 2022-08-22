#' @title Calculate consensus Independent Component Analysis 
#' @description calculate consensus independent component analysis (ICA) 
#' @param X a `SummarizedExperiment` object. Assay used as data matrix with 
#' features in rows and samples in columns. 
#' See \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
#' @param ncomp number of components. 
#' @param ntry number of consensus runs. Default value is 1
#' @param show.every numeric logging period in iterations (disabled for 
#' `ncore`s > 1). Default value is 1
#' @param filter.thr Filter out genes (rows) with values lower than this value 
#' from `X` 
#' @param ncores number of cores to be set for parallel calculation. Default 
#' value is 1
#' @param bpparam parameters from the `BiocParallel`
#' @param reduced If TRUE returns reduced result (no `X`, `i.best`, 
#' see 'return')
#' @param exclude are samples excluded during multiple run? If TRUE one random 
#' sample will be excluded per run. Default is TRUE 
#' @param fun the functional form of the G function used in the approximation 
#' to neg-entropy in fastICA. Default value is "logcosh"
#' @param alg.typ parameter for fastICA(). If alg.typ == "deflation" the 
#' components are extracted one at a time. If alg.typ == "parallel" the 
#' components are extracted simultaneously. Default value is "deflation"
#' @param verbose logic TRUE or FALSE. Use TRUE for print process steps. 
#'     Default value is FALSE
#' @return a list with
#'         \item{X}{input `SummarizedExperiment` object}
#'         \item{nsamples, nfeatures}{dimension of X}
#'         \item{S, M}{consensus metagene and weight matrix}
#'         \item{ncomp}{number of components}
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

#' @importFrom fastICA fastICA
#' @importFrom BiocParallel bplapply SnowParam MulticoreParam bpparam
#' @import topGO org.Hs.eg.db 
#' @import GO.db 
#' @importFrom graph nodeData
#' @importFrom sm sm.density
#' @importFrom graphics axis barplot box boxplot legend lines par plot.new 
#' plot.window points polygon rect text title
#' @importFrom grDevices dev.off pdf
#' @importFrom stats aov cor density mad median p.adjust pnorm quantile rexp
#' @importFrom SummarizedExperiment assay colData
#' @importFrom methods new

consICA <- function(X, 
                    ncomp=10,
                    ntry=1,
                    show.every=1,
                    filter.thr=NULL, 
                    ncores=1,
                    bpparam=NULL,
                    reduced=FALSE,
                    exclude=TRUE,
                    fun="logcosh",
                    alg.typ="deflation",
                    verbose=FALSE){
  
    if(!inherits (X, "SummarizedExperiment")){
      message("X must be a SummarizedExperiment object")
      return(NULL)
    }
  
    Xse <- X
    X <- as.matrix(assay(X))
    if (!is.null(filter.thr)) X <- X[apply(X,1,max)>filter.thr, , drop=FALSE]

    Res <- list() # output
    Res$X <- Xse
    S <- list()
    M <- list()
    
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
    idx.excludeSamples <- sample(seq.int(1,ncol(X)), 
                                 ntry, replace = (ntry > ncol(X)))
    if (ntry==1 | (!exclude) ) idx.excludeSamples <- integer(0)
    itry <- 1
    
    ###############################
    ## Parallel section starts
    ##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    if(verbose){
      message("*** Starting ",ifelse(ncores>1,"parallel ",""),"calculation on ",
          ncores,"core(s)...\n")
      message("*** System: ",.Platform$OS.type,"\n")
      message("***",ncomp," components, ",ntry," runs, ",nrow(X)," features, ",
          ncol(X),"samples.\n")
      message("*** Start time: ",as.character(Sys.time()),"\n")
    }
    t0 <- Sys.time()
    ## multi run ICA
    if (ncores > 1) {
      
      par_fica <- function(x=x, Res=Res,Z,ncomp=ncomp,alg.typ=alg.typ,fun=fun,
                           idx.excludeSamples=idx.excludeSamples){
        suppressPackageStartupMessages({
          requireNamespace("fastICA") 
        })

        SP <- Res$S + NA
        MP <- Res$M + NA

        if (length(idx.excludeSamples) == 0){
          ic <- fastICA(Z,n.comp=ncomp,alg.typ=alg.typ,fun=fun)
          SP[,] <- ic$S
          MP[,] <- ic$A
        }else{
          z <- Z[,-idx.excludeSamples[x]]
          ic <- fastICA(z,n.comp=ncomp,alg.typ=alg.typ,fun=fun)
          SP[,] <- ic$S
          MP[,-idx.excludeSamples[x]] <- ic$A
          MP[is.na(MP)] <- 0
        }
        return(list(S=SP,M=MP))
      }
      
      bp_param <- set_bpparam(ncores, BPPARAM = bpparam)
      bp_param$progressbar = TRUE
      seqntry <- seq.int(ntry)
      MRICA <- bplapply(X=seqntry, FUN = par_fica, BPPARAM = bp_param,
                        Res=Res,Z=X,ncomp=ncomp,alg.typ=alg.typ,fun=fun,
                        idx.excludeSamples=idx.excludeSamples)

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
            if (length(idx.excludeSamples) == 0){
                ic <- fastICA(X, n.comp=ncomp, alg.typ=alg.typ,fun=fun)
                MRICA[[itry]]$S[,] <- ic$S
                MRICA[[itry]]$M[,] <- ic$A
            }else{
                x <- X[,-idx.excludeSamples[itry]]
                ic <- fastICA(x, n.comp=ncomp, alg.typ=alg.typ,fun=fun)
                MRICA[[itry]]$S[,] <- ic$S
                MRICA[[itry]]$M[,-idx.excludeSamples[itry]] <- ic$A
                MRICA[[itry]]$M[is.na(MRICA[[itry]]$M)] <- 0
            }
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
        Res$i.best <- NULL
    }
    
    Res$ncomp  <- ncomp
    Res$nsamples <- ncol(X)
    Res$nfeatures <- nrow(X)
    if(verbose) message("consensus ICA done\n")
    return(Res)
}

#' @title Runs \code{\link[fastICA]{fastICA}}
#' @description Runs \code{\link[fastICA]{fastICA}} once and store in a 
#' consICA manner
#' @param X a `SummarizedExperiment` object. Assay used as data matrix with 
#' features in rows and samples in columns. 
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
                  alg.typ="deflation"){ 
  
    if(!inherits (X, "SummarizedExperiment")){
      message("X must be a SummarizedExperiment object")
      return(NULL)
    }
  
    Xse <- X
    X <- as.matrix(assay(X))    
    Res <- list()
    
    ## create containers
    X <- as.matrix(X)
    if (!is.null(filter.thr)) X <- X[apply(X,1,max)>filter.thr, ,drop=FALSE]
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
    
    return(Res)
}
