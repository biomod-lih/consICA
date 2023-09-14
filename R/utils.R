#' @title Is the object is consensus ICA compliant?
#' @description Check if the object is a list in the same format as the result 
#' of `consICA()`
#' @param cica list
#' @return TRUE or FALSE
#' @usage is.consICA(cica)
#' @examples 
#' # returns TRUE
#' is.consICA(list("ncomp" = 2, "nsples" = 2, "nfeatures" = 2, 
#'                 "S" = matrix(0,2,2),"M" = matrix(0,2,2)))
#' @export
is.consICA <- function(cica){
    if(!is.list(cica))  return(FALSE)
    if(is.null(cica))   return(FALSE)
    
    if(is.null(cica$ncomp)) return(FALSE)
    if(!is.numeric(cica$ncomp)) return(FALSE)
    if(cica$ncomp <= 1) return(FALSE)
    
    if(is.null(cica$nsamples)) return(FALSE)
    if(!is.numeric(cica$nsamples)) return(FALSE)
    if(cica$nsamples <= 1) return(FALSE)
    
    if(is.null(cica$nfeatures)) return(FALSE)
    if(!is.numeric(cica$nfeatures)) return(FALSE)
    if(cica$nfeatures <= 1) return(FALSE)
    
    if(is.null(cica$S)) return(FALSE)
    if(!is.matrix(cica$S)) return(FALSE)
    if(nrow(cica$S) != cica$nfeatures | ncol(cica$S) != cica$ncomp)return(FALSE)
    
    if(is.null(cica$M)) return(FALSE)
    if(!is.matrix(cica$M)) return(FALSE)
    if(nrow(cica$M) != cica$ncomp | ncol(cica$M) != cica$nsamples)return(FALSE) 
    
    return(TRUE)
}

#' @title Sort dataframe
#' @description Sort dataframe, adapted from  
#' http://snippets.dzone.com/user/r-fanatic
#' @param x a data.frame
#' @param key sort by this column 
#' @param ... other parameters for `order` function (e. g. `decreasing`)
#' @usage sortDataFrame(x,key, ...)
#' @return sorted dataframe
#' @examples
#' df <- data.frame("features" = c("f1", "f2", "f3"), fdr = c(0.02, 0.002, 1))
#' sortDataFrame(df, "fdr")
#' @export
sortDataFrame <- function(x, key, ...) {
    if (missing(key)) {
        rn <- rownames(x)
        if (all(rn %in% seq.int(1,nrow(x)))) rn <- as.numeric(rn)
        x[order(rn, ...), , drop=FALSE]
    } else {
        x[do.call("order", c(x[key], ...)), , drop=FALSE]
    }
}

#' @title Sort Genes of consICA object
#' @description Sort Genes for independent components
#' @param Genes list compilant to `getFeatures`output
#' @return sorted list
#' @examples
#' #features <- list("ic1" = list(
#' #              "pos" = data.frame("features" = c("f1", "f2", "f3"),
#' #                                 "fdr" = c(0.0043, 0.4, 0.04)), 
#' #              "neg" = data.frame("features" = c("f1", "f2", "f3"),
#' #                                 "fdr" = c(0, 0.1, 0.9))))
#' #sortFeatures(features)
#' @keywords internal
sortFeatures <- function(Genes){
  if(length(Genes) < 1){
    message("Empty Genes vector\n")
    return(NULL)
  }
  ncomp <- length(Genes)
  for (icomp in seq.int(1,ncomp)){
    for (direct in c("neg","pos"))
      if (nrow(Genes[[icomp]][[direct]])>1)
        Genes[[icomp]][[direct]]<-sortDataFrame(Genes[[icomp]][[direct]],"fdr")
  }
  return(Genes)
}

getTopIdx <- function(x,n){
  return(order(x,na.last=TRUE,decreasing=TRUE)[seq.int(1,n)])
}

## Draw a table to graphical window
drawTable <- function(data,x0=0,y0=1,dx=0.2,dy=0.05,row.names=TRUE,
                      cex=1,col=1,new=TRUE, bg=NA, float.format = "%.2e"){
  if (is.matrix(data)) data <- data.frame(data,check.names = FALSE,
                                          stringsAsFactors = FALSE)
  ## factor->character
  i <- vapply(data, is.factor, FUN.VALUE = logical(1))
  data[i] <- lapply(data[i], as.character)  
  ## Make character table Tab from data 
  Tab <- data.frame(matrix(nrow = nrow(data)+1, 
                        ncol = ifelse(row.names,ncol(data)+1,ncol(data))))
  if (row.names) { # -/mc
    Tab[-1,1] <- rownames(data)
    Tab[1,1] <- ""
    Tab[1,-1]<-names(data) #
    Tab[-1,-1]<-data
    if (!is.null(float.format)){
      for (i in seq.int(1,ncol(data))){
        if (is.numeric(data[[i]])){
          if (is.data.frame(Tab[-1,-1])){
            Tab[-1,-1][,i] <- sprintf(float.format,data[[i]])
          }else{
            Tab[-1,-1] <- sprintf(float.format,data[[i]])
          }
        }
      }
    }
  } else {
    Tab[1,]<-names(data)
    Tab[-1,]<-data
    if (!is.null(float.format))
      for (i in seq.int(1,ncol(data)))
        if (is.numeric(data[[i]]))
          Tab[-1,][,i] <- sprintf(float.format,data[[i]])
  }
  
  x <- double(ncol(Tab))
  y <- double(nrow(Tab))
  ## check and elongate dx (if needed)
  if (length(dx) < length(x[-1])) 
    x <- c(0,dx[-length(dx)],rep(dx[length(dx)],ncol(Tab)-length(dx)))
  if (length(dx) > length(x[-1])) 
    x <- c(0,dx[seq.int(1,(length(x)-1))])
  if (length(dx) == length(x[-1])) 
    x <- c(0,dx)
  ## check and elongate dy (if needed)
  if (length(dy) < nrow(Tab)-1) 
    y <- c(0,dy[-length(dy)],rep(dy[length(dy)],nrow(Tab)-length(dy)))
  if (length(dy) > length(y[-1])) 
    y[-1] <- dy[seq.int(1,length(y[-1]))]
  if (length(dy) == length(y[-1]))
    y[-1] <- dy
  ## check and elongate (in 2 dim) colors
  if (length(col) < ncol(Tab)*nrow(Tab)) {
    col <- rep(col,ncol(Tab)*nrow(Tab))
  }
  
  x <- cumsum(x) + x0
  y <- y0 - cumsum(y)
  
  if (!new) 
    plot.new()
  k<-1;iy<-1;ix<-1
  for (iy in seq.int(1,nrow(Tab))){
    if (y[iy]<mean(dy)){
      y[seq.int(iy,nrow(Tab))] <- y[seq.int(iy,nrow(Tab))] + 1 - y[iy]
      plot.new()
    }
    for (ix in seq.int(1,ncol(Tab))){
      if (!is.na(bg)) rect(x[ix],y[iy], 1.05,y[iy]-mean(dy),col=bg,border=NA)
      text(x[ix],y[iy],Tab[iy,ix],adj=c(0,1),cex=cex,col=col[k],
           font = ifelse((ix==1 && row.names)|(iy==1),2,1))
      k<-k+1
    }
  }
}

violinplot<-function (x, ...,xlab="",ylab="", main="", range = 1.5, h = NULL, 
                      ylim = NULL, names = NULL, 
                      horizontal = FALSE, col = "magenta", 
                      border = "black", lty = 1, 
                      lwd = 1, rectCol = "black", colMed = "white", 
                      pchMed = 19, 
                      at, add = FALSE, wex = 1, drawRect = TRUE, cex.axis = 1){
  if(!is.list(x)){
    datas <- list(x, ...)
  } else{
    datas<-x
  }
  n <- length(datas)
  if (missing(at)) 
    at <- seq.int(1,n)
  upper <- vector(mode = "numeric", length = n)
  lower <- vector(mode = "numeric", length = n)
  q1 <- vector(mode = "numeric", length = n)
  q3 <- vector(mode = "numeric", length = n)
  med <- vector(mode = "numeric", length = n)
  base <- vector(mode = "list", length = n)
  height <- vector(mode = "list", length = n)
  baserange <- c(Inf, -Inf)
  args <- list(display = "none")
  if (!(is.null(h))) 
    args <- c(args, h = h)
  for (i in seq.int(1,n)) {
    data <- datas[[i]]
    data.min <- min(data)
    data.max <- max(data)
    q1[i] <- quantile(data, 0.25)
    q3[i] <- quantile(data, 0.75)
    med[i] <- median(data)
    iqd <- q3[i] - q1[i]
    upper[i] <- min(q3[i] + range * iqd, data.max)
    lower[i] <- max(q1[i] - range * iqd, data.min)
    est.xlim <- c(min(lower[i], data.min), max(upper[i], 
                                               data.max))
    smout <- do.call("sm.density", c(list(data, xlim = est.xlim), 
                                     args))
    hscale <- 0.4/max(smout$estimate) * wex
    base[[i]] <- smout$eval.points
    height[[i]] <- smout$estimate * hscale
    t <- range(base[[i]])
    baserange[1] <- min(baserange[1], t[1])
    baserange[2] <- max(baserange[2], t[2])
  }
  if (!add) {
    xlim <- if (n == 1) 
      at + c(-0.5, 0.5)
    else range(at) + min(diff(at))/2 * c(-1, 1)
    if (is.null(ylim)) {
      ylim <- baserange
    }
  }
  if (is.null(names)) {
    label <- names(x)
  }
  else {
    label <- names
  }
  boxwidth <- 0.05 * wex
  if (!add) 
    plot.new()
  if (!horizontal) {
    if (!add) {
      plot.window(xlim = xlim, ylim = ylim)
      axis(2,cex.axis =cex.axis )
      axis(1, at = at, labels = label, las = 2,cex.axis=cex.axis  )
    }
    box()
    for (i in seq.int(1,n)) {
      polygon(c(at[i] - height[[i]], rev(at[i] + height[[i]])), 
              c(base[[i]], rev(base[[i]])), col = col, border = border, 
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(at[c(i, i)], c(lower[i], upper[i]), lwd = lwd, 
              lty = lty)
        rect(at[i] - boxwidth/2, q1[i], at[i] + boxwidth/2, 
             q3[i], col = rectCol)
        points(at[i], med[i], pch = pchMed, col = colMed)
      }
    }
  }
  else {
    if (!add) {
      plot.window(xlim = ylim, ylim = xlim)
      axis(1,cex.axis=cex.axis  )
      axis(2, at = at, labels = label,cex.axis=cex.axis )
    }
    box()
    for (i in seq.int(1,n)) {
      polygon(c(base[[i]], rev(base[[i]])), c(at[i] - height[[i]], 
                                              rev(at[i] + height[[i]])), 
              col = col, border = border, 
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(c(lower[i], upper[i]), at[c(i, i)], lwd = lwd, 
              lty = lty)
        rect(q1[i], at[i] - boxwidth/2, q3[i], at[i] + 
               boxwidth/2, col = rectCol)
        points(med[i], at[i], pch = pchMed, col = colMed)
      }
    }
  }
  invisible(list(upper = upper, lower = lower, median = med, 
                 q1 = q1, q3 = q3))
  title(main,xlab=xlab,ylab=ylab)
}

num2fact <- function(x, nlev = 4, digits=1){
  if (!is.numeric(x)) return(factor(x))
  f <- rep(ifelse(digits==1,sprintf("q%d",1),sprintf("q%02d",1)),length(x))
  thr <- quantile(x,seq(1/nlev,1,by = 1/nlev),na.rm=TRUE)
  for (ilev in seq.int(2,nlev))
    f[x > thr[ilev-1]] <- ifelse(digits==1,sprintf("q%d",ilev),
                                sprintf("q%02d",ilev))
  f[is.na(x)] <- NA
  f <- factor(f)
  return(f)
}

#' Set up for the parallel computing for biocParallel
#' Adapt from `FEAST`
#' This function sets up the environment for parallel computing.
#' @param ncores number of processors
#' @param BPPARAM bpparameter from bpparam
#' @return BAPPARAM settings
#' @keywords internal
set_bpparam <- function(ncores = 0, BPPARAM = NULL){
  if (is.null(BPPARAM)) {
    if (ncores != 0) {
      if (.Platform$OS.type == "windows") {
        result <- SnowParam(workers = ncores)
      }
      else {
        result <- MulticoreParam(workers = ncores)
      }
    }
    else {
      result <- bpparam()
    }
    return(result)
  }
  else {
    return(BPPARAM)
  }
}

#' @title Convert input object as numeric matrix
#' @param obj input data with features in rows and samples in columns. 
#' Could be a `SummarizedExperiment` object, matrix or `Seurat` object. 
#' For `SummarizedExperiment` with multiple assays or `Seurat` pass the name 
#' with `assay_string` parameter, otherwise the first will be taken.
#' See \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
#' @param assay_string name of assay for `SummarizedExperiment` or `Seurat` 
#' input object `obj`. Default value is NULL
#' @return matrix
#' @keywords internal
get_X_num <- function(obj, assay_string = NULL){
  
  if(! (inherits (obj, "SummarizedExperiment") | inherits (obj, "matrix") | 
        inherits (obj, "Seurat")) ){
    message("Input must be a matrix, SummarizedExperiment or Seurat object")
    return(NULL)
  }
  
  if(inherits (obj, "SummarizedExperiment")){
    
    if(is.null(assay_string)){
      X_num <- as.matrix(assay(obj))
    } else {
      if(length(which(names(SummarizedExperiment::assays(obj)) == assay_string)) == 0){
        message(
          "Given name of assay was not found in X. Check `assay_string`")
        return(NULL)
      }
      X_num <- as.matrix(SummarizedExperiment::assays(obj)[[assay_string]])
    }
  }
  
  if(inherits (obj, "matrix")){
    X_num <- obj
  }

  if(inherits (obj, "Seurat")){
    suppressPackageStartupMessages({
      requireNamespace("Seurat")
    })
    
    if(is.null(assay_string)) assay_string <- "data"
    
    if(length(which(c("data","counts","scale.data") == assay_string)) == 0){
      message(
        "Given name of assay was not found in X. Check `assay_string`")
      #return(NULL)
    }
    X_num <- as.matrix(GetAssayData(object = obj, slot = assay_string))
  }
  return(X_num)
}
