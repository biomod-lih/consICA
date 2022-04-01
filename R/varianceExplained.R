#' @title Estimate the variance explained by the model
#' @description This is an adaptation of 
#' \code{\link[MOFA2]{calculate_variance_explained}}.\cr The method estimates 
#' the variance explained by the model and by each independent component.\cr
#' As in the original \code{\link[MOFA2]{calculate_variance_explained}}, 
#' we used the coefficient\cr of determination (R2) between the normalized 
#' input (X-mean(X)) and (S*M)
#' @param cica list compliant to `consICA()` result
#' @param X input matrix used for the model. Will be used if consICA$X is NULL,
#' ignore otherwise.
#' @return a list of:\cr
#'     \item{R2}{total variance explained by the model}
#'     \item{R2_ics}{Amount of variance explained by the each independent 
#'      component}
#' @examples
#' data("samples_data")
#' cica <- consICA(samples_data$X, ncomp=15, ntry=10, show.every=0)
#' var_ic <- estimateVarianceExplained(cica)
#' @export
estimateVarianceExplained <- function(cica, X=NULL) {
    if(!is.consICA(cica))  return (NULL)
    if(!is.null(cica$X))  X<-cica$X
    if(is.null(X))  return (NULL)
    if(!is.matrix(X))  return (NULL)
    
    S <- cica$S # G x C
    M <- cica$M # C x S
    
    if(nrow(X) != nrow(S) | ncol(X) != ncol(M) | ncol(S)!= nrow(M))return(NULL)
    
    Xn <- X - mean(X,na.rm=TRUE) ##### check consICA script
    denom <-  sum( (Xn - mean(Xn,na.rm=TRUE))^2,na.rm=TRUE)
    
    R2 <- 1 - (sum(( Xn - S %*% M )^2, na.rm=TRUE)/denom)
    
    ### could be parallelized 
    R2_ic <- 1 - vapply(FUN.VALUE = 0, #/mc
        1:ncol(S),
        function(ic)
            sum((Xn - matrix(S[,ic],ncol=1) %*% matrix(M[ic,],nrow=1) )^2, 
                na.rm=TRUE)/denom
    )
    names(R2_ic) <- colnames(S)
    return(list(R2=R2,R2_ics = R2_ic))
}

#' @title Barplot variance explained by each IC
#' @description Method to plot variance explained (R-squared) by the MOFA model 
#' for each view and latent factor. \cr
#' As a measure of variance explained for gaussian data we adopt the 
#' coefficient of determination (R2). \cr
#' For details on the computation see the help of the 
#' \code{\link[consICA]{estimateVarianceExplained}} function
#' @param cica consICA compliant list with 
#' @param sort specify the arrangement as 'asc'/'desc'. No sorting if NULL
#' @param las orientation value for the axis labels (0 - always parallel to the 
#' axis, 1 - always horizontal, 2 - always perpendicular to the axis, 
#' 3 - always vertical)
#' @param title character string with title of the plot
#' @param x.cex specify the size of the tick label numbers/text with a numeric 
#' value of length 1
#' @param ... extra arguments to be passed to \code{\link{barplot}}
#' @return A numeric vector compliant to `barplot` output
#' @export
#' @examples
#' data("samples_data")
#' cica <- consICA(samples_data$X, ncomp=15, ntry=10, show.every=0)
#' p <- plotICVarianceExplained(cica, sort = "asc")
plotICVarianceExplained <- function(cica,
                                    sort=NULL,
                                    las=2,
                                    title = "Variance explained per IC",
                                    x.cex = NULL,
                                    ...){
    
    if(!is.consICA(cica)) return (NULL)
    if(is.null(cica$eVarExpl$R2_ics)) 
        cica$eVarExpl = estimateVarianceExplained(cica)
    if(is.null(cica$eVarExpl$R2_ics)) return (NULL)
    if(!is.numeric(x.cex)) x.cex = NULL
    
    r2_ics <- cica$eVarExpl$R2_ics
    
    if(!is.null(sort)){
        if(sort == "asc")
            r2_ics = sort(r2_ics,decreasing = FALSE)
        if(sort == "desc")
            r2_ics = sort(r2_ics,decreasing = TRUE)
    }
    
    if(is.null(x.cex)){
        p <- barplot(r2_ics,main=title,las=las, ...)
    } else{
        p = barplot(r2_ics,main=title,las=las,xaxt="n", ...)
        axis(side = 1,at=p,tick = FALSE,labels = names(r2_ics),
             las=2,cex.axis=x.cex)
    }
    return(p)
}
