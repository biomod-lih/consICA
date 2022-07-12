#' @title Get features from consICA deconvolution result
#' @description Extract names of features (rows in `X` and `S` matrices) and
#' their false discovery rates values
#' @param cica  list compliant to `consICA()` result
#' @param alpha value in [0,1] interval. Used to filter features with 
#' FDR < `alpha`. Default value is 0.05
#' @param sort sort features decreasing FDR. Default is FALSE
# @usage getFeatures(cica,0.05)
#' @return list of dataframes `pos` for positive and `neg` for negative 
#' affecting features with columns:
#'     \item{features}{names of features}
#'     \item{fdr}{false discovery rate value}
#' @author Petr V. Nazarov
#' @examples
#' data("samples_data")
#' # Get deconvolution of X matrix
#' cica <-  consICA(samples_data, ncomp=10, ntry=1, show.every=0)
#' # Get features names and FDR for each component
#' features <- getFeatures(cica)
#' # Positive affecting features for first components are
#' ic1_pos <- features$ic.1$pos
#' @export
getFeatures <- function(cica, alpha = 0.05, sort = FALSE){
    if(!is.consICA(cica)) return (NULL)
    if(alpha >= 1 | alpha <= 0) return(NULL)
    
    Features <- lapply(colnames(cica$S), function(icomp, S = cica$S) {
      z  <- S[,icomp]
      z  <- (z - median(z)) / mad(z)
      pv <- pnorm(z)
      pv[pv>0.5] <- 1-pv[pv>0.5]
      fdr <- p.adjust(pv,method="BH")
      fdr.pos <- fdr
      fdr.neg <- fdr
      fdr.pos[fdr > alpha | z<0]=1
      fdr.neg[fdr > alpha | z>0]=1
      features <- rownames(cica$S)
      list(
        pos = data.frame(
          features = features[fdr.pos<1],
          fdr = fdr[fdr.pos<1],
          stringsAsFactors=FALSE
        ),
        neg = data.frame(
          features = features[fdr.neg<1],
          fdr = fdr[fdr.neg<1],
          stringsAsFactors=FALSE
        )
      )
    })
    names(Features) <- colnames(cica$S)
    
    if(sort){
      Features <- sortFeatures(Features)
    }
    return(Features)
}
