#' @title Set orientation for independent components 
#' @description Set orientation for independent components as positive in most 
#' enriched direction. Use first element of `GOs` for direction establishment.
#' @param cica list compliant to `consICA()` result. Must contain GO, 
#' see `getGO()`
#' @param verbose logic TRUE or FALSE. Use TRUE for print process steps. 
#' Default is FALSE
# @usage setOrientation(cica)
#' @return cica object after rotation, with rotated `S`, `M` and 
#' added `compsign` which is vector defined rotation: 
#' `S_rot = S * compsign, M_rot = M * compsign, GO_rot = GO * compsign`
#' @note Implemented inside `getGO()` in version >= 1.1.1.
#' @author Petr V. Nazarov
#' @examples
#' \dontrun{
#' data("samples_data")
#' # Get deconvolution of X matrix
#' #cica <-  consICA(samples_data, ncomp=10, ntry=1, show.every=0)
#' cica <- consICA(samples_data, ncomp=2, ntry=1, show.every=0) # timesaving 
#' example
#' GOs <- getGO(cica, db = "BP")
#' # Get already rotated S matrix and Gene Ontologies
#' cica <- getGO(cica, db = "BP")
#' 
#' # Get Gene Ontologies without rotation (actually, you don't need to do this)
#' # This may used for GO generated with version < 1.1.1. Add GO to cica list.
#' cica <- getGO(cica, db = "BP", rotate = FALSE)
#' # Rotate components
#' cica <- setOrientation(cica, verbose = T)
#' # Which components was rotated
#' which(cica$compsign == -1)
#' }
#' @export

setOrientation <- function(cica, verbose = FALSE){
  
  if(!is.consICA(cica)) {
    message("First parameter should be compliant to `consICA()` result\n")
    return (NULL)
  }
  if(is.null(cica$GO)) {
    message("cica must cotain GO. Use getGO().\n")
  }
  
  Res <- cica
  Res$compsign <- double(ncol(cica$S))
  
  dbs <- names(Res$GO)
  
  for(ic in seq.int(1, ncol(Res$S))){
    pos <- min(Res$GO[[dbs[1]]][[ic]]$pos$FDR) 
    neg <- min(Res$GO[[dbs[1]]][[ic]]$neg$FDR)
    if(neg < pos) {
      if(verbose){
        message("Rotating IC",ic,"\n")
      }
      Res$compsign[ic] <- -1
      Res$S[,ic] <- -Res$S[,ic]
      Res$M[ic,] <- -Res$M[ic,]
      for (db in dbs){
        tmp <- Res$GO[[db]][[ic]]$pos
        Res$GO[[db]][[ic]]$pos <- Res$GO[[db]][[ic]]$neg
        Res$GO[[db]][[ic]]$neg <- tmp
      }
    } else {
      Res$compsign[ic] <- 1
    }
  }
  return(Res)
}
