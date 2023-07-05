#' @title Outside part of multiple run Independent Component Analysis 
#' @description Calculate a common part for consensus Independent Component 
#' Analysis (ICA)
#' @param X matrix with features in rows and samples in columns
#' @param n.comp number of components
#' @param row.norm rows normalization flag. Default value is FALSE
#' @param verbose logic TRUE or FALSE. Use TRUE for print process steps. 
#'     Default value is FALSE
#' @return a list with 
#'         \item{X}{input matrix}
#'         \item{X1}{interim calculated matrix}
#'         \item{K}{pre-whitening matrix that projects data onto the first 
#'        `n.comp` principal components}
#'        
#' @author Maryna Chepeleva
# @examples
# Z <- matrix(rnorm(2000),nrow=100)
# preICA <- outICA(Z, n.comp=5)

outICA <- function(X, n.comp, row.norm = FALSE, verbose = FALSE){
  
  dd <- dim(X)
  d <- dd[dd != 1L]
  if (length(d) != 2L) 
    stop("data must be matrix-conformal")
  X <- if (length(d) != length(dd)) {
    matrix(X, d[1L], d[2L])
  }  else as.matrix(X)
  
  n <- nrow(X)
  p <- ncol(X)
  if (n.comp > min(n, p)) {
    message("'n.comp' is too large: reset to ", min(n, p))
    n.comp <- min(n, p)
  }
  
  if (verbose) 
    message("Centering")
  X <- scale(X, scale = FALSE)
  X <- if (row.norm){
    t(scale(X, scale = row.norm))
  } else t(X)
  
  if (verbose) 
    message("Whitening")
  
  #V <- X %*% t(X)/n
  V <- mat.mult(X, t(X))/n
  s <- La.svd(V)
  
  D <- diag(c(1/sqrt(s$d)))
  #K <- D %*% t(s$u)
  K <- mat.mult(D, t(s$u))
  
  K <- matrix(K[1:n.comp, ], n.comp, p)
  #X1 <- K %*% X
  X1 <- mat.mult(K, X)
  
  return(list(X=X, X1=X1, K=K))
}


#' @title  Fast Independent Component Analysis for multi-run mode
#' @description Adaptation of \code{\link[fastICA]{fastICA}} for quick 
#' multiple-run calculations for consensus Independent Component Analysis (ICA)
#' @param X matrix with features in rows and samples in columns
#' @param n.comp number of components.
#' @param preICA output of `outICA()`. Default is NULL
#' @param alg.typ parameter for fastICA(). If alg.typ == "deflation" the 
#' components are extracted one at a time. If alg.typ == "parallel" the 
#' components are extracted simultaneously. Default value is "deflation"
#' @param fun the functional form of the G function used in the approximation 
#' to neg-entropy in fastICA. Default value is "logcosh"
#' @param w.init initial weights
#' @param alpha default is 1
#' @param row.norm set TRUE if the normalization by rows is needed. 
#' Default is FALSE
#' @param maxit default is 200
#' @param tol default is 1e-04
#' @param verbose logic TRUE or FALSE. Use TRUE for print process steps. 
#'     Default value is FALSE
#' @param verbose logic TRUE or FALSE. Use TRUE for print process steps. 
#'     Default value is FALSE
#' @return a list with (compliant to `fastICA()`output)
#'         \item{X}{pre-processed data matrix}
#'         \item{K}{pre-whitening matrix that projects data onto the first 
#'        `n.comp` principal components}
#'         \item{W}{estimated un-mixing matrix}
#'         \item{A}{estimated mixing matrix}
#'         \item{S}{estimated source matrix}
#'        
#' @author Maryna Chepeleva
# @examples
# Z <- matrix(rnorm(2000),nrow=100)
# ic1 <- coreICA(Z, n.comp=5, alg.typ="deflation",fun="logcosh")
# preICA <- outICA(Z, n.comp=5)
# ic2 <- coreICA(Z, n.comp=5, preICA=preICA, alg.typ="deflation",fun="logcosh")

coreICA <- function (X, n.comp, preICA = NULL, 
                     alg.typ = c("parallel", "deflation"), 
                     fun = c("logcosh", "exp"),  w.init = NULL,
                     alpha = 1, row.norm = FALSE, maxit = 200, 
                     tol = 1e-04, verbose = FALSE) {
  
  if (alpha < 1 || alpha > 2) 
    stop("alpha must be in range [1,2]")
  
  if(is.null(preICA)){
    if (verbose) message("Out ICA part")
    preICA <- outICA(X = X, n.comp = n.comp, 
                     row.norm = row.norm , verbose = verbose)
    if (verbose) message("End Out ICA part")
  } 
  X <- preICA$X
  X1 <- preICA$X1
  K <- preICA$K

  if (is.null(w.init)) {
    w.init <- matrix(rnorm(n.comp^2), n.comp, n.comp)
  } else {
    if (!is.matrix(w.init) || length(w.init) != (n.comp^2)) 
      stop("w.init is not a matrix or is the wrong size")
  }
  
  a <- if (alg.typ == "deflation") 
    ica.R.def(X1, n.comp, tol = tol, fun = fun, alpha = alpha, 
              maxit = maxit, verbose = verbose, w.init = w.init)
  else if (alg.typ == "parallel") 
    ica.R.par(X1, n.comp, tol = tol, fun = fun, alpha = alpha, 
              maxit = maxit, verbose = verbose, w.init = w.init)
  w <- a %*% K

  S <- mat.mult(w, X)
  A <- t(w) %*% solve(w %*% t(w))
  #A <- mat.mult(t(w), solve(mat.mult(w, t(w))))
  return(list(X = t(X), K = t(K), W = t(a), A = t(A), S = t(S)))
}
