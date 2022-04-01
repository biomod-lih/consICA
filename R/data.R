#' Samples of gene expression
#'
#' A dataset containing the expression of 2454 genes for 472 samples form skin 
#' cutaneous melanoma (SKCM) TCGA cohort, their metadata such as age, gender, 
#' cancer type etc. and survival time-to-event data
#'
#' @format A list with tree matrix:
#' \describe{
#'   \item{X}{expression matrix with genes by rows and samples by columns}
#'   \item{Var}{data frame with sample metadata (clinical variables)}
#'   \item{Sur}{survival data frame: time and event (0-alive, 1-dead)}
#' }
#' @usage data(samples_data)
"samples_data"