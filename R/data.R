#' Samples of gene expression
#'
#' A dataset containing the expression of 2454 genes for 472 samples from skin 
#' cutaneous melanoma (SKCM) TCGA cohort, their metadata such as age, gender, 
#' cancer type etc. and survival time-to-event data
#'
#' @format A SummarizedExperiment object:
#' \describe{
#'   \item{assay}{expression matrix with genes by rows and samples by columns}
#'   \item{colData}{data frame with sample metadata (clinical variables)}
#' }
#' @usage data(samples_data)
"samples_data"