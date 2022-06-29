test_that("getFeatures test", {
  
  require(BiocParallel)
  require("fastICA", include.only = c('fastICA'))
  data("samples_data")
  
  cica <- consICA(samples_data, ncomp=15, ntry=4)
  features <- getFeatures(cica)
  
  expect_equal(length(features), cica$ncomp)
  
  for(ic in seq_len(cica$ncomp)){
    expect_true(!is.null(features[[ic]]$pos$features))
    expect_true(!is.null(features[[ic]]$pos$fdr))
    expect_true(features[[ic]]$pos$fdr >= 0)
    expect_true(!is.null(features[[ic]]$neg$features))
    expect_true(!is.null(features[[ic]]$neg$fdr))
    expect_true(features[[ic]]$neg$fdr >= 0)
  }
})
