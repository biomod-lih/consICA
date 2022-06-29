test_that("estimateVarianceExplained test", {
  
  require(BiocParallel)
  require("fastICA", include.only = c('fastICA'))
  data("samples_data")
  
  cica <- consICA(samples_data, ncomp=15, ntry=10, show.every=0)
  var_ic <- estimateVarianceExplained(cica)
  
  expect_true(var_ic$R2 >= -1 & var_ic$R2 <= 1)
  expect_true(sum(!(var_ic$R2_ics >= -1 & var_ic$R2_ics <= 1)) == 0)
})

test_that("plotICVarianceExplained test", {
  
  require(BiocParallel)
  require("fastICA", include.only = c('fastICA'))
  data("samples_data")
  
  cica <- consICA(samples_data, ncomp=15, ntry=10, show.every=0)
  p <- plotICVarianceExplained(cica, sort = "asc")
  
  expect_true(is.numeric(p))
  expect_true(length(p) == cica$ncomp)
})
