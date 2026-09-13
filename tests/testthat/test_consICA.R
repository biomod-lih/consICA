test_that("consICA test", {
  
  require(BiocParallel)
  require("fastICA", include.only = c('fastICA'))
  data("samples_data")
  samples_data <- samples_data[1:400, 1:60]
  
  cica10 <- consICA(samples_data, ncomp=6, ntry=4, 
                 show.every=2, filter.thr=4)
  
  cicaNA <- consICA(c(1), ncomp=6, ntry=4)
  
  cica_redu <- consICA(samples_data, ncomp=6, ntry=4, 
                    show.every=2, filter.thr=4, reduced = TRUE)
  
  expect_equal(length(colnames(cica10$S)), 6)
  expect_equal(length(rownames(cica10$M)), 6)
  expect_true(is.null(cicaNA))
  expect_true(is.null(cica_redu$X))
})

test_that("oneICA test", {
  
  require("fastICA", include.only = c('fastICA'))
  data("samples_data")
  samples_data <- samples_data[1:400, 1:60]
  
  ica10 <- oneICA(samples_data, ncomp=6, filter.thr=4)
  icaNA <- oneICA(c(1), ncomp=6)
  ica_redu <- oneICA(samples_data, ncomp=6, filter.thr=4, reduced = TRUE)

  expect_equal(length(colnames(ica10$S)), 6)
  expect_equal(length(rownames(ica10$M)), 6)
  expect_true(is.null(icaNA))
  expect_true(is.null(ica_redu$X))
})

