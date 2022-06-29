test_that("consICA test", {
  
  data("samples_data")
  
  cica <- consICA(samples_data, ncomp=10, ntry=4, show.every=0)
  
  surv <- survivalAnalysis(cica, 
        surv = SummarizedExperiment::colData(samples_data)[,c("time", "event")]) 

  expect_true(!is.null(surv$cox.model))
  expect_true(!is.null(surv$hazard.score))
  expect_true(is.numeric(surv$hazard.score))
  expect_true(length(surv$hazard.score) <= cica$ncomp)
})
