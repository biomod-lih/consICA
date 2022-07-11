test_that("getGO test", {
  
  require(BiocParallel)
  require("fastICA", include.only = c('fastICA'))
  data("samples_data")
  
  cica2 <- consICA(samples_data, ncomp=2, ntry=1, show.every=0) #exp timesave
  GOs <- getGO(cica2, db = "BP")
  
  expect_true(!is.null(GOs$GOBP))
  expect_true(is.numeric(GOs$GOBP$ic01$pos$FDR))
  expect_equal(length(GOs$GOBP), cica2$ncomp)
  expect_true(sum(GOs$GOBP$ic01$pos$Score < 0) == 0)
})
