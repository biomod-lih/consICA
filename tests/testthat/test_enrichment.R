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


test_that("get_score test", {
  genes_test <- c("A2ML1","AACSP1","AADACL2","AARD","ABCA12","ABCA13","ABCA6",    
             "ABCA8","ABCA9","ABCB4","ABCB5","ABCC11","ABCC2","ABCD2",   
             "A2ML1","AACSP1","AADACL2","AARD","ABCA12","ABCA13","ABCA6",        
             "ABCA8","ABCA9","ABCB4","ABCB5","ABCC11","ABCC2","ABCD2")
  
  fdr <- c(1, 1.448701e-02, 1, 1, 1, 
           0,  1.052181e-07, 1, 1, 1, 4.335426e-02,
           1, 1, 0, 3.786880e-02, 1, 1, 2.193489e-02,
           0, 0, 1, 1, 4.383766e-02, 1, 1, 1, 1,1)
  names(fdr) <- genes_test
  fc <- NULL
  ntop <- NA
  thr.fdr <- 0.05
  thr.fc <- NA
  score1 <- get_score(genes = genes_test, fc, thr.fc, fdr, thr.fdr, ntop)
  score2 <- get_score(genes = genes_test, fc, thr.fc, fdr, thr.fdr, ntop = 20)
  
  expect_true(is.numeric(score1))
  expect_true(is.numeric(score2))
  expect_true(sum(score1 < 0) == 0)
  expect_true(sum(score2 < 0) == 0)
})
