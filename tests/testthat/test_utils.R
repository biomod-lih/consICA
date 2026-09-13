test_that("is.consICA validator accepts valid and rejects invalid objects", {
  
  expect_false(is.consICA(2026))
  expect_false(is.consICA(NULL))
  expect_false(is.consICA(list(a = 1)))
  
  bad <- list(ncomp = 3, nsamples = 5, nfeatures = 10,
              S = matrix(0, 10, 3), M = matrix(0, 3, 4))
  expect_false(is.consICA(bad))
  
  good <- list(ncomp = 3, nsamples = 5, nfeatures = 10,
               S = matrix(0, 10, 3), M = matrix(0, 3, 5))
  expect_true(is.consICA(good))
})

test_that("sortDataFrame orders rows by key", {
  
  df <- data.frame(a = c(3, 1, 2), b = c("x", "y", "z"),
                   stringsAsFactors = FALSE)
  
  asc <- sortDataFrame(df, "a")
  expect_equal(asc$a, c(1, 2, 3))
  expect_equal(asc$b, c("y", "z", "x"))
  
  desc <- sortDataFrame(df, "a", decreasing = TRUE)
  expect_equal(desc$a, c(3, 2, 1))
})