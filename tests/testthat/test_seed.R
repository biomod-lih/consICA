test_that("consICA(seed=) is reproducible on a single core", {

    require("fastICA", include.only = c('fastICA'))
    data("samples_data")
    X <- samples_data[1:400, 1:60]

    a1 <- consICA(X, ncomp = 3, ntry = 2, ncores = 1, seed = 1, show.every = 0)
    a2 <- consICA(X, ncomp = 3, ntry = 2, ncores = 1, seed = 1, show.every = 0)
    expect_equal(a1$S, a2$S)
})

test_that("consICA(seed=) is reproducible in parallel", {

    skip_on_bioc()

    require("fastICA", include.only = c('fastICA'))
    data("samples_data")
    X <- samples_data[1:400, 1:60]

    p1 <- consICA(X, ncomp = 3, ntry = 2, ncores = 2, seed = 1, show.every = 0)
    p2 <- consICA(X, ncomp = 3, ntry = 2, ncores = 2, seed = 1, show.every = 0)
    expect_equal(p1$S, p2$S)
    expect_equal(p1$M, p2$M)

    # different seed
    p3 <- consICA(X, ncomp = 3, ntry = 2, ncores = 2, seed = 99, show.every = 0)
    expect_false(isTRUE(all.equal(p1$S, p3$S)))
})
