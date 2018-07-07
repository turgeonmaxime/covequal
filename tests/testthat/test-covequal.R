context("covequal")

n1 <- 100
n2 <- 100
p <- 200

X <- matrix(rnorm(n1*p), ncol = p)
Y <- matrix(rnorm(n2*p), ncol = p)

test_that("no error in computing lambda", {
    output1 <- try(test_covequal(X, Y, inference = "TW",
                                 nperm = 10),
                   silent = TRUE)
    output2 <- try(test_covequal(X, Y, inference = "permutation",
                                 nperm = 10),
                   silent = TRUE)

    expect_false(inherits(output1, "try-error"))
    expect_false(inherits(output2, "try-error"))
})
