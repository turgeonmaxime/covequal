#' Test for equality of covariance matrices
#'
#' Uses Roy's union-intersection principle for testing for equality of covariance matrices between
#' two samples. Also provides p-values.
#'
#' @param X matrix of size n1 x p
#' @param Y matrix of size n2 x p
#' @param inference Method for computing p-value.
#' @param nperm Number of permutations. See details.
#' @return A list containing the test statistic and the p-value.
#' @export
#' @importFrom corpcor fast.svd
#' @importFrom RMTstat ptw
test_covequal <- function(X, Y, inference = c("TW", "permutation"), nperm) {
    stopifnot(is.matrix(X), is.matrix(Y),
              ncol(X) == ncol(Y))

    inference <- match.arg(inference, c("TW", "permutation"))
    if (missing(nperm)) {
        nperm <- switch(inference,
                        "TW" = 25,
                        "permutation" = 500)
    }
    # stopifnot(is.integer(nperm))

    # p <- ncol(X)
    n1 <- nrow(X)
    n2 <- nrow(Y)

    # Compute test statistic
    test_stat <- if (n1 < n2) compute_largest_root(X, Y) else compute_largest_root(Y, X)

    # Perform permutations
    null_dist <- permutation(X, Y, nperm)

    # Compute p-value
    p_value <- switch(inference,
                      "TW" = fit_locationScale(test_stat, null_dist),
                      "permutation" = mean(null_dist < test_stat))

    return(list("test_statistic" = test_stat,
                "pvalue" = p_value))

}
