# Fit location-scale family----
fit_locationScale <- function(lambda, dist) {
    # Use method of moments
    muTW <- -1.2065335745820
    sigmaTW <- sqrt(1.607781034581)

    muS <- mean(log(dist))
    sigmaS <- stats::sd(log(dist))

    sigma1 <- sigmaS/sigmaTW
    mu1 <- muS - sigma1 * muTW

    TW <- (log(lambda) - mu1)/sigma1
    pvalue <- ref_dist(TW)

    return(pvalue)
}

# This may change to some other function in the future
ref_dist <- function(x) RMTstat::ptw(x, beta = 1, lower.tail = FALSE, log.p = FALSE)

# Compute test statistic----
compute_largest_root <- function(first_mat, second_mat) {
    first_mat_c <- scale(first_mat, center = TRUE, scale = FALSE)
    W <- crossprod(scale(second_mat, center = TRUE, scale = FALSE))

    svdRes <- corpcor::fast.svd(first_mat_c)
    Xp <- svdRes$v %*% diag(1/svdRes$d)
    C <- crossprod(Xp, W %*% Xp)
    svdC <- corpcor::fast.svd(C)
    Xpp <- svdC$u
    singWeights <- Xp %*% Xpp
    largestRoot <- max(crossprod(singWeights,  W %*% singWeights))

    return(largestRoot)
}

# Permutations----
permutation <- function(X, Y, nperm) {
    if (nrow(X) < nrow(Y)) {
        first_mat <- X
        second_mat <- Y
    } else {
        first_mat <- Y
        second_mat <- X
    }
    combined_data <- rbind(first_mat, second_mat)
    n1 <- nrow(first_mat)
    n2 <- nrow(second_mat)
    replicate(nperm, {
        indices <- sample(1:(n1 + n2), n1)
        first_mat_perm <- combined_data[indices,]
        second_mat_perm <- combined_data[which(!1:(n1 + n2) %in% indices),]
        compute_largest_root(first_mat_perm, second_mat_perm)
    })
}
