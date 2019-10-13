source("sparseCov.r")
ldas= function (x, grouping, prior = proportions,alpha = 0.1, type = 'hard', tol = 1e-04, method='mle', ...) 
{
    if (is.null(dim(x))) stop("'x' is not a matrix")
    x <- as.matrix(x)
    if (any(!is.finite(x))) stop("infinite, NA or NaN values in 'x'")
    n <- nrow(x)
    p <- ncol(x)
    if (n != length(grouping)) stop("nrow(x) and length(grouping) are different")
    g <- as.factor(grouping)
    lev <- lev1 <- levels(g)
    counts <- as.vector(table(g))
    if (!missing(prior)) {
        if (any(prior < 0) || round(sum(prior), 5) != 1) stop("invalid 'prior'")
        if (length(prior) != nlevels(g)) stop("'prior' is of incorrect length")
        prior <- prior[counts > 0L]
    }
    if (any(counts == 0L)) {
        empty <- lev[counts == 0L]
        warning(sprintf(ngettext(length(empty), "group %s is empty", 
            "groups %s are empty"), paste(empty, collapse = " ")), domain = NA)
        lev1 <- lev[counts > 0L]
        g <- factor(g, levels = lev1)
        counts <- as.vector(table(g))
    }
    proportions <- counts/n
    ng <- length(proportions)
    names(prior) <- names(counts) <- lev1
    method <- match.arg(method)
    group.means <- tapply(c(x), list(rep(g, p), col(x)), mean)
    f1 <- sqrt(diag(var(x - group.means[g, ])))
    if (any(f1 < tol)) {
        const <- format((1L:p)[f1 < tol])
        stop(sprintf(ngettext(length(const), "variable %s appears to be constant within groups", 
            "variables %s appear to be constant within groups"), 
            paste(const, collapse = " ")), domain = NA)
    }
    scaling <- diag(1/f1, , p)
    X <- (x - group.means[g, ]) %*% scaling
    
    #sparse covariance matrix used here.
    varsp=sparseCov(dat=X, alf=alpha, THRSH = type)
    X.s <- svd(varsp, nu = 0L)
    rank <- sum(X.s$d > tol)
    if (rank == 0L) stop("rank = 0: variables are numerically constant")
    if (rank < p) warning("variables are collinear")
    scaling <- scaling %*% X.s$v[, 1L:rank] %*% diag(1/sqrt(X.s$d[1L:rank]), , rank)
    
    #fisher method
    fac <- if (method == "moment") 1/(n - ng) else 1/n
    xbar <- colSums(prior %*% group.means)
    fac <- if (method == "mle") 1/ng else 1/(ng - 1)
    X <- sqrt((n * prior) * fac) * scale(group.means, center = xbar, scale = FALSE) %*% scaling
    X.s <- svd(X, nu = 0L)
    rank <- sum(X.s$d > tol * X.s$d[1L])
    if (rank == 0L) 
        stop("group means are numerically identical")
    scaling <- scaling %*% X.s$v[, 1L:rank]
    if (is.null(dimnames(x))) 
        dimnames(scaling) <- list(NULL, paste("LD", 1L:rank, 
            sep = ""))
    else {
        dimnames(scaling) <- list(colnames(x), paste("LD", 1L:rank, sep = ""))
        dimnames(group.means)[[2L]] <- colnames(x)
    }
    cl <- match.call()
    cl[[1L]] <- as.name("lda")
    structure(list(prior = prior, 
                   counts = counts, 
                   means = group.means, 
                   scaling = scaling, 
                   lev = lev, 
                   svd = X.s$d[1L:rank], 
                   N = n, 
                   call = cl), 
              class = "ldas")
}


predict.ldas = function (object, newdata, prior = object$prior, dimen, ...) 
{
    if (!inherits(object, "ldas")) stop("object not of class \"ldas\"")
    if (missing(newdata)) stop('no newdata input')
    if (is.null(dim(newdata))) 
        dim(newdata) <- c(1L, length(newdata))
    x <- as.matrix(newdata)
    
    if (ncol(x) != ncol(object$means)) stop("wrong number of variables")
    if (length(colnames(x)) > 0L && any(colnames(x) != dimnames(object$means)[[2L]])) 
      warning("variable names in 'newdata' do not match those in 'object'")
    ng <- length(object$prior)
    if (!missing(prior)) {
        if (any(prior < 0) || round(sum(prior), 5) != 1) 
            stop("invalid 'prior'")
        if (length(prior) != ng) 
            stop("'prior' is of incorrect length")
    }
    
    means <- colSums(prior %*% object$means)
    scaling <- object$scaling
    x <- scale(x, center = means, scale = FALSE) %*% scaling
    dm <- scale(object$means, center = means, scale = FALSE) %*% scaling
    dimen <- if (missing(dimen)) length(object$svd) else min(dimen, length(object$svd))
    N <- object$N
    
    dist <- matrix(0, nrow = nrow(x), ncol = ng)
    p <- ncol(object$means)
    X <- x * sqrt(N/(N - ng))
    for (i in 1L:ng) {
        nk <- object$counts[i]
        dev <- scale(X, center = dm[i, ], scale = FALSE)
        dev <- 1 + rowSums(dev^2) * nk/(N * (nk + 1))
        dist[, i] <- prior[i] * (nk/(nk + 1))^(p/2) * dev^(-(N - ng + 1)/2)
    }
    posterior <- dist/drop(dist %*% rep(1, ng))
    nm <- names(object$prior)
    cl <- factor(nm[max.col(posterior)], levels = object$lev)
    dimnames(posterior) <- list(rownames(x), nm)
    list(class = cl, 
         posterior = posterior, 
         x = x[, 1L:dimen, drop = FALSE])
}