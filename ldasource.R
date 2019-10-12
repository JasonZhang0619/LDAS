lda.default = function (x, grouping, prior = proportions, tol = 1e-04, method = c("moment", "mle", "mve", "t"), CV = FALSE, nu = 5, 
          ...) 
{
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  g <- as.factor(grouping)
  lev <- lev1 <- levels(g)
  counts <- as.vector(table(g))
  if (!missing(prior)) {
    if (any(prior < 0) || round(sum(prior), 5) != 1) 
      stop("invalid 'prior'")
    if (length(prior) != nlevels(g)) 
      stop("'prior' is of incorrect length")
    prior <- prior[counts > 0L]
  }
  if (any(counts == 0L)) {
    empty <- lev[counts == 0L]
    warning(sprintf(ngettext(length(empty), "group %s is empty", 
                             "groups %s are empty"), paste(empty, collapse = " ")), 
            domain = NA)
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

  scaling <- diag(1/f1, , p)
  if (method == "mve") {
    cov <- n/(n - ng) * cov.rob((x - group.means[g, ]) %*% 
                                  scaling)$cov
    sX <- svd(cov, nu = 0L)
    rank <- sum(sX$d > tol^2)
    if (rank == 0L) 
      stop("rank = 0: variables are numerically constant")
    if (rank < p) 
      warning("variables are collinear")
    scaling <- scaling %*% sX$v[, 1L:rank] %*% diag(sqrt(1/sX$d[1L:rank]), 
                                                    , rank)
  }
  else if (method == "t") {
    if (nu <= 2) 
      stop("'nu' must exceed 2")
    w <- rep(1, n)
    repeat {
      w0 <- w
      X <- x - group.means[g, ]
      sX <- svd(sqrt((1 + p/nu) * w/n) * X, nu = 0L)
      X <- X %*% sX$v %*% diag(1/sX$d, , p)
      w <- 1/(1 + drop(X^2 %*% rep(1, p))/nu)
      print(summary(w))
      group.means <- tapply(w * x, list(rep(g, p), col(x)), 
                            sum)/rep.int(tapply(w, g, sum), p)
      if (all(abs(w - w0) < 0.01)) 
        break
    }
    X <- sqrt(nu/(nu - 2) * (1 + p/nu)/n * w) * (x - group.means[g, 
                                                                 ]) %*% scaling
    X.s <- svd(X, nu = 0L)
    rank <- sum(X.s$d > tol)
    if (rank == 0L) 
      stop("rank = 0: variables are numerically constant")
    if (rank < p) 
      warning("variables are collinear")
    scaling <- scaling %*% X.s$v[, 1L:rank] %*% diag(1/X.s$d[1L:rank], 
                                                     , rank)
  }
  
  
  else {
    fac <- if (method == "moment") 
      1/(n - ng)
    else 1/n #mle
    
    X <- sqrt(fac) * (x - group.means[g, ]) %*% scaling
    X.s <- svd(X, nu = 0L)
    rank <- sum(X.s$d > tol)
    if (rank == 0L) 
      stop("rank = 0: variables are numerically constant")
    if (rank < p) 
      warning("variables are collinear")
    scaling <- scaling %*% X.s$v[, 1L:rank] %*% diag(1/X.s$d[1L:rank], 
                                                     , rank)
  }
  
  xbar <- colSums(prior %*% group.means)
  fac <- if (method == "mle") 1/ng 
  else 1/(ng - 1)
  
  X <- sqrt((n * prior) * fac) * scale(group.means, center = xbar, 
                                       scale = FALSE) %*% scaling
  X.s <- svd(X, nu = 0L)
  rank <- sum(X.s$d > tol * X.s$d[1L])
  
  if (rank == 0L) stop("group means are numerically identical")
  scaling <- scaling %*% X.s$v[, 1L:rank]
  
  if (is.null(dimnames(x))) 
    dimnames(scaling) <- list(NULL, paste("LD", 1L:rank,sep = ""))
  else {
    dimnames(scaling) <- list(colnames(x), paste("LD", 1L:rank, sep = ""))
    dimnames(group.means)[[2L]] <- colnames(x)
  }
  cl <- match.call()
  cl[[1L]] <- as.name("lda")
  structure(list(prior = prior, counts = counts, means = group.means, 
                 scaling = scaling, lev = lev, svd = X.s$d[1L:rank], N = n, 
                 call = cl), class = "lda")
}



predict.lda = function (object, newdata, prior = object$prior, dimen, method = c("plug-in","predictive", "debiased"), ...) 
{
  x <- as.matrix(newdata)

  ng <- length(object$prior)
  means <- colSums(prior * object$means)#overall mean
  
  scaling <- object$scaling
  x <- scale(x, center = means, scale = FALSE) %*% scaling
  dm <- scale(object$means, center = means, scale = FALSE) %*% scaling
  
  method <- match.arg(method)
  dimen <- if (missing(dimen)) 
    length(object$svd)
  else min(dimen, length(object$svd))
  N <- object$N
  
  if (method == "plug-in") {
    dm <- dm[, 1L:dimen, drop = FALSE]
    dist <- matrix(0.5 * rowSums(dm^2) - log(prior), nrow(x), 
                   length(prior), byrow = TRUE) - x[, 1L:dimen, drop = FALSE] %*% 
      t(dm)
    dist <- exp(-(dist - apply(dist, 1L, min, na.rm = TRUE)))
  }
  else if (method == "debiased") {
    dm <- dm[, 1L:dimen, drop = FALSE]
    dist <- matrix(0.5 * rowSums(dm^2), nrow(x), ng, byrow = TRUE) - 
      x[, 1L:dimen, drop = FALSE] %*% t(dm)
    dist <- (N - ng - dimen - 1)/(N - ng) * dist - matrix(log(prior) - 
                                                            dimen/object$counts, nrow(x), ng, byrow = TRUE)
    dist <- exp(-(dist - apply(dist, 1L, min, na.rm = TRUE)))
  }
  else {
    dist <- matrix(0, nrow = nrow(x), ncol = ng)
    p <- ncol(object$means)
    X <- x * sqrt(N/(N - ng))
    for (i in 1L:ng) {
      nk <- object$counts[i]
      dev <- scale(X, center = dm[i, ], scale = FALSE)
      dev <- 1 + rowSums(dev^2) * nk/(N * (nk + 1))
      dist[, i] <- prior[i] * (nk/(nk + 1))^(p/2) * dev^(-(N - 
                                                             ng + 1)/2)
    }
  }
  posterior <- dist/drop(dist %*% rep(1, ng))
  nm <- names(object$prior)
  cl <- factor(nm[max.col(posterior)], levels = object$lev)
  dimnames(posterior) <- list(rownames(x), nm)
  list(class = cl, posterior = posterior, x = x[, 1L:dimen, 
                                                drop = FALSE])
}