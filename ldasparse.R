library(MASS)
source("sparseCov.r")
source("matGen.r")

#lda model with sparsecov estimate 
ldas = function (x, grouping, prior = proportions, alpha = 0.1, type = 'hard', tol = 1e-04, method='mle')
{
  #input data matrix : x, labels : grouping, prior, alpha for sparsecov estimate 
  ##if(is.null(dim(x))) stop('x is not a matrix')
  x <- as.matrix(x)
  ##if(any(!is.finite(x))) stop('infinite, NA or NAN in x')
  n <- nrow(x)
  p <- ncol(x)
  ##if (n != length(grouping)) stop("nrow(x) and length(grouping) are different")
  g <- as.factor(grouping)#all labels
  lev <- lev1 <- levels(g)#all groups
  counts <- as.vector(table(g))
  
  # if(!missing(prior))
  #   {
  #   if (any(prior < 0) || round(sum(prior), 5) != 1) stop("invalid 'p rior'")
  #   if (length(prior) != nlevels(g)) stop("'prior' is of incorrect length")
  #   prior = prior[counts > 0L]
  # }
  
  # if (any(counts == 0L)) 
  #   {
  #   empty = lev[counts == 0L]
  #   warning(sprintf(ngettext(length(empty), "group %s is empty",
  #                            "groups %s are empty"), paste(empty, collapse = " ")), domain = NA)
  #   lev1 = lev[counts > 0L]
  #   g = factor(g, levels = lev1)
  #   counts = as.vector(table(g))
  #   }

  proportions <- counts/n
  ng <- length(proportions)
  names(prior) <- names(counts) <- lev1
  
  #calculate means for groups
  group.means <- tapply(c(x), list(rep(g, p), col(x)), mean) ###neat code with "tapply"
  colnames(group.means) = colnames(x)
  
  f1 <- sqrt(diag(var(x - group.means[g, ])))
  # if (any(f1 < tol)) {
  #   const <- format((1L:p)[f1 < tol])
  #   stop(sprintf(ngettext(length(const), "variable %s appears to be constant within groups", 
  #                         "variables %s appear to be constant within groups"), 
  #                paste(const, collapse = " ")), domain = NA)
  # }
  scaling = diag(1/f1, ,p)
  #fac = if (method == "moment") 1/(n - ng) else 1/n #mle
  x=(x-group.means[g,])%*% scaling
  #use sparsecov function to estimate the cov
  varsp=sparseCov(dat=x, alf=alpha, THRSH = type)
  
  X.s = svd(varsp, nu = 0L)
  rank <- sum(X.s$d > tol)
  if (rank == 0L) stop("rank = 0: variables are numerically constant")
  if (rank < p) warning("variables are collinear")
  scaling <- scaling %*% X.s$v[, 1L:rank] %*% diag(1/sqrt(X.s$d[1L:rank]), , rank)
  
  structure(list(prior = prior, counts = counts, means = group.means, 
                 scaling = scaling, lev = lev, svd = X.s$d[1L:rank], N = n, 
                 call = cl), class = "lda")
}

#predict function for ldas model
predict.ldas = function(object, newdata, prior = object$prior, dimen)
{
  if (is.null(dim(newdata))) stop('no newdata')
  
  x <- as.matrix(newdata)
  ng <- length(object$prior)
  dimen <- if (missing(dimen)) length(object$svd) else min(dimen, length(object$svd)) # dimensionality is ranks as default
  N <- object$N #sample zise
  
  dist <- matrix(0, nrow = nrow(x), ncol = ng)
  dm <- object$means
  
  scaling <- object$scaling
  means <- colSums(prior * object$means) #overall mean
  #project newdata x and groups means dm into lower dimensional space.
  x <- scale(x, center = means, scale = FALSE) %*% scaling
  dm <- scale(object$means, center = means, scale = FALSE) %*% scaling
  
  # for (i in 1L:ng) {
  #   dev <- scale(x, center =  dm[i, ], scale = FALSE)
  #   dist[,i] = -2*log(prior[i])+ diag(dev%*%ginv(varsp)%*%t(dev))
  # }
  
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
  cl <- factor(nm[max.col(posterior)], levels = object$lev)#choose the one with maxposterior
  dimnames(posterior) <- list(rownames(x), nm)
  list(class = cl, posterior = posterior, x = x[, 1L:dimen, drop = FALSE])
}



#predict function given covariance matrix
predict.ldacov = function (newdata, data, grouping, prior = proportions, ccov = diag(diag(cov(x-group.means[g,]))))
{
  if (is.null(dim(newdata))) dim(newdata) <- c(1L, length(newdata))
  x <- as.matrix(data)
  n <- nrow(x)
  p <- ncol(x)
  g <- as.factor(grouping)#all labels
  counts <- as.vector(table(g))
  proportions <- counts/n
  ng <- length(proportions)
  names(prior) <- names(counts) <- levels(g)
  
  group.means <- tapply(c(x), list(rep(g, p), col(x)), mean) ###neat code with "tapply"
  ccov=ccov
  colnames(group.means) = colnames(x)
  dist <- matrix(0, nrow = nrow(newdata), ncol = ng)
  for (i in 1L:ng) 
  {
    dev <- scale(newdata, center =group.means[i, ], scale = FALSE)
    dist[,i] = -2*log(prior[i])+ diag(dev%*%solve(ccov)%*%t(dev))
  }
  
  nm <- levels(g)
  cl <- factor(nm[max.col(-dist)], levels = levels(g))#choose the one with min distance
  list(class = cl)
}





