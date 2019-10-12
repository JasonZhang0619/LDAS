library(MASS)
source("sparseCov.r")
source("matGen.r")

#lda model with sparsecov estimate 
ldas = function (x, grouping, prior = proportions, alpha = 0.1, type = 'hard', tol = 1e-04)
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
  
  # f1 <- sqrt(diag(var(x - group.means[g, ])))
  # if (any(f1 < tol)) {
  #   const <- format((1L:p)[f1 < tol])
  #   stop(sprintf(ngettext(length(const), "variable %s appears to be constant within groups", 
  #                         "variables %s appear to be constant within groups"), 
  #                paste(const, collapse = " ")), domain = NA)
  # }
  scaling = diag(1/f1, , p)
  
  x=x-group.means[g,]
  #use sparsecov function to estimate the cov
  varsp=sparseCov(dat=x, alf=alpha, THRSH = type)
  
  list(prior = prior, 
       means = group.means, 
       lev = lev, 
       varsp = varsp)
}

#predict function for ldas model
predict.ldas = function(object, newdata, prior = object$prior)
{
  if (is.null(dim(newdata))) dim(newdata) <- c(1L, length(newdata))
  x <- as.matrix(newdata)
  ng <- length(object$prior)
  
  dist <- matrix(0, nrow = nrow(x), ncol = ng)
  varsp = object$varsp
  dm <- object$means
  
  for (i in 1L:ng) {
    dev <- scale(x, center =  dm[i, ], scale = FALSE)
    dist[,i] = -2*log(prior[i])+ diag(dev%*%ginv(varsp)%*%t(dev))
  }

  nm <- names(object$prior)
  cl <- factor(nm[max.col(-dist)], levels = object$lev)#choose the one with min distance
  list(class = cl)
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





