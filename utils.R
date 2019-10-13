# install.packages("HDclassif")
# install.packages("rda")
# install.packages("PDSCE")
# install.packages("e1071")
library(MASS)
library(Matrix)
library(rda)#SCRDA

library(e1071)#SVM
library(HDclassif)#HDDA
library(sparseLDA)#SDA
library(penalizedLDA)#PenalizedLDA
source("matGen.r")
source("ldaS.r")

data_generate = function(N,means, cov , sd = 1)
{
  
  p = ncol(means)
  ng = nrow(means)
  groupmeans = means
  
  n = N*ng # total number of observations each group has N observations
  data = matrix(t(groupmeans), n, p, byrow = TRUE) #data without noise
  groups = 1:ng-1 #start from 0
  grouping = rep(groups, N)#labels of data 0 1 0 1 0 1...
  
  noise = matrix(rnorm(n * p , sd), n, p) #noises are uncorrelated by default
  if (!missing(cov))
  {
    cov.root = chol(cov)
    noise = apply(noise, 1, function(x)
      cov.root %*% x)
    noise = t(noise)
  }
  data = data + noise
  list(data=data,grouping=grouping)
}

#unify the form of different methods
#used in both simulation and real data***
format.methods=function(trainx,trainy,testx,name,FUN)
{
  if (!is.matrix(testx)){testx=t(as.matrix(testx))}
  switch(name,
         SVM={
           fit=svm(trainx,trainy)
           result=as.vector(predict(fit,testx))
         },
         NB={
           fit=naiveBayes(trainx,trainy)
           result=as.vector(predict(fit,testx))
         },
         RDA={#use CV on trainset and use the best alpha and delta 
           fit=rda(t(trainx),trainy) #The columns are sample observations and the rows are variables
           cv=rda.cv(fit,t(trainx),trainy)$cv.err
           index = which(cv==min(cv),arr.ind = T)#(alpha,delta)
           result = as.vector(predict(fit,t(trainx),trainy,t(testx), 
                                      alpha=fit$alpha[index[1]],
                                      delta=fit$delta[index[2]])) ##annoying issue: this step will pop "Fold i :" every time it runs
           result=levels(trainy)[result]
         },
         PenalizedLDA={#use CV on trainset and use the best lambda and K
           trainy_numeric=as.numeric(trainy)# pLda only accept numeric labesls
           cv.out <-PenalizedLDA.cv(trainx,trainy_numeric,lambdas=c(1e-4,1e-3,1e-2,.1,1,10))#,nfold = 10) ##annoying issue: this step will pop "Fold i :" every time it runs
           lambda=cv.out$bestlambda
           K=cv.out$bestK
           fit=PenalizedLDA(trainx,trainy_numeric,lambda=lambda,K=cv.out$bestK)
           result=predict(fit,testx)$ypred[,K]
           result=levels(trainy)[result]
         },
         {#HDDA, LDA and all other methods that follow this format
           fit=FUN(trainx,trainy)
           result=as.vector(predict(fit,testx)$class)}
  )
  result
}

