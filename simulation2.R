setwd("C:/Users/jason/Dropbox/Thesis/sparse")
source("utils.r")

#setting for generate data
setting <-list(p=100, #dimensions
               T=3, #the distance of group means ranges from 0 to T
               ngrids=20, #n even grids in 0-T
               train=20, #for each grid, take the average of 25 traningsets
               test=20, #for each trainstet, take the average of 20 test obervations
               N=10, #sample size for each group
               balance=TRUE,#training set includes N-1:N-1 or N-1: N
               ng=2) #number of group

covtypes=c("TriDiag", "ArMat", "Banded", "Random")
covtype=1 #which type of covariance we use

types = c('hard', 'soft', 'scad', 'adpt')
alfs=c(0.05,0.01,0.0075,0.005,0.0025,0.001)

methods=c('LDA','SVM','NB','HDDA','RDA','PenalizedLDA')#,'SDA')
funcs=c(lda,svm,naiveBayes,hdda,rda,PenalizedLDA)#,sda)

cat('dimension:',setting$p,'\n',
    'difference of two centroids in one dimension: 0 -', setting$T,' including',setting$n,'grids','\n',
    'covariance type:',covtypes[covtype],'\n',
    'samplesize of each group:',setting$N,'\n')
k = 0.5                            #(must less than 0.5?)
cov = switch(covtype, genTriDiag(setting$p,k), genArMat(setting$p,k),genBandedMat(setting$p,5),genRandSpMat(setting$p))

result=data.frame()
for(t in (0:setting$ngirds)/setting$ngrids*setting$T)
{
  #ngrids different distance of means
  means1 = rep(0, setting$p)
  means2 = t * rep(1, setting$p)
  means = rbind(means1, means2)
  result.t=data.frame()
  for (ntrain in 1:setting$train)
  {
    traindata=data_generate(setting$N*setting$ng,means,cov)
    testdata=data_generate(setting$test/setting$ng,means,cov)
    ytrue=testdata$grouping
    #LDAS
    for (alf in alfs)
      for (type in types)
      {
        fit=ldas(traindata$data,traindata$grouping,alpha = alf, type=type)
        ypred=predict.ldas(fit,testdata$data)$class
        result.t=rbind(result,data.frame(ytrue=ytrue,ypred=ypred,method='LDAS',alf=alf,type=type))
        }
    #other methods
    for(i in 1:length(methods))
    {
      ypred=format.methods(traindata$data,traindata$grouping,testdata$data,methods[i],funcs[i])
      result.t=rbind(result,data.frame(ytrue=ytrue,ypred=ypred,method=methods[i],alf=NA,type=NA))
    }
  }
  result.t$t=t
  result=rbind(result,result.t)
}
result$prec=result$ytrue==result$ypred
