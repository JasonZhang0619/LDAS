#install.packages("cvTools")
library(cvTools)
library(e1071)
library(HDclassif)
library(rda)
library(penalizedLDA)
source("ldasparse.r")
source("functions.r")

data='MP'
switch(data,
       MP={#loading MP data
         # source("MP/MP.Library.R")
         # gct.file='MP/boston.plus.cell.lines.maxed.2.gct'
         # MPdata=MP.Gct2Frame(gct.file)
         # datax=t(as.matrix(MPdata$ds))
         # cls.file='MP/boston.plus.cell.lines.2.cls'
         # CLS=MP.ReadClsFile(cls.file)
         # datay=as.factor(CLS$class.list)
         load('data/MP.rda')
         },
       
       Khan={#loading Khan data
         # load("khan.rda")
         # datax=t( cbind(as.matrix(khan$train),as.matrix(khan$test)) )
         # datay=c(as.vector(khan$train.classes),as.vector(khan$test.classes))
         # datay=as.factor(datay)
         load('data/Khan_combined.rda')
       }
       )

cat('The',data,'data includes',nrow(datax),'cases in',ncol(datax),'dimensions space covering',length(levels(datay)),'classes\n')
#thresholding types and false positive rate for sparse estimate
types = c('hard', 'soft', 'scad', 'adpt')
alfs=c(0.1,0.05,0.01,0.001)

#methods to compare with
methods=c('LDA','SVM','NB','RDA','PenalizedLDA','HDDA')#,'sda',)
funcs=c(lda,svm,NA,rda,PenalizedLDA,hdda)#,sda,)
names(funcs)=methods

datax1=datax[-1,]
datay1=datay[-1]
samplesize=nrow(datax)

#k=10 folds CV
k <- 10### # of folds
folds <- cvFolds(samplesize, K=k)
pred.sparse=array(0,dim=c(length(alfs),length(types),samplesize),dimnames=list(paste('alfs =',alfs),types))
pred.methods=array(0,dim=c(length(methods),samplesize),dimnames = list(methods))
result=data.frame()
for (p in 1:k){
  result.k=data.frame()
  cat("working on",p,'/',k,"-th fold:\n")
  trainx <- datax[folds$subsets[folds$which != p], ] #Set the training set
  trainy <- datay[folds$subsets[folds$which != p]]
  validationx <- datax[folds$subsets[folds$which == p], ]
  validationy <- datay[folds$subsets[folds$which == p]]
  cat('LDAS\n')
  for (alf in alfs)
    for (type in types){
      cat('alfs=',alf,'types=',type,'\n')
      ypred = predict.ldas( ldas(trainx, trainy, alpha = alf ,type = type ),validationx)$class
      result.k=rbind(result.k,data.frame(ytrue=validationy,ypred=ypred,method='LDAS',alf=alf,type=type))
    }
  cat('other methods\n')
  for(name in methods){
    cat(name,'\n')
    ypred = format.methods(trainx,trainy,validationx,name,funcs[[name]])
    result.k=rbind(result.k,data.frame(ytrue=validationy,ypred=ypred,method=name,alf=NA,type=NA))
  }
  result.k$fold=p
  result=rbind(result,result.k)
}

save(pred.sparse,pred.methods,folds,file=paste(data,'result_dimension.rda',sep=''))

exsit=file.exists('real/MP.csv')
write.table(result,"real/MP.csv",row.names = FALSE, col.names = !exsit, sep = ",", append=TRUE)
