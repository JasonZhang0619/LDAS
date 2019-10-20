#install.packages("cvTools")
library(cvTools)
library(e1071)
library(HDclassif)
library(rda)
library(penalizedLDA)
source("ldaS.r")
source("utils.r")

data='Khan'
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

#feature selection
library(HiDimDA)
SelectedV <- SelectV(datax,datay,maxp=1000)
datax=datax[,SelectedV$vkpt]

types = c('hard', 'soft', 'scad', 'adpt')
alfs=c(0.1,0.05,0.01,0.001)

#methods to compare with
methods=c('LDA','SVM','NB','RDA','PenalizedLDA','HDDA','Dlda', 'Mlda', 'Slda', 'RFlda')#,'sda',)
funcs=c(lda,svm,naiveBayes,rda,PenalizedLDA,hdda,Dlda, Mlda, Slda, RFlda)#,sda,)
names(funcs)=methods

# datax1=datax[-1,]
# datay1=datay[-1]
samplesize=nrow(datax)

#k=10 folds CV
k <- 10### # of folds
folds <- cvFolds(samplesize, K=k)
# pred.sparse=array(0,dim=c(length(alfs),length(types),samplesize),dimnames=list(paste('alfs =',alfs),types))
# pred.methods=array(0,dim=c(length(methods),samplesize),dimnames = list(methods))
result=data.frame()
for (p in 1:k){
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
      result=rbind(result,data.frame(ytrue=validationy,ypred=ypred,method='LDAS',alf=alf,type=type,index=folds$subsets[folds$which == p]))
    }
  cat('other methods\n')
  for(name in methods){
    cat(name,'\n')
    ypred = format.methods(trainx,trainy,validationx,name,funcs[[name]])
    result=rbind(result,data.frame(ytrue=validationy,ypred=ypred,method=name,alf=NA,type=NA,index=folds$subsets[folds$which == p]))
  }
}

exsit=file.exists('real/MP.csv')
write.table(result,"real/MP.csv",row.names = FALSE, col.names = !exsit, sep = ",", append=TRUE)

cv=data.frame(index=folds$subsets,fold=folds$which)
write.table(cv,"real/MP_CV.csv",row.names = FALSE, sep = ",", append=TRUE)

# exsit=file.exists('real/khan_CV.csv')
# if(exsit) cv=read.csv('real/khan_CV.csv')
# folds=data.frame(subsets=cv$index,which=cv$fold)
# methods=c('Dlda', 'Mlda', 'Slda', 'RFlda')
# result=data.frame()
# funcs=c(Dlda, Mlda, Slda, RFlda)
# names(funcs)=methods
# k=max(folds$which)
# for (p in 1:k){
#   cat("working on",p,'/',k,"-th fold:\n")
#   trainx <- datax[folds$subsets[folds$which != p], ] #Set the training set
#   trainy <- datay[folds$subsets[folds$which != p]]
#   validationx <- datax[folds$subsets[folds$which == p], ]
#   validationy <- datay[folds$subsets[folds$which == p]]
#   cat('other methods\n')
#   for(name in methods){
#     cat(name,'\n')
#     ypred = format.methods(trainx,trainy,validationx,name,funcs[[name]])
#     ypred = levels(trainy)[as.numeric(ypred)+1]
#     result=rbind(result,data.frame(ytrue=validationy,ypred=ypred,method=name,alf=NA,type=NA,index=folds$subsets[folds$which == p]))
#   }
# }
# 
# exsit=file.exists('real/khan.csv')
# write.table(result,"real/khan.csv",row.names = FALSE, col.names = !exsit, sep = ",", append=TRUE)




