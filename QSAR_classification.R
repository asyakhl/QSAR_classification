rm(list = ls())
total.time=proc.time()
library(glmnet)
library(ggplot2)
qsar=read.table(file="/biodeg.csv", sep=";")
qsar$V42=ifelse(qsar$V42=="RB", 1,0)
qsar=qsar[,-c(19)]

X             =    model.matrix(V42~., qsar)[, -1]
y             =    factor(qsar$V42)
n             =    dim(X)[1] # sample size
p             =    dim(X)[2] # number of predictors/features
S             =    100

X             =    scale(X)

# Err is S x 6 matrix
# column 1 of Err = total train error
# column 2 of Err = total test error

# column 3 of Err = false positive train error 
# column 4 of Err = false negative train error

# column 5 of Err = false positive test error 
# column 6 of Err = false negative test error
# column 7 of Err = minimum CV error


Err.rf     =    matrix(0, nrow = S, ncol = 6) 
Err.rad.svm     =    matrix(0, nrow = S, ncol = 6) 
Err.log        =    matrix(0, nrow = S, ncol = 6) 
Err.lasso        =    matrix(0, nrow = S, ncol = 7)
Err.ridge        =    matrix(0, nrow = S, ncol = 7)
Err.svm = matrix(0, nrow = S, ncol = 7) 
time.method=matrix(0, nrow=S, ncol=7)

for (s in 1:S) {
  
  # randomly splitting the data into test and train set
  random_order  =  sample(n)
  n.train       =  floor(n*0.9)
  n.test        =  n-n.train
  trainSet      =  random_order[1:n.train]
  testSet       =  random_order[(1+n.train):n]
  y.test        =  y[testSet]
  X.test        =  X[testSet, ]
  y.train       =  y[trainSet]
  X.train       =  X[trainSet, ]
  y.os.train    =  y.train   # initialize the over-sampled (os) set to train the models
  X.os.train    =  X.train   # initialize the over-sampled (os) set to train the models
  
  # to take into account the imbalance
  # below we over-sample (with replacement) the  data so that the data is balanced
  imbalance     =     TRUE   
  if (imbalance == TRUE) {
    index.yis0      =      which(y.train==0)  # idetify the index of  points with label 0
    index.yis1      =      which(y.train==1) # idetify the index of  points with label 1
    n.train.1       =      length(index.yis1)
    n.train.0       =      length(index.yis0)
    if (n.train.1 > n.train.0) {     # we need more 0s in out training set, so we over sample with replacement
      more.train    =      sample(index.yis0, size=n.train.1-n.train.0, replace=TRUE)
    }         else {    # we need more 1s in out training set, so we over sample with replacement          
      more.train    =      sample(index.yis1, size=n.train.0-n.train.1, replace=TRUE)
    }
    
    y.os.train        =       as.factor(c(y.train, y.train[more.train])-1) 
    X.os.train        =       rbind2(X.train, X.train[more.train,]) 
    
    
  }
  log.os.train.data=data.frame(X.os.train, as.factor(y.os.train))
  log.train.data=data.frame(X.train, as.factor(y.train))
  log.test.data=data.frame(X.test, as.factor(y.test))
  colnames(log.os.train.data)[41]="y"
  colnames(log.train.data)[41]="y"
  colnames(log.test.data)[41]="y"
  log.fit=glm(y~.,data=log.os.train.data, family=binomial)
  y.train.prob=predict(log.fit, log.train.data,type="response")
  y.train.hat=rep(0, length(y.train.prob))
  y.train.hat[y.train.prob>0.5]=1
  
  y.test.prob=predict(log.fit, log.test.data,type="response")
  y.test.hat=rep(0, length(y.test.prob))
  y.test.hat[y.test.prob>0.5]=1
  
  Err.log[s,1]    =     mean(y.train != y.train.hat)
  Err.log[s,2]    =     mean(y.test != y.test.hat)
  
  # column 3 of Err = false positive train error 
  # column 4 of Err = false negative train error
  Err.log[s,3]  =  mean(1 == y.train.hat[y.train==0]) # false positive
  Err.log[s,4]  =  mean(0 == y.train.hat[y.train==1]) # false negative
  
  # column 5 of Err = false positive test error 
  # column 6 of Err = false negative test error
  Err.log[s,5]  =  mean(1 == y.test.hat[y.test==0]) # false positive
  Err.log[s,6]  =  mean(0 == y.test.hat[y.test==1]) # false negative
  
  # optimize lasso logistic regression using cross validation
  m                 =     10
  ptm.lasso.cv  = proc.time()
  lasso.cv          =     cv.glmnet(X.os.train, y.os.train, family = "binomial", alpha = 1,  
                                    standardize = FALSE, nfolds = 10, type.measure="class")
  ptm.lasso.cv=proc.time()-ptm.lasso.cv
  time.method[s,1] = ptm.lasso.cv[3]
  
  ptm.lasso=proc.time()
  lasso.fit         =     glmnet(X.os.train, y.os.train, lambda = lasso.cv$lambda.min, family = "binomial", alpha = 1)
  ptm.lasso= proc.time()-ptm.lasso
  time.method[s,2]=ptm.lasso[3]
  
  y.train.hat       =     predict(lasso.fit, newx = X.train, type = "class")
  y.test.hat        =     predict(lasso.fit, newx = X.test, type = "class")
  Err.lasso[s,1]    =     mean(y.train != y.train.hat)
  Err.lasso[s,2]    =     mean(y.test != y.test.hat)
  
  # column 3 of Err = false positive train error 
  # column 4 of Err = false negative train error
  Err.lasso[s,3]  =  mean(1 == y.train.hat[y.train==0]) # false positive
  Err.lasso[s,4]  =  mean(0 == y.train.hat[y.train==1]) # false negative
  
  # column 5 of Err = false positive test error 
  # column 6 of Err = false negative test error
  Err.lasso[s,5]  =  mean(1 == y.test.hat[y.test==0]) # false positive
  Err.lasso[s,6]  =  mean(0 == y.test.hat[y.test==1]) # false negative
  Err.lasso[s,7]=min(lasso.cv$cvm)
  
  
  # optimize ridge logistic regression using cross validation
  m                 =     10
  ptm.ridge.cv = proc.time()
  ridge.cv          =     cv.glmnet(X.os.train, y.os.train, family = "binomial", alpha = 0,  nfolds = 10, type.measure="class")
  ptm.ridge.cv = proc.time() -ptm.ridge.cv
  time.method[s,3]=ptm.ridge.cv[3]
  ptm.ridge= proc.time()
  ridge.fit         =     glmnet(X.os.train, y.os.train, lambda = ridge.cv$lambda.min, family = "binomial", alpha = 0)
  ptm.ridge =proc.time()-ptm.ridge
  time.method[s,4]=ptm.ridge[3]
  
  y.train.hat       =     predict(ridge.fit, newx = X.train, type = "class")
  y.test.hat        =     predict(ridge.fit, newx = X.test, type = "class")
  Err.ridge[s,1]    =     mean(y.train != y.train.hat)
  Err.ridge[s,2]    =     mean(y.test != y.test.hat)
  
  # column 3 of Err = false positive train error 
  # column 4 of Err = false negative train error
  Err.ridge[s,3]  =  mean(1 == y.train.hat[y.train==0]) # false positive
  Err.ridge[s,4]  =  mean(0 == y.train.hat[y.train==1]) # false negative
  
  # column 5 of Err = false positive test error 
  # column 6 of Err = false negative test error
  Err.ridge[s,5]  =  mean(1 == y.test.hat[y.test==0]) # false positive
  Err.ridge[s,6]  =  mean(0 == y.test.hat[y.test==1]) # false negative
  Err.ridge[s,7]  =  min(ridge.cv$cvm)
  # random forrest
  # alternative way of breaking data into train and test
  os.train.data           =      data.frame(X.os.train, as.factor(y.os.train))
  train.data              =      data.frame(X.train, as.factor(y.train))
  test.data               =      data.frame(X.test, as.factor(y.test))
  names(os.train.data)[41]=      "y"
  names(train.data)[41]   =      "y"
  names(test.data)[41]    =      "y"
  
  library(randomForest)
  ptm.rf=proc.time()
  rf.fit            =     randomForest(y~., data = os.train.data, mtry = sqrt(p), ntree=300)
  ptm.rf=proc.time()-ptm.rf
  time.method[s,5]=ptm.rf[3]
  y.train.hat       =     predict(rf.fit, newdata = train.data)
  y.test.hat        =     predict(rf.fit, newdata = test.data)
  Err.rf[s,1]       =     mean(y.train != y.train.hat)
  Err.rf[s,2]       =     mean(y.test != y.test.hat)
  # column 3 of Err = false positive train error 
  # column 4 of Err = false negative train error
  Err.rf[s,3]  =  mean(1 == y.train.hat[y.train==0]) # false positive
  Err.rf[s,4]  =  mean(0 == y.train.hat[y.train==1]) # false negative
  
  # column 5 of Err = false positive test error 
  # column 6 of Err = false negative test error
  Err.rf[s,5]  =  mean(1 == y.test.hat[y.test==0]) # false positive
  Err.rf[s,6]  =  mean(0 == y.test.hat[y.test==1]) # false negative
  
  
  
  library(e1071)
  svm.os.train.data=data.frame(X.os.train, y=as.factor(y.os.train))
  svm.train.data       =   data.frame(X.train, y = as.factor(y.train))
  svm.test.data       =   data.frame(X.test, y = as.factor(y.test))
  
  ptm.svm.cv=proc.time()
  tune.svm  =   tune(svm, y~., data=svm.os.train.data, kernel = "radial",
                     ranges = list(cost = 10^seq(-2,2,length.out = 5),
                                   gamma = 10^seq(-2,2,length.out = 5)))
  svm.fit = tune.svm$best.model
  ptm.svm.cv=proc.time()-ptm.svm.cv
  time.method[s,6]=ptm.svm.cv[3]
  
  ptm.svm=proc.time()
  y.train.hat=predict(svm.fit,svm.train.data)
  ptm.svm=proc.time()-ptm.svm
  time.method[s,7]=ptm.svm[3]
  y.test.hat=predict(svm.fit,svm.test.data)
  
  
  Err.svm[s,1]       =     mean(y.train != y.train.hat)
  Err.svm[s,2]       =     mean(y.test != y.test.hat)
  
  # column 3 of Err = false positive train error 
  # column 4 of Err = false negative train error
  Err.svm[s,3]  =  mean(1 == y.train.hat[y.train==0]) # false positive
  Err.svm[s,4]  =  mean(0 == y.train.hat[y.train==1]) # false negative
  
  # column 5 of Err = false positive test error 
  # column 6 of Err = false negative test error
  Err.svm[s,5]  =  mean(1 == y.test.hat[y.test==0]) # false positive
  Err.svm[s,6]  =  mean(0 == y.test.hat[y.test==1]) # false negative
  Err.svm[s,7]  = summary(tune.svm)$best.performance
  
  cat(sprintf("s=%1.f| Test:log=%.2f| lasso=%.2f,ridge=%.2f| rf=%.2f| svm=%.2f| ||| Train:log=%.2f| lasso=%.2f,ridge=%.2f| rf=%.2f| svm=%.2f|\n", 
              s, Err.log[s,2], Err.lasso[s,2],Err.ridge[s,2], Err.rf[s,2], Err.svm[s,2], Err.log[s,1], Err.lasso[s,1],Err.ridge[s,1], Err.rf[s,1], Err.svm[s,1]))
  
}


err.train           =     data.frame(c(rep("log",S), rep("logistic lasso", S),  rep("logistic ridge", S),rep("rf", S),  rep("svm", S)) , 
                                     c(Err.log[, 1], Err.lasso[, 1], Err.ridge[, 1], Err.rf[,1], Err.svm[,1]))
err.test           =     data.frame(c(rep("log",S), rep("logistic lasso", S),  rep("logistic ridge", S),rep("rf", S),  rep("svm", S)) , 
                                    c(Err.log[, 2], Err.lasso[, 2], Err.ridge[, 2], Err.rf[,2], Err.svm[,2]))
# train false positive (fp)
err.train.fp        =     data.frame(c(rep("log",S), rep("logistic lasso", S),  rep("logistic ridge", S),rep("rf", S),  rep("svm", S)) , 
                                     c(Err.log[, 3], Err.lasso[, 3], Err.ridge[, 3], Err.rf[,3], Err.svm[,3]))
# train false negative (fn)
err.train.fn        =     data.frame(c(rep("log",S), rep("logistic lasso", S),  rep("logistic ridge", S),rep("rf", S),  rep("svm", S)) , 
                                     c(Err.log[, 4], Err.lasso[, 4], Err.ridge[, 4], Err.rf[,4], Err.svm[,4]))

# test false positive (fp)
err.test.fp        =     data.frame(c(rep("log",S), rep("logistic lasso", S),  rep("logistic ridge", S),rep("rf", S),  rep("svm", S)) , 
                                    c(Err.log[, 5], Err.lasso[, 5], Err.ridge[, 5], Err.rf[,5], Err.svm[,5]))
# train false negative (fn)
err.test.fn        =     data.frame(c(rep("log",S), rep("logistic lasso", S),  rep("logistic ridge", S),rep("rf", S),  rep("svm", S)) , 
                                    c(Err.log[, 6], Err.lasso[, 6], Err.ridge[, 6], Err.rf[,6], Err.svm[,6]))


colnames(err.train)    =     c("method","err")
colnames(err.test)     =     c("method","err")
colnames(err.train.fp) =     c("method","err")
colnames(err.train.fn) =     c("method","err")
colnames(err.test.fp)  =     c("method","err")
colnames(err.test.fn)  =     c("method","err")

p1 = ggplot(err.train)   +     aes(x=method, y = err, fill=method) +   geom_boxplot()  +
  ggtitle("train errors") +
  theme( axis.title.x = element_text(size = 12, face     = "bold", family = "Courier"),
         plot.title          = element_text(size = 12, family   = "Courier"), 
         axis.title.y        = element_text(size = 12, face     = "bold", family = "Courier"), 
         axis.text.x         = element_text(angle= 45, hjust    = 1, size = 11, face = "bold", family = "Courier"), 
         axis.text.y         = element_text(angle= 45, vjust    = 0.7, size = 11, face = "bold", family = "Courier"),
         legend.text = element_text(size=10, face="bold"),
         legend.title = element_text(size=10, face="bold"))+
  
  ylim(0, 0.4)  

p2 = ggplot(err.test)   +     aes(x=method, y = err, fill=method) +   geom_boxplot()   +
  ggtitle("test errors") +
  theme( axis.title.x = element_text(size = 12, face     = "bold", family = "Courier"),
         plot.title          = element_text(size = 12, family   = "Courier"), 
         axis.title.y        = element_text(size = 12, face     = "bold", family = "Courier"), 
         axis.text.x         = element_text(angle= 45, hjust    = 1, size = 11, face = "bold", family = "Courier"), 
         axis.text.y         = element_text(angle= 45, vjust    = 0.7, size = 11, face = "bold", family = "Courier"),
         legend.text = element_text(size=10, face="bold"),
         legend.title = element_text(size=10, face="bold"))+ ylim(0, 0.4)  
p3 = ggplot(err.train.fp)   +     aes(x=method, y = err, fill=method) +   geom_boxplot()  +  
  ggtitle("train fp errors") +
  theme( axis.title.x = element_text(size = 12, face     = "bold", family = "Courier"),
         plot.title          = element_text(size = 12, family   = "Courier"), 
         axis.title.y        = element_text(size = 12, face     = "bold", family = "Courier"), 
         axis.text.x         = element_text(angle= 45, hjust    = 1, size = 10, face = "bold", family = "Courier"), 
         axis.text.y         = element_text(angle= 45, vjust    = 0.7, size = 10, face = "bold", family = "Courier"),
         legend.text = element_text(size=10, face="bold"),
         legend.title = element_text(size=10, face="bold"))+ylim(0, 0.4)  

p4 = ggplot(err.train.fn)   +     aes(x=method, y = err, fill=method) +   geom_boxplot()  +  
  ggtitle("train fn errors") +
  theme( axis.title.x = element_text(size = 12, face     = "bold", family = "Courier"),
         plot.title          = element_text(size = 12, family   = "Courier"), 
         axis.title.y        = element_text(size = 12, face     = "bold", family = "Courier"), 
         axis.text.x         = element_text(angle= 45, hjust    = 1, size = 10, face = "bold", family = "Courier"), 
         axis.text.y         = element_text(angle= 45, vjust    = 0.7, size = 10, face = "bold", family = "Courier"),
         legend.text = element_text(size=10, face="bold"),
         legend.title = element_text(size=10, face="bold"))+  ylim(0, 0.4)  

p5 = ggplot(err.test.fp)   +     aes(x=method, y = err, fill=method) +   geom_boxplot()  +  
  ggtitle("test fp errors") +
  theme( axis.title.x = element_text(size = 12, face     = "bold", family = "Courier"),
         plot.title          = element_text(size = 12, family   = "Courier"), 
         axis.title.y        = element_text(size = 12, face     = "bold", family = "Courier"), 
         axis.text.x         = element_text(angle= 45, hjust    = 1, size = 10, face = "bold", family = "Courier"), 
         axis.text.y         = element_text(angle= 45, vjust    = 0.7, size = 10, face = "bold", family = "Courier"),
         legend.text = element_text(size=10, face="bold"),
         legend.title = element_text(size=10, face="bold"))+  ylim(0, 0.4)  

p6 = ggplot(err.test.fn)   +     aes(x=method, y = err, fill=method) +   geom_boxplot()  + 
  ggtitle("test fn errors") +
  theme( axis.title.x = element_text(size = 12, face     = "bold", family = "Courier"),
         plot.title          = element_text(size = 12, family   = "Courier"), 
         axis.title.y        = element_text(size = 12, face     = "bold", family = "Courier"), 
         axis.text.x         = element_text(angle= 45, hjust    = 1, size = 10, face = "bold", family = "Courier"), 
         axis.text.y         = element_text(angle= 45, vjust    = 0.7, size = 10, face = "bold", family = "Courier"),
         legend.text = element_text(size=10, face="bold"),
         legend.title = element_text(size=10, face="bold"))+  ylim(0, 0.4)  



library("gridExtra")
grid.arrange(p1, p3, p4, p2, p5, p6, ncol=3)


# plotting minimum CV errors for logistic lasso, logistic ridge, and svm
err.cv.min           =     data.frame(c(rep("logistic lasso", S),  rep("logistic ridge", S),  rep("svm", S)) , 
                                      c(Err.lasso[, 7], Err.ridge[, 7], Err.svm[,7]))

colnames(err.cv.min)= c("method", "min.cv.error")
ggplot(err.cv.min)   +     aes(x=method, y = min.cv.error, fill=method) +   geom_boxplot()  +
  ggtitle("minimum CV errors") +
  theme( axis.title.x = element_text(size = 12, face     = "bold", family = "Courier"),
         plot.title          = element_text(size = 12, family   = "Courier"), 
         axis.title.y        = element_text(size = 12, face     = "bold", family = "Courier"), 
         axis.text.x         = element_text(angle= 45, hjust    = 1, size = 11, face = "bold", family = "Courier"), 
         axis.text.y         = element_text(angle= 45, vjust    = 0.7, size = 11, face = "bold", family = "Courier"),
         legend.text = element_text(size=10, face="bold"),
         legend.title = element_text(size=10, face="bold"))+
  
  ylim(0, 0.4)  

library(glmnet)
all.lambdas=c(lasso.cv$lambda,0)
lasso.betas=glmnet(X.os.train, y.os.train, alpha=1, lambda=all.lambdas, family = "binomial")
lasso.cv.with0lammbda          =     cv.glmnet(X.os.train, y.os.train, family = "binomial", alpha = 1,  
                                               standardize = FALSE, lambda=all.lambdas, nfolds = 10, type.measure="class")
beta.l2norms=c(1:length(all.lambdas))
for (i in 1:length(all.lambdas)){
  beta.l2norms[i]=(coef(lasso.betas)[,i]%*%coef(lasso.betas)[,i])/(coef(lasso.betas)[,length(all.lambdas)]%*%coef(lasso.betas)[,length(all.lambdas)])
}

to.graph.lasso=data.frame(method=rep("lasso", length(all.lambdas)), l2norm=beta.l2norms, cv.error=lasso.cv.with0lammbda$cvm)

all.lambdas=c(ridge.cv$lambda,0)
ridge.betas=glmnet(X.os.train, y.os.train, alpha=0, lambda=all.lambdas, family = "binomial")
ridge.cv.with0lammbda          =     cv.glmnet(X.os.train, y.os.train, family = "binomial", alpha = 0,  
                                               standardize = FALSE, lambda=all.lambdas, nfolds = 10, type.measure="class")
beta.l2norms=c(1:length(all.lambdas))
for (i in 1:length(all.lambdas)){
  beta.l2norms[i]=(coef(ridge.betas)[,i]%*%coef(ridge.betas)[,i])/(coef(ridge.betas)[,length(all.lambdas)]%*%coef(ridge.betas)[,length(all.lambdas)])
}

to.graph.ridge=data.frame(method=rep("ridge", length(all.lambdas)), l2norm=beta.l2norms, cv.error=ridge.cv.with0lammbda$cvm)

combined=rbind(to.graph.lasso, to.graph.ridge)


ggplot(combined)   +     aes(x=l2norm, y = cv.error, group=method) +   geom_line(aes(color=method))  +
  geom_point(aes(color=method))+
  xlab(expression(paste("||",hat(beta)[lambda], "||"[2], "/", hat(beta)[lambda==0], "||"[2])))+
  scale_color_manual(values=c("deeppink3","forestgreen"))+
  ggtitle("Cross-Validation Error Rates: \n LASSO and RIDGE") +
  theme( axis.title.x = element_text(size = 15, face     = "bold", family = "Courier"),
         plot.title          = element_text(size = 12, family   = "Courier",face     = "bold"), 
         axis.title.y        = element_text(size = 15, face     = "bold", family = "Courier"), 
         axis.text.x         = element_text(angle= 45, hjust    = 1, size = 12, face = "bold", family = "Courier"), 
         axis.text.y         = element_text(angle= 45, vjust    = 0.7, size = 12, face = "bold", family = "Courier"),
         legend.text = element_text(size=12, face="bold"),
         legend.title = element_text(size=12, face="bold"),
         legend.position="top")


lasso.coef=as.vector(coef(lasso.fit))[-1]
ridge.coef=as.vector(coef(ridge.fit))[-1]

lasso.coef.data=data.frame(predictor=rownames(coef(lasso.fit))[-1], lasso.coefs=lasso.coef)
ridge.coef.data=data.frame(predictor=rownames(coef(ridge.fit))[-1], ridge.coefs=ridge.coef)

ggplot(lasso.coef.data, aes(reorder(predictor, lasso.coefs), lasso.coefs)) +
  theme(axis.title.y = element_blank(),
        axis.title.x =element_text(size=11,face="bold"),
        axis.text.y  = element_text(vjust=0.7, size=9,face="bold"),
        plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text.x =element_text(size=11,face="bold"))+
  geom_bar(stat = "identity", fill = "deeppink3") +
  coord_flip() +
  geom_text(aes(label = round(lasso.coefs,2)),
            nudge_y = ifelse(lasso.coef.data$lasso.coefs > 0, 0.2, -0.2) )+
  ggtitle("Importance of Logistic.Lasso Parameters")


ggplot(ridge.coef.data, aes(reorder(predictor, ridge.coefs), ridge.coefs)) +
  theme(axis.title.y = element_blank(),
        axis.title.x =element_text(size=11,face="bold"),
        axis.text.y  = element_text(vjust=0.7, size=9,face="bold"),
        plot.title = element_text(hjust=0.5, size=13, face="bold"),
        axis.text.x =element_text(size=11,face="bold"))+
  geom_bar(stat = "identity", fill = "forestgreen") +
  coord_flip() +
  geom_text(aes(label = round(ridge.coefs,2)),
            nudge_y = ifelse(ridge.coef.data$ridge.coefs > 0, 0.05, -0.05) )+
  ggtitle("Importance of Logistic.Ridge Parameters")


rf.var.importance=data.frame(Variable=rownames(rf.fit$importance),
                             MeanDecreaseGini=as.vector(rf.fit$importance))

ggplot(rf.var.importance, aes(reorder(Variable, MeanDecreaseGini), MeanDecreaseGini)) +
  theme(axis.title.y = element_blank(),
        axis.title.x =element_text(size=11,face="bold"),
        axis.text.y  = element_text(vjust=0.7, size=9,face="bold"),
        axis.text.x  = element_text(vjust=0.7, size=9,face="bold"),
        plot.title = element_text(hjust=0.5, size=13, face="bold"))+
  geom_bar(stat = "identity", fill = "forestgreen") +
  coord_flip() +
  geom_text(aes(label = round(MeanDecreaseGini,2)),
            nudge_y = 2 )+
  ggtitle("Variable Importance of Random Forest")

ggplot(data = tune.svm$performances, aes(x = as.factor(cost),
                                         y = as.factor(gamma), fill=error)) +geom_tile()+
  scale_fill_distiller(palette = "RdPu")+
  theme(axis.title.y = element_text(size=12,face="bold"),
        axis.title.x =element_text(size=12,face="bold"),
        axis.text.y  = element_text(vjust=0.7, size=11,face="bold"),
        axis.text.x  = element_text(vjust=0.7, size=11,face="bold"),
        plot.title = element_text(hjust=0.5, size=13, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        legend.title = element_text(size=12, face="bold"))+
  ggtitle("Heatmap of SVM CV Error")+ xlab("Cost")+ylab("Gamma")

time.method_n0.9=data.frame(Method=c("Lasso.CV", "Lasso.Fit", "Ridge.CV", "Ridge.Fit", "RF", "SVM.CV", "SVM.Fit"), Time_sec=colMeans(time.method))
time.method_n0.9
total.time=proc.time()-total.time
