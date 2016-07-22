### HW3: Shrinkage Estimation: Ridge Regression Versus Bayesian

In LAB 2 of the handout we discussed Ridge Regression (2.1 and 2.2) and the Bayesian view of it (2.3). This material shows that
the Ridge Regression estimate can be viewed as the posterior mode of a Bayesian model with a Gaussian likelihood and a Gaussian
prior for effects. The equivalence holds for lambda=varE/varB.

We have seen (see example 2.1 of the handout) that the regularization parameter (labmda) controls the extent of srhinkage of estimates of effects
and with this: (i) the ability of the model to fit the training data, (ii) model degree of freedom, and (iii) prediction accuracy in testing data sets.

**Choosing the regularization parameter using cross-validation**

Example 2.1. shows how to search for a value of lambda using a training-testing partition. Here is a simplified version of the code.

```R
rm(list=ls())
 ##### DATA #############################################
 library(BGLR); data(wheat); X=scale(wheat.X)/sqrt(ncol(wheat.X)); Y=wheat.Y;y<-Y[,1]; N<-nrow(X) ; p<-ncol(X)
 set.seed(10001010)
 tst<-sample(1:N,size=150,replace=FALSE)
 XTRN<-X[-tst,]; yTRN<-y[-tst]; XTST<-X[tst,]; yTST<-y[tst]

 ## FITTING MODEL OVER A GRID OF VALUES OF lambda
 h2=seq(from=.05,to=.95,by=.05)
 lambda<-(1-h2)/h2
 
 sqCorTST<-numeric()
 C0<-crossprod(XTRN)	
 rhs<-crossprod(XTRN,yTRN)

 for(i in 1:length(lambda)){ #loop over values of lambda
   C<-C0
   diag(C)=diag(C0)+lambda[i]
   CInv<-chol2inv(chol(C))
   sol<-crossprod(CInv, rhs)
   yHatTST<-XTST%*%sol
   sqCorTST[i]<- cor(yTST,yHatTST)^2
   print(i)
  }  
  
  plot(sqCorTST~h2,type='o',col=4,ylab='R2-TST',xlab='h2')
  abline(v=h2[which(sqCorTST==max(sqCorTST))],col=2,lty=2)

```

