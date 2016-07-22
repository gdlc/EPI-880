### HW3: Shrinkage Estimation: Ridge Regression Versus Bayesian
**Due**: Monday July 25th 10:00 am by e-mail to gustvoc@msu.edu


In LAB 2 of the handout we discussed Ridge Regression (2.1 and 2.2) and the Bayesian view of it (2.3). This material shows that
the Ridge Regression estimate can be viewed as the posterior mode of a Bayesian model with a Gaussian likelihood and a Gaussian
prior for effects. The equivalence holds for lambda=varE/varB.

We have seen (see example 2.1 of the handout) that the regularization parameter (labmda) controls the extent of srhinkage of estimates of effects
and with this: (i) the ability of the model to fit the training data, (ii) model degree of freedom, and (iii) prediction accuracy in testing data sets.

**Choosing the regularization parameter using cross-validation**

Example 2.1. shows how to search for a value of lambda using a training-testing partition. Here is a simplified version of the code.

*Box 1: Ridge regression with search of lambda using training-testing data*

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
  
  plot(sqCorTST~h2,type='o',col=4,ylab='R2-TST',xlab='h2',ylim=range(sqCorTST)*c(.9,1.1))
  abline(v=h2[which(sqCorTST==max(sqCorTST))],col=2)
  abline(h=sqCorTST[which(sqCorTST==max(sqCorTST))],col=2)

```

**The Bayesian approach**

Alternatively, we can infer variance parameters from data using a Bayesian model. This is discussed in the handout in Example 2.2. Here is a simplified version of a script that we can use to estimate variance components using BGLR.

*Box 1: Bayesian Ridge Regression*

```R
## Continued from the code provided above ##
 library(BGLR)
 fm=BGLR(y=yTRN,ETA=list(list(X=XTRN,model='BRR')),nIter=6000,burnIn=1000)
 varE=scan('varE.dat')[-c(1:100)]
 varB=scan('ETA_1_varB.dat')[-c(1:100)]
 postMeanLambda=mean(varE/varB)
 postMeanH2=mean(varB/(varB+varE))
 abline(v=postMeanH2,col='green',lty=2)
 yHatTST.Bayes=XTST%*%fm$ETA[[1]]$b
 sqCorTST.Bayes=cor(yTST,yHatTST.Bayes)^2
 abline(h=sqCorTST.Bayes,col='green')
```


The above example shows that the Bayesian approach leads to results (in terms of extent of shrinkage and prediction accuracy) similar to the ones you would obtain if you were to choose lambda to mximize prediction accuracy. However, as we have discussed this is just a single training-testing partition and, as with any other estimate, there is uncertainty about it.

**Assigment**

Repeat the examples of Box 1 and Box 2 for 50 different partitions of the data into training and testing. For this you will need to:  (i) comment the line where the seed is set (this will lead to different partitions), and (ii) run the example 50 times, and for each time record the optimal value of h2 for Box 1, the estimated h2 for Box 2 and the sqaured correlation for Box 1 and 2. Provide a summary of your findings with one or two plots.

Add one sentence commenting your results.




