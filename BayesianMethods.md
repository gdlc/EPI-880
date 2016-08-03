### Bayesian Whole Genome Regression with Genotypes of Nominally Unrelated Individuals

Data can be downloaded from the following [link](https://www.dropbox.com/s/8mfk0dh2oj3ch8j/Z.RData?dl=0).

**Loading the data and computing marker summaries**

```R
 load('~/Dropbox/GENOMIC_PREDICTION_COURSE/data/examples/data/Z.RData')
 library(BGData)
 tmp=colMeans(Z,na.rm=T)/2
 MAF=ifelse(tmp<.5,tmp,1-tmp)
 minMAF=.01
 Z=Z[,MAF>=minMAF]
```

**Computing LD and comparing it with that of in the mice data set**
```R
load('~/Dropbox/GENOMIC_PREDICTION_COURSE/data/examples/data/Z.RData')

p=ncol(Z)
COR.NU=matrix(nrow=p,ncol=500,NA)

for(i in 1:500){
	print(i)
    for(j in 1:(p-i)){
         COR.NU[j,i]=cor(Z[,j],Z[,j+i])
    }    
}

library(BGLR)
data(mice)
Z=mice.X
p=ncol(Z)
COR.MICE=matrix(nrow=p,ncol=500,NA)

for(i in 1:500){
	print(i)
    for(j in 1:(p-i)){
         COR.MICE[j,i]=cor(Z[,j],Z[,j+i])
    }    
}

R2.MICE=colMeans(COR.MICE^2,na.rm=T)
R2.NU=colMeans(COR.NU^2,na.rm=T)
plot(R2.MICE,col=2,ylim=c(1,0))
points(R2.NU,col=4)


```
**Computing genomic relationships**

```R
 Z=scale(Z)/sqrt(ncol(Z))
 G=tcrossprod(Z)
 plot(diag(G))
 tmp=rowSums(G>.1)-1 
 plot(tmp)
```


**Assessment of population structure**
```R
 EVD=eigen(G)
 plot(c(0,cumsum(EVD$values)),x=0:nrow(G),type='o',col=2,cex=.1)
 plot(EVD$vectors[,1:2])
```

**Toy simulation**
```R
 h2=.25
 nQTL=10
 QTL=floor(seq(from=50,to=ncol(Z),length=nQTL))
 n=nrow(Z); p=ncol(Z)
 b=rep(0,p)
 b[QTL]=runif(nQTL,min=.1,max=.5)*rep(c(-1,1),ceiling(nQTL/2))[1:nQTL]
 signal=Z%*%b
 K=sqrt(h2)/sd(signal)
 b=b*K
 signal=Z%*%b
 error=rnorm(n,sd=sqrt(1-h2))
 y=signal+error
 var(signal)/var(y)
```

**Various Bayesian Models**
```R
 library(BGLR)

 # GBLUP
  fmGBLUP=BGLR(y=y,ETA=list(list(K=G,model='RKHS')),nIter=12000,burnIn=2000,saveAt='GBLUP_')
 
 # BRR
  ETA=list(list(X=Z,model='BRR'))
  fmBRR=BGLR(y=y,ETA=ETA,nIter=12000,burnIn=2000,saveAt='BRR_')

 # BayesA
  ETA=list(list(X=Z,model='BayesA'))
  fmBA=BGLR(y=y,ETA=ETA,nIter=12000,burnIn=2000,saveAt='ByesA_')

 # BayesB
  ETA=list(list(X=Z,model='BayesB'))
  fmBB=BGLR(y=y,ETA=ETA,nIter=12000,burnIn=2000,saveAt='ByesB_')
  
 # BayesC
  ETA=list(list(X=Z,model='BayesC'))
  fmBC=BGLR(y=y,ETA=ETA,nIter=12000,burnIn=2000,saveAt='ByesC_')
  
  
 # Bayesian Lasso
  ETA=list(list(X=Z,model='BL'))
  fmBL=BGLR(y=y,ETA=ETA,nIter=12000,burnIn=2000,saveAt='BL_')
``` 


### Proposed tasks
 - Run the simulation above for nQTL=5,10,20,50,100,500,p  and estimate the proportion of non-zero effects (`probIn`) with BayesB
 - Modify the above code by introducing a testing set, compare the prediction accuracy of `BRR`, `BayesA` and `BayesB` over a number of TRN-TST partitions.
 - Modify the simultation by including dominance, estimate using G-BLUP the heritability for an additive only moodel and an additive+dominance model using G-BLUP.
