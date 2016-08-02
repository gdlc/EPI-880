
```R
 load('~/Dropbox/GENOMIC_PREDICTION_COURSE/data/examples/data/Z.RData')
 h2=.25
 nQTL=10
 QTL=floor(seq(from=50,to=ncol(Z),length=nQTL))


 n=nrow(Z)
 p=ncol(Z)
 b=rep(0,p)
 b[QTL]=runif(nQTL,min=.1,max=.5)*rep(c(-1,1),ceiling(nQTL/2))[1:nQTL]
 signal=Z%*%b
 
 K=sqrt(h2)/sd(signal)
 b=b*K
 signal=Z%*%b
 
 error=rnorm(n,sd=sqrt(1-h2))
 y=signal+error
 
```

```R
 library(BGData)
 G=getG(Z)
 

``` 
 
