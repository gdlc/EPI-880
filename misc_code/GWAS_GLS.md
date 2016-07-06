
#### A function to evaluate the -2logLikelihood

```R
neg2LogLik<-function(logVar,V,y,d,n=length(y)){
  y<-y-mean(y)
  Vy<-crossprod(V,y)
  Vy2<-as.vector(Vy)^2
  varE<-exp(logVar[1])
  varU<-exp(logVar[2])
  lambda<-varU/varE
  dStar<-(d*lambda+1)
  sumLogD<-sum(log(dStar))
  neg2LogLik_1<- ( n*log(varE) + sumLogD )
  neg2LogLik_2<- (sum(Vy2/dStar))/varE
  out<- neg2LogLik_1+neg2LogLik_2
  return(out)
}
```

#### A function to get GLS estimates
```R
GWAS.GLS<-function(y,W=NULL,X,G=NULL,V=NULL,d=NULL,vU,vE,verbose=F, getContrasts=function(x){return(as.matrix(x)) }, bruteForce=T){
  ##  GLS eith eigenvectors
  #X'[VDV'D*vG+IvE]^-1X=X'V[D*k+I]^-1V'X/vE
  # Z=SV'X  where S=diag(1/sqrt(d*vG+vE))
  # Z'Z= X'VS'SV'X
  if(is.null(V)){ EVD=eigen(G) }
  V=EVD$vectors
  d=EVD$values
  n=length(y)
  p=ncol(X)
  if(is.null(W)){ 
    W=matrix(nrow=n,ncol=1,1) 
  }else{
    W=cbind(1,W)
  }
  dStar=d*vU+vE
  OUT=matrix(nrow=p,ncol=3,NA)
  colnames(OUT)<-c('RSS-Dif','df','pValue')
  for(i in 1:ncol(V)){ V[,i]=V[,i]/sqrt(dStar[i]) }  # V%*%diag(1/sqrt(d*vG+vE))
  W=crossprod(V,W)
  y=crossprod(V,y)
  C0=crossprod(W)
  rhs0<-crossprod(W,y)
  CInv0=chol2inv(chol(C0))
  sol0<-crossprod(CInv0,rhs0)
  error0<-y-W%*%sol0
  WRSS0<-sum(error0^2)
 # TMP<-getContrasts(x=X[,i])
  TMP<-getContrasts(X=X[,i], model="step")
  dimZ<-ncol(TMP)
  dimC<-ncol(C0)+dimZ
  for(i in 1:p){
    TMP<-getContrasts(X=X[,i],model="step")
    Z=crossprod(V,TMP)
    if(bruteForce){
      Z=cbind(Z,W) 
      C<-crossprod(Z)
      rhs=crossprod(Z,y)
    }else{
      C=matrix(nrow=dimC,ncol=dimC,0)
      C[(dimZ+1):dimC,(dimZ+1):dimC]<-C0
      C[1:dimZ,1:dimZ]=crossprod(Z)
      C[1:dimZ,(dimZ+1):dimC]<-crossprod(Z,W)
      C[(dimZ+1):dimC,1:dimZ]<-C[1:dimZ,(dimZ+1):dimC]
      rhs<-c(crossprod(Z,y),rhs0)
    }
    
    CInv=chol2inv(chol(C))
    sol<-crossprod(CInv,rhs)
    OUT[i,2]<-dimZ
    if(bruteForce){
      error=y-Z%*%sol
    }else{
      error=y-Z%*%sol[1:dimZ]-W%*%sol[(dimZ+1):dimC]
    }
    WRSS1=sum(error^2)      
    OUT[i,1]<-WRSS0-WRSS1 # LRT-statistic
    
    if(verbose){ print(i) }
  }
  OUT[,3]=1-pchisq(df=1,q=OUT[,1])
  return(OUT)
}
```

#### Example comparing p-values from OLS, GLS and PC

```R


```
