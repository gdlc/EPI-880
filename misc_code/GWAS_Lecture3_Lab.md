Download data from [here](https://www.dropbox.com/s/a40mje6bdkuxwnm/DATA.rda?dl=0)

> LAB 1 : SINGLE MARKER REGRESSIONS WITH REAL HUMAN DATA

```R
rm(list=ls())

library(BGLR)
library(BGData)
library(qqman)

load("~/Downloads/DATA.rda")

### X : genotype matrix for first 3000 SNPs from a genome of size 500,000 SNPs
### X2 : genotype matrix for 2764 equally spaced SNPs from a genome of size 500000 SNPs
### pheno : matrix of phenotypes (height adjusted for age, gender, and race)

dev.off()
y=pheno$AdjustedHeight
TMP=GWAS(y~1,data=new('BGData',geno=X,pheno=data.frame(y=y)),method='lm')
TMP2=GWAS(y~1,data=new('BGData',geno=X2,pheno=data.frame(y=y)),method='lm')

head(TMP)

################  Manhattan plots ##################

plot(-log10(TMP[,4]),cex=.5)
abline(h=1.3,col="red")

plot(-log10(TMP2[,4]),cex=.5)
abline(h=1.3,col="red")

```

> LAB 2 : MULTIPLE COMPARISONS CORRECTION USING SIMULATED PHENOTYPES

```R
dev.off()
rm(list=ls())
load("~/Downloads/DATA.rda")

set.seed(50)
grp=kmeans(X2,centers=2,nstart=30)$cluster
col = matrix(,nrow(pheno),1)
col[which(grp==1)]="red"
col[which(grp==2)]="blue"
plot(pheno$PC1,pheno$PC2,col=col)

### Create X matrix for group 2
X_1 = X[which(grp==2),]

if(length(which(colSums(X_1)==0))>0) X_1 = X_1[,-which(colSums(X_1)==0)]
X.new=scale(X_1); ### Center and scale X matrix
nQTL=20 ### Number of QTL
h2.full = c(0,0.5) ### h2 at which phenotype is simulated
pval_unadjusted=pval_fdr=y=vector('list',length(h2.full))


for(i in 1:length(h2.full)) pval_unadjusted[[i]]=pval_fdr[[i]]=matrix(,ncol(X.new),1)
for(i in 1:length(h2.full)){
  h2=h2.full[i]
  QTL=floor(seq(from=50,to=ncol(X.new),length=nQTL))
  nQTL=length(QTL); n=nrow(X.new)
  b=rnorm(n=nQTL,0,sqrt(h2))
  u<-as.vector(X.new[,QTL]%*%b)
  varE = 1-h2
  error<-rnorm(n=n,sd=sqrt(varE))
  y[[i]]<-u+error
  y_1 = y[[i]]
  TMP<-GWAS(y_1~1,data=new('BGData',geno=X_1,pheno=data.frame(y_1)),method='lm')
  
  pval_unadjusted[[i]]=TMP[,4]
  pval_fdr[[i]]=p.adjust(TMP[,4],"fdr")
  message("h2=",h2)
}

#### qq plots

par(mfrow=c(1,2))
for(i in 1:length(h2.full)){
  qq(pval_unadjusted[[i]])
}

par(mfrow=c(1,2))
for(i in 1:length(h2.full)){
  expected_p=sort(seq(1/ncol(X.new),1,by=(1/ncol(X.new))),decreasing = FALSE)
  observed_p=sort(pval_unadjusted[[i]],decreasing = FALSE)
  plot(-log10(expected_p),-log10(observed_p),cex=.7)
}

#####  Manhattan plots

dev.off()
par(mfrow=c(2,2))
for(i in 1:2){
  plot(-log10(pval_unadjusted[[i]]),cex=.5,col="gray",ylab="-log10P",main=paste0("h2=",h2.full[i]," pval=unadjusted"))
  points(x=QTL,col=2,cex=.7,pch=19,y=-log10(pval_unadjusted[[i]])[QTL])
  abline(h=1.3,col="red")
  abline(h=4.8,col="blue")
}
for(i in 1:2){
  plot(-log10(pval_fdr[[i]]),cex=.5,col="gray",ylab="-log10P",main=paste0("h2=",h2.full[i]," pval=fdr"))
  points(x=QTL,col=2,cex=.7,pch=19,y=-log10(pval_fdr[[i]])[QTL])
  abline(h=1.3,col="red")
}

### Permutation test

n.permute=500
pval_permute = matrix(,n.permute,1)
pval_permute2 = matrix(,ncol(X_1),1)

for(j in 1:n){
  y_1 = sample(y[[2]])
  TMP<-GWAS(y_1~1,data=new('BGData',geno=X_1,pheno=data.frame(y_1)),method='lm')
  pval_permute[j,]=min(TMP[,4])
  names(pval_permute) = names(which(TMP[,4]==min(TMP[,4])))
  message("h2=",h2," rep=", j)
}

for(k in 1:ncol(X_1)){
  pval_permute2[k,] = (length(which(pval_permute < pval_unadjusted[[2]][k]))+1)/(n.permute+1)
}

#### Plots for all four methods

par(mfrow=c(2,2))
plot(-log10(pval_unadjusted[[2]]),cex=.5,col="gray",ylab="-log10P",main="unadjusted")
points(x=QTL,col=2,cex=.7,pch=19,y=-log10(pval_unadjusted[[2]])[QTL])
abline(h=1.3,col="red")

plot(-log10(pval_unadjusted[[2]]),cex=.5,col="gray",ylab="-log10P",main="unadjusted/bonferroni")
points(x=QTL,col=2,cex=.7,pch=19,y=-log10(pval_unadjusted[[2]])[QTL])
abline(h=4.8,col="blue")

plot(-log10(pval_fdr[[2]]),cex=.5,col="gray",ylab="-log10P",main="FDR")
points(x=QTL,col=2,cex=.7,pch=19,y=-log10(pval_fdr[[2]])[QTL])
abline(h=1.3,col="red")

plot(-log10(pval_permute2),cex=.5,col="gray",ylab="-log10P",main="permutation")
points(x=QTL,col=2,cex=.7,pch=19,y=-log10(pval_permute2[QTL]))
abline(h=1.3,col="red")

```
> LAB 3: POPULATION STRUCTURE CORRECTION USING SIMULATED PHENOTYPES

```R
dev.off()
rm(list=ls())
load("~/Downloads/DATA.rda")

set.seed(50)
PC = cbind(pheno$PC1,pheno$PC2)
X_1 = X
X.new=scale(X_1); 
nQTL=20
h2 = 0.5
reps=10

pval_fdr1=pval_fdr2=matrix(,ncol(X.new),reps)
for(r in 1:reps){
  for(i in 1:length(h2)){
    set.seed(10*r)
    QTL=floor(seq(from=50,to=ncol(X.new),length=nQTL))
    nQTL=length(QTL); n=nrow(X.new)
    b=rnorm(n=nQTL,0,sqrt(h2))
    u<-as.vector(X.new[,QTL]%*%b)
    varE = 1-h2
    error<-rnorm(n=n,sd=sqrt(varE))
    y<-u+error
    
    TMP1=GWAS(y~PC,data=new('BGData',geno=X_1,pheno=data.frame(y=y),map=data.frame()),method='lm')
    TMP2=GWAS(y~1,data=new('BGData',geno=X_1,pheno=data.frame(y=y),map=data.frame()),method='lm')
    
    pval_fdr1[,r]=p.adjust(TMP1[,4],"fdr")
    pval_fdr2[,r]=p.adjust(TMP2[,4],"fdr")
    message("h2=",h2," rep=",r)
  }
}

#### qq plots

par(mfrow=c(1,2))
qq(pval_unadjusted1[-QTL,5])
qq(pval_unadjusted2[-QTL,5])

#####  Manhattan plots

par(mfrow=c(2,1))
plot(-log10(pval_fdr1[,5]),cex=.5,col="gray",ylab="-log10P",main="PC corrected")
points(x=QTL,col=2,cex=.7,pch=19,y=-log10(pval_fdr1[,5])[QTL])
abline(h=1.3,col="red")
plot(-log10(pval_fdr2[,5]),cex=.5,col="gray",ylab="-log10P",main="Unadjusted")
points(x=QTL,col=2,cex=.7,pch=19,y=-log10(pval_fdr2[,5])[QTL])
abline(h=1.3,col="red")


### Type I errors (% false positives)

noncausalSNPs = ncol(X)-nQTL

type1error1=type1error2=matrix(,reps,length(h2.full))
for(i in 1:length(h2.full)){
  for(j in 1:reps){
    type1error1[j,i]=length(which(-log10(pval_fdr1[-QTL,j])>-log10(0.05)))
    type1error2[j,i]=length(which(-log10(pval_fdr2[-QTL,j])>-log10(0.05)))
  }
}

mean(type1error1/noncausalSNPs)
mean(type1error2/noncausalSNPs)

#### Power

power1=power2=matrix(,reps,length(h2.full))
for(j in 1:reps){
    power1[j,i]=length(which(-log10(pval_fdr1[QTL,j])>-log10(0.05)))
    power2[j,i]=length(which(-log10(pval_fdr2[QTL,j])>-log10(0.05)))
}

mean(power1/nQTL)
mean(power2/nQTL)

```
