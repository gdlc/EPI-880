> LAB 1 : SINGLE MARKER REGRESSIONS

```R
rm(list=ls())

library(BGLR)
library(BGData)
library(qqman)
### Download data from [here](https://www.dropbox.com/s/a40mje6bdkuxwnm/DATA.rda?dl=0)

#####################################################
###### LAB 1: SINGLE MARKER REGRESSIONS #############
#####################################################


##### Data ############################################

load("~/Dropbox/UAB/Research/BST880_Lecture3/DATA.rda")

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


