

###Modeling Genetic By Environmental Interactions. 

In the course we have focused on models for a single trait/environment. These models can be extended for the analysis of multiple-traits 
as well as for analysis of multi-environment data. This is described in [Lopez-Cruz et al. (2015)](http://www.ncbi.nlm.nih.gov/pubmed/25660166). 

The following entry in the BGLR webpage also provides information about the model [link](https://github.com/gdlc/BGLR-R/blob/master/inst/md/GxE_usingInteractions.md).


**Objectives**: Compare stratified analyses and interaction model based on:

  - Within-environment estimated genomic heritability
  - DIC/ pD (see `fm$fit` in BGLR)
  - Within-environment prediction accuracy

**Data**: Wehat data set in BGLR (environments 1-3)

---------------------------------------------------------------------------------------------------------------------------------------

###Modeling Genetic Heterogneity using Interactions

In the models covered in the course we have assumed that the regression of phenotypes on markers is homogeneous across subjects.

Human, plant and animal genomes usually cluster into groups that reflect the eselction/migration history of the population. This can lead
to heterogeneous effects. Effect heterogeneity can be accounted for using interactions. This is described in [de los Campos et al. (2015)](http://link.springer.com/article/10.1007%2Fs13253-015-0222-5). 

The following entry in the BGLR webpage also provides information about the model [link](https://github.com/gdlc/BGLR-R/blob/master/inst/md/heterogeneity_interactions.md).

**Objectives**:  Compare stratified analyses and interaction model based on:

  - Within-environment estimated genomic heritability
  - DIC/ pD (see `fm$fit\` in BGLR)
  - Within-environment prediction accuracy

**Data**: Wehat data set in BGLR (all traits, conduct analysis one trait at a time).

---------------------------------------------------------------------------------------------------------------------------------------

###Estimating the proportion of variance of phenotypes explained by principal components

Using the [singular-value decomposition](https://en.wikipedia.org/wiki/Singular_value_decomposition) of the scaled/centered genotype matrix we can extract the eigenvectors that span the row-space of the genotype matrix. The eigen-vectors are mutually orthogonal; therefore we can decompose the genomic variance into components explained by eigen-vectors. Further details about this can be found in the following article by [Janss et al. (2012)](http://www.genetics.org/content/192/2/693.short).


**Objectives**:  Estimate the proportion of variance explained by the 598 eigenvectors with positive eigenvalues.
 
**Data**: Wehat data set in BGLR (all traits, conduct analysis one trait at a time).  

**Methods**: to estimate the fraction of variance explained by PCs, use the `saveEffects=TRUE`, see the following entry for further details: [link-1](https://github.com/gdlc/BGLR-R/blob/master/inst/md/example_saveEffects.md), [link-2](https://github.com/gdlc/BGLR-R/blob/master/inst/md/heritability.md).


---------------------------------------------------------------------------------------------------------------------------------------

### Modeling network data using Reproducing Kernel Hilbert Spaces

In LAB 4 we have argued that RKHS methods can be used to perfom regressions on almost any type of input sets, including network data. In this project we we will learn how to perform regressions using network data.

**Objectives**: to estimate the proportion of variance of an outcome that can be explained by a network and to evaluate the predictive power that we can achieve using network regressions.

**Data**: you could obtain data sets from the following [website](https://snap.stanford.edu/data/).

**Methods**:
  - The network can be maped into a Kernel using the so-called [diffusion kernel](https://www.semanticscholar.org/paper/Diffusion-Kernels-on-Graphs-and-Other-Discrete-Kondor-Lafferty/6320770fe216ebbba769b9f0a006669b616a03d0/pdf). Since the kernel involves a bandwidht paramter you may want to compare multiple kernels using different values for the bandwidth parameter.

  - Once the kernel is computed we can fit models using [BGLR](https://github.com/gdlc/bglr-r)
  

