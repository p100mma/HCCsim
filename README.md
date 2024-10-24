## HCCsim - Hierarchical Clustering based Correlation Structure Simulator

For in depth description of the method, see:

[1] Piotr Stomma, Witold R Rudnicki, HCS—hierarchical algorithm for simulation of omics datasets, 
Bioinformatics, Volume 40, Issue Supplement_2, September 2024, Pages ii98–ii104, 
https://doi.org/10.1093/bioinformatics/btae392

### Current status

Development version of the package. 
- all important component functions of the package are present with tested examples and documentation. Persistent user could use them to replicate HCR and HCS [1].

To do:
- friendlier interface functions to the main components:
    - main HCR process function, possibly HCD class and methods for visualizing decomposition and assesing its completness, with preset default values. Current versions `initial_clusterNreconstruct()` and `subclusterNreconstruct` are very general.
    - main HCS process function acting on output of HCR.
    - tools for visualization and fit assesment of HCS.
    - automation of finding good fit of metalog distributions to generating PCs, possibly by AIC or BIC criteria.
- prepare vignettes on how to use it with several datasets and how to tune parameters.
- prepare vignettes on how to use MCL [2] with this package to produce suitable initial clustering.

### Known issues

**Big issue to resolve**:
- In [1], we used MCL (Markov clustering algortihm [2]) to perform initial clustering, in particular the implementation of original author (https://micans.org/mcl/).
- Port of this to R for example by Rcpp is tricky because of platform dependent elements present in the implementation (which inclusion could possibly block package from CRAN)
- Existing R implementations are unsuitable for real life size datasets.
- Right now, best course of action is to use this method outside of R and to import the clusters.

**Fit of distributions of individual variables**
- main focus of the simulator is to replicate correlations with reduced number of parameters
- adding random noise from normal distribution works excellent for the purpose of covariance replication
- but in some cases leads to disturbance of shapes of marginal distributions of simulated variables
- this can possibly be resolved by more careful choice of distribution for added noise, for example to adjust for bimodality by using nosie from mixtures of gaussians, etc.

  ### Other references
  [2] Van Dongen, Stijn, Graph clustering via a discrete uncoupling process, Siam Journal on Matrix Analysis and Applications 30-1, p121-141, 2008
