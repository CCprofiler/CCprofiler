# CCprofiler
*CCprofiler* is an open-source R-package for the analysis of co-fractionation MS datasets, for example SEC-SWATH-MS. It offers a collection of R functions for quality control and filtering, protein quantificaion, protein- and complex-centric data analysis and vsualization. For detailed description of the functions please see the vignette that comes with the software.


## install
Install the dependencies:

```
cran_deps <- c('proxy', 'doSNOW', 'igraph', 'pracma', 'Rmpfr') 
install.packages(cran_deps)
BiocManager::install("qvalue")

```

