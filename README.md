# CCprofiler
*CCprofiler* is an open-source R-package for the analysis of co-fractionation MS datasets, for example SEC-SWATH-MS. It offers a collection of R functions for quality control and filtering, protein quantificaion, protein- and complex-centric data analysis and vsualization. For detailed description of the functions please see the vignette that comes with the software.


## Install
Install the dependencies:

```{r}
cran_deps <- c('proxy', 'doSNOW', 'igraph', 'pracma', 'Rmpfr') 
install.packages(cran_deps)
bioc_deps <- c('ensembldb','evobiR','preprocessCore','seqinr','Biostrings', 'qvalue')
BiocManager::install(bioc_deps)

```

## Citing CCprofiler
If you use CCprofiler in your research, please cite our publication in Molecular Systems Biology describing the workflow.
Moritz Heusel*, Isabell Bludau*, George Rosenberger, Robin Hafen, Max Frank, Amir Banaei-Esfahani, Ben Collins, Matthias Gstaiger, Ruedi Aebersold 
Complex-centric proteome profiling by SEC-SWATH-MS 
Mol Syst Biol (2019)15:e8438 
DOI: https://doi.org/10.15252/msb.20188438 
