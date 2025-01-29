# consICA: Consensus ICA R-package for multiomics data analysis
consICA implements a data-driven deconvolution method – consensus independent component analysis (ICA) to decompose heterogeneous omics data and extract features suitable for patient diagnostics and prognostics. The method separates biologically relevant transcriptional signals from technical effects and provides information about cellular composition and biological processes [1]. The implementation of parallel computing in the package ensures the efficient analysis on the modern multicore systems.

## Installation
Package is available from bioconductor 3.15 (R version >= 4.1.0)
```r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("consICA")
```

Package is also available on github
```r
library(devtools)
install_github("biomod-lih/consICA")
```  
### For Linux users
As for fast data calculation we use the [RcppGSl](https://github.com/eddelbuettel/rcppgsl) package in dependencies, for Linux it requires a separately installed [libgsl0-dev](https://packages.debian.org/sid/libgsl0-dev) library. Please install it **before** the `consICA` installation:
```bash
sudo apt-get update -y
sudo apt-get install -y libgsl0-dev
```

## Quick start
Read vignette
```r
browseVignettes("consICA")
``` 

## Contact
petr.nazarov@lih.lu

## References
1. Nazarov, P.V., Wienecke-Baldacchino, A.K., Zinovyev, A. et al. Deconvolution of transcriptomes and miRNomes by independent component analysis provides insights into biological processes and clinical outcomes of melanoma patients. BMC Med Genomics 12, 132 (2019). https://doi.org/10.1186/s12920-019-0578-4
2. Chepeleva M, Kaoma T, Muller A et al. сonsICA: Multimodal data deconvolution, integration and elucidation of biological processes in cancer research [version 1; not peer reviewed]. F1000Research 2023, 12:1260 ([Poster](https://f1000research.com/posters/12-1260)). https://doi.org/10.7490/f1000research.1119635.1
