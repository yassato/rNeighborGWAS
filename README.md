# rNeighborGWAS    
This is a developer version of the rNeighborGWAS package. Please see *vignette("rNeighborGWAS")* for usage.  
CRAN version is available at https://cran.r-project.org/package=rNeighborGWAS.  

## Installation
Please install the package via GitHub using the devtools library as *devtools::install_github("yassato/rNeighborGWAS", repo="master")*.  

## Dependency
Note that the rNeighborGWAS requires the following R packages.  
- gaston
- parallel

## Release Notes
version 1.2.3 (developer version): recreated using R version 4.0.3; asymmetric neighbor effects are implemented.    
version 1.2.2 (CRAN version): partial PVEs provided by calc_PVEnei(); nei_lm() added.  
version 1.2.1: testthat files fixed.  
version 1.2.0: nei_lmm() and gaston2neiGWAS() added; nei_coval() and neiGWAS() refactored.  
version 1.0.0: Initial version registered in CRAN.  
