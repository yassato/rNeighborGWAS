# rNeighborGWAS    
This is a developer version of the rNeighborGWAS package. Please see *vignette("rNeighborGWAS")* for usage.  
CRAN version is available at https://cran.r-project.org/package=rNeighborGWAS.  

## Installation
Please install the package via GitHub using the devtools library as *devtools::install_github("yassato/rNeighborGWAS")*.  

## Dependency
Note that the rNeighborGWAS requires the following R packages.  
- gaston
- parallel

## Release Notes
version 1.2.2 (developer version): partial PVEs separated in calc_PVEnei(); nei_lm() separated from neiGWAS().  
version 1.2.1 (CRAN release version): testthat files fixed.  
version 1.2.0: nei_lmm() and gaston2neiGWAS() added; nei_coval() and neiGWAS() refactored.  
version 1.0.0: Initial version registered in CRAN.  
