#' @title rNeigborGWAS: Testing Neighbor Effects in Marker-based Regressions
#'
#' @description
#' This package provides a set of functions to test neighbor effects in genome-wide association studies.
#' The neighbor effects are estimated using the Ising model of ferromagnetism.
#' See Sato et al. (2021) for motivation and modeling.
#'
#' @docType package
#' @name rNeighborGWAS-package
#' @author Yasuhiro Sato (\email{sato.yasuhiro.36c@kyoto-u.jp})
#' @details
#' The flow of neighbor GWAS consists of two steps, (i) variation partitioning and (ii) association mapping.
#' In the first step, we compute proportion of phenotypic variation explained by neighbor effects, and estimate their effective area.
#' In the second step, we test neighbor effects, and map their association score on a genome.
#' In addition to standard GWAS inputs, spatial information of individuals is required to run these analyses.
#' See \code{vignette("rNeighborGWAS")} for how to use this package.
#' @note
#' A developer version of rNeighborGWAS is available at the GitHub reporsitory (\url{https://github.com/yassato/rNeighborGWAS}).
#' If the CRAN version is out of work or you want to use the latest methods, you may install the package via GitHub.
#' @references
#' Sato Y, Yamamoto E, Shimizu KK, Nagano AJ (2021) Neighbor GWAS: incorporating neighbor genotypic identity into genome-wide association studies of field herbivory. Heredity 126(4):597-614. https://doi.org/10.1038/s41437-020-00401-w
#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
