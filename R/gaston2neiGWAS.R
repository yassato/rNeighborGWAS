#' Convert [gaston]'s bed.matrix data to rNeighborGWAS genotype data.
#'
#' A function convert a bed.matrix dataset to rNeighborGWAS genotype data.
#' @param x A [bed.matrix] created using the \code{gaston} package (Perdry & Dandine-Roulland 2020).
#' @return A list including an individual x marker matrix, a data.frame including chromosome numbers in the first column, and SNP positions in the second column, and a numeric vector including phenotypes for individuals.
#' \itemize{
#'  \item{\code{geno}} {Individual x marker matrix}
#'  \item{\code{gmap}} {Data.frame including chromosome numbers in the first column, and SNP positions in the second column}
#'  \item{\code{pheno}} {Numeric vector including phenotypes for individuals}
#' }
#' @references
#' Perdry H, Dandine-Roulland C. (2020) gaston: Genetic Data Handling (QC, GRM, LD, PCA) & Linear Mixed Models. https://CRAN.R-project.org/package=gaston
#' @author Yasuhiro Sato (\email{sato.yasuhiro.36c@kyoto-u.jp})
#' @import gaston
#' @details 
#' This function converts genotype data into -1, 0, or 1 digit as the rNeighborGWAS format. Zero indicates heterozygotes.
#' @examples
#' data("TTN", package="gaston")
#' x <- gaston::as.bed.matrix(TTN.gen, TTN.fam, TTN.bim)
#' g <- gaston2neiGWAS(x)
#' @export
gaston2neiGWAS <- function(x) {
  gmap <- x@snps[,c("chr","pos")]
  if (all(is.na(gmap))) {
    gmap <- NULL
  }
  pheno <- x@ped$pheno
  if (!is.numeric(pheno) || all(is.na(pheno))) {
    pheno <- NULL
  }
  list(geno = gaston::as.matrix(x) - 1L, gmap = gmap, pheno = pheno)
}
