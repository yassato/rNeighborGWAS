#' Simulating phenotypes with self and neighbor effects
#'
#' A function to simulate phenotypes caused by self and neighbor effects, with the proportion of phenotypic variation explained (PVE) by fixed and random effects controlled.
#' @param geno An individual x marker matrix. Bialleles (*i.e.*, A or a) must be converted into -1 or 1 digit.
#' @param smap A matrix showing a spatial map for individuals. The first and second column include spatial points along a x-axis and y-axis, respectively.
#' @param scale A numeric scalar indicating the maximum spatial distance between a focal individual and neighbors to define neighbor effects.
#' @param alpha Distance decay coefficient \eqn{\alpha} in a dispersal kernel. Default is set at Inf, meaning no distance decay.
#' @param grouping A integer vector assigning each individual to a group. This argument can be useful when a "smap" contains different experimental replicates. Default setting means that all individuals are belong to a single group.
#' @param kernel An option to select a negative exponential kernel \code{"exp"} or Gaussian kernel \code{"gaussian"}.
#' @param n_causal No. of causal markers in a simulated phenotype
#' @param pveB Proportion of phenotypic variation explained by fixed effects.
#' @param pve Proportion of phenotypic variation explained by fixed and random effects.
#' @param b_ratio A vector composed of two numeric scalars that control the ratio of contributions of self or neighbor effects to a phenotype. The first and second element are for self and neighbor effects, respectively.
#' @return A vector of simulated phenotype values for all individuals
#' @author Yasuhiro Sato (\email{sato.yasuhiro.36c@kyoto-u.jp})
#' @examples
#' set.seed(1)
#' g <- matrix(sample(c(-1,1),100*1000,replace = TRUE),100,1000)
#' gmap <- cbind(c(rep(1,nrow(g)/2),rep(2,nrow(g)/2)),c(1:ncol(g)))
#' x <- runif(nrow(g),1,100)
#' y <- runif(nrow(g),1,100)
#' smap <- cbind(x,y)
#' grouping <- c(rep(1,nrow(g)/2), rep(2,nrow(g)/2))
#' pheno <- nei_simu(geno=g,smap=smap,scale=44,grouping=grouping,n_causal=50,pveB=0.4,pve=0.8)
#'
#' fake_nei <- list()
#' fake_nei[[1]] <- g
#' fake_nei[[2]] <- gmap
#' fake_nei[[3]] <- smap
#' fake_nei[[4]] <- data.frame(pheno,grouping)
#' names(fake_nei) <- c("geno","gmap","smap","pheno")
#' @export
nei_simu = function(geno, smap, scale, alpha=Inf, grouping=grouping, kernel="exp", n_causal, pveB, pve, b_ratio=c(1,1)) {
  g_nei <- nei_coval(geno=geno, smap=smap, scale=scale, alpha=alpha, grouping=grouping, kernel="exp")
  g_nei <- (g_nei-mean(g_nei))/stats::sd(g_nei)

  #pheno simu
  q <- ncol(geno)
  b_self <- rep(0,q)
  b_self[sample(1:q,n_causal)] = 1

  b_nei <- rep(0,q)
  b_nei[sample(1:q,n_causal/2)] <- 1
  b_nei[sample(1:q,n_causal/2)] <- -1

  K_self <- tcrossprod(geno)
  K_self <- (q/2+K_self/2)/q
  K_nei <- tcrossprod(g_nei)/(q-1)

  b_self <- geno%*%b_self
  b_nei <- g_nei%*%b_nei
  eigenK_self <- eigen(K_self)
  eigenK_nei <- eigen(K_nei)

  pheno <- qtl_pheno_simu(b_self, b_nei, eigenK_self, eigenK_nei, b_ratio = b_ratio, pveB=pveB, pve=pve)
  while(((round(pheno$res_pveB[1],2)==pveB)&(round(pheno$res_pve[1],2)==pve))==FALSE) {
    pheno <- qtl_pheno_simu(b_self, b_nei, eigenK_self, eigenK_nei, b_ratio = b_ratio, pveB=pveB, pve=pve)
  }
  return(pheno$y)
}
