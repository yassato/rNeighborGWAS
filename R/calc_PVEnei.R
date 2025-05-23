#' Calculating phenotypic variation explained by neighbor effects
#'
#' A function to calculate PVE by neighbor effects for a series of neighbor distance using a mixed model.
#' @param pheno A numeric vector including phenotypes for individuals
#' @param geno An individual x marker matrix. Bialleles (i.e., A or a) must be converted into -1 or 1 digit.
#' @param smap A matrix showing a spatial map for individuals. The first and second column include spatial points along an x-axis and y-axis, respectively.
#' @param scale_seq A numeric vector including a set of the maximum spatial distance between a focal individual and neighbors to define neighbor effects. A scalar is also allowed.
#' @param addcovar An optional matrix including additional non-genetic covariates. It contains no. of individuals x no. of covariates.
#' @param grouping A positive integer vector assigning each individual to a group. This argument can be useful when a "smap" contains different experimental replicates. Default setting means that all individuals are belong to a single group.
#' @param response An option to select if the phenotype is a \code{"quantitative"} trait subject to linear models, or a \code{"binary"} trait subject to logistic models.
#' @param n_core No. of cores for a multi-core computation. This does not work for Windows OS. Default is a single-core computation.
#' @return A numeric matrix including a given spatial scale, PVE by neighbor effects, and p-values.
#' \itemize{
#'  \item scale Maximum neighbor distance given as an argument
#'  \item PVEself Proportion of phenotypic variation explained (PVE) by self effects. RVE is returned when \code{response = "binary"}
#'  \item  PVEnei Proportion of phenotypic variation explained (PVE) by neighbor effects. RVE is returned when \code{response = "binary"}
#'  \item  p-value p-value by a likelihood ratio test between models with or without neighbor effects (when s is not zero); or between a null model and model with self effects alone (when s = 0). NA when \code{response = "binary"}
#' }
#' @details
#' This function uses mixed models via the \code{gaston} package (Perdry & Dandine-Roulland 2020).
#' If \code{"binary"} is selected, \code{logistic.mm.aireml()} is called via the \code{gaston} package.
#' In such a case, \code{PVEnei} below is replaced by the ratio of phenotypic variation explained (RVE) by neighbor effects as RVE_nei = \eqn{\sigma_2^2}/\eqn{\sigma_1^2} and p-values are not provided.
#' @references
#' Perdry H, Dandine-Roulland C. (2020) gaston: Genetic Data Handling (QC, GRM, LD, PCA) & Linear Mixed Models. https://CRAN.R-project.org/package=gaston
#' @author Yasuhiro Sato (\email{sato.yasuhiro.36c@kyoto-u.jp})
#' @import Matrix gaston parallel
#' @examples
#' set.seed(1)
#' g <- matrix(sample(c(-1,1),100*1000,replace = TRUE),100,1000)
#' gmap <- cbind(c(rep(1,nrow(g)/2),rep(2,nrow(g)/2)),c(1:ncol(g)))
#' x <- runif(nrow(g),1,100)
#' y <- runif(nrow(g),1,100)
#' smap <- cbind(x,y)
#' grouping <- c(rep(1,nrow(g)/2),rep(2,nrow(g)/2))
#' pheno <- nei_simu(geno=g,smap=smap,scale=44,grouping=grouping,n_causal=50,pveB=0.4,pve=0.8)
#'
#' fake_nei <- list()
#' fake_nei[[1]] <- g
#' fake_nei[[2]] <- gmap
#' fake_nei[[3]] <- smap
#' fake_nei[[4]] <- data.frame(pheno,grouping)
#' names(fake_nei) <- c("geno","gmap","smap","pheno")
#'
#' min_s <- min_dist(fake_nei$smap, fake_nei$pheno$grouping)
#' scale_seq <- c(min_s, quantile(dist(fake_nei$smap),c(0.2*rep(1:5))))
#'
#' pve_out <- calc_PVEnei(geno=fake_nei$geno, pheno=fake_nei$pheno[,1],
#'                        smap=fake_nei$smap, scale_seq=scale_seq,
#'                        addcovar=as.matrix(fake_nei$pheno$grouping),
#'                        grouping=fake_nei$pheno$grouping)
#' delta_PVE(pve_out)
#' @export
calc_PVEnei = function(pheno, geno, smap, scale_seq, addcovar=NULL, grouping=rep(1,nrow(smap)), response=c("quantitative", "binary"), n_core=1L) {
  scale_seq <- sort(scale_seq)
  response <- match.arg(response)

  geno0 <- geno
  geno[is.na(geno)] <- 0L

  q <- ncol(geno)
  K_self <- tcrossprod(geno)
  K_self <- ((q-1)/2+K_self/2)/(q-1)
  K_self <- as.matrix(Matrix::nearPD(K_self,maxit=10^6)$mat)

  if(is.null(addcovar)) {
    X <- matrix(1, nrow=length(pheno))
  } else {
    X <- cbind(1, addcovar)
  }

  if(response=="quantitative") {
    res0 <- gaston::lmm.aireml(Y=pheno, X=X, K=list(K_self), verbose=FALSE)
    PVEs <- res0$tau[1]/sum(res0$tau, res0$sigma2)
    p_val <- stats::pchisq(2*(res0$logL-res0$logL0), 1, lower.tail=FALSE)
    resList <- c(0, PVEs, 0, p_val)
  } else { # if(response=="binary"){
    res0 <- gaston::logistic.mm.aireml(Y=pheno, X=X, K=list(K_self), verbose=FALSE)
    PVEs <- res0$tau
    p_val <- NA
    resList <- c(0, PVEs, 0, p_val)
  }

  for(s in scale_seq) {
    if(is.numeric(s)==TRUE) { message("scale = ", round(s,3)) }
    g_nei <- nei_coval(geno=geno0, smap=smap, scale=s, alpha=Inf, kernel="exp", grouping=grouping, n_core=n_core)

    K_nei <- tcrossprod(g_nei)/(q-1)
    K_nei <- as.matrix(Matrix::nearPD(K_nei,maxit=10^6)$mat)

    if(response=="quantitative") {
      res <- gaston::lmm.aireml(Y=pheno, X=X, K=list(K_self,K_nei), verbose=FALSE)
	  PVEself <- res$tau[1]/sum(res$tau, res$sigma2)
      PVEnei <- res$tau[2]/sum(res$tau, res$sigma2)
      p_val <- stats::pchisq(2*(res$logL-res0$logL), 1, lower.tail=FALSE)
    } else {
      res <- gaston::logistic.mm.aireml(Y=pheno, X=X, K=list(K_self,K_nei), verbose=FALSE)
	  PVEself <- res$tau[1]/res$tau[2]
      PVEnei <- res$tau[2]/res$tau[1]
      p_val <- NA
    }
    resList <- rbind(resList, c(s, PVEself, PVEnei, p_val))
  }
  colnames(resList) <- c("scale", "PVEself", "PVEnei", "p_val")
  resList <- as.data.frame(resList)
  rownames(resList) <- NULL
  return(resList)
}
