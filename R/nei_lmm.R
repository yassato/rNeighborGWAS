#' Mixed models for testing self and neighbor effects
#'
#' A function to provide coefficients and p-values of self and neighbor effects for each marker.
#' @param geno An individual x marker matrix. Bialleles (i.e., A or a) must be converted into -1 or 1 digit.
#' @param g_nei An output of \code{nei_coval()} object, namely an individual x marker matrix including neighbor genotypic identity.
#' @param pheno A numeric vector including phenotypes for individuals
#' @param addcovar An optional matrix including additional non-genetic covariates. It contains no. of individuals x no. of covariates.
#' @param response An option to select if the phenotype is a \code{"quantitative"} trait subject to linear models, or a \code{"binary"} trait subject to logistic models.
#' @param n_core No. of cores for a multi-core computation. This does not work for Windows OS. Default is a single-core computation.
#' @return A data.frame including coefficients and p-values of self and neighbor effects, without the chromosome numbers and marker position.
#' \itemize{
#'  \item{\code{beta_self}} {coefficient for self effects}
#'  \item{\code{beta_self}} {coefficient for neighbor effects}
#'  \item{\code{p_self}} {p-value for self effects by a likelihood ratio test between a null and standard GWAS model}
#'  \item{\code{p_nei}} {p-value for neighbor effects by a likelihood ratio test between models with or without neighbor effects}
#' }
#' @details This function is a subset of \code{neiGWAS()}. \code{nei_lmm()} gives detailed results but requires more computational time.
#' @author Yasuhiro Sato (\email{sato.yasuhiro.36c@kyoto-u.jp})
#' @import Matrix gaston parallel
#' @seealso \code{\link{neiGWAS}}
#' @export
nei_lmm = function(geno, g_nei, pheno, addcovar=NULL, response=c("quantitative","binary"), n_core=1L) {
  response <- match.arg(response)

  RcppParallel::setThreadOptions(numThreads=n_core)

  geno[is.na(geno)] <- 0L

  q <- ncol(geno)
  K_self <- tcrossprod(geno)
  K_self <- ((q-1)/2+K_self/2)/(q-1)
  K_nei <- tcrossprod(g_nei)/(q-1)

  K_self <- as.matrix(Matrix::nearPD(K_self,maxit=10^6)$mat)
  K_nei <- as.matrix(Matrix::nearPD(K_nei,maxit=10^6)$mat)

  if(is.null(addcovar)) {
    X <- matrix(1, nrow=length(pheno))
  } else {
    X <- cbind(1, addcovar)
  }

  eiKs <- eigen(K_self)
  res00 <- gaston::lmm.diago(Y=pheno, X=X, eigenK=eiKs, verbose=FALSE)
  LL00 <- gaston::lmm.diago.profile.likelihood(tau=res00$tau, s2=res00$sigma2, Y=pheno, X=X, eigenK=eiKs)[1,1]

  aireml <- gaston::lmm.aireml(Y=pheno, X=X, K=list(K_self,K_nei), verbose=FALSE)
  PVEself <- aireml$tau[1]/sum(aireml$tau, aireml$sigma2)
  PVEnei <- aireml$tau[2]/sum(aireml$tau, aireml$sigma2)
  message("PVE_self = ", PVEself)
  message("PVE_nei = ", PVEnei)

  K <- aireml$tau[1]*K_self + aireml$tau[2]*K_nei
  eiK <- eigen(K)
  res01 <- gaston::lmm.diago(Y=pheno, X=X, eigenK=eiK, verbose=FALSE)

  test_marker_i = function(i) {
    if(response=="quantitative") {
      X0 <- cbind(1, addcovar, geno[,i])
      X1 <- cbind(X0, g_nei[,i])

      LL_self0 <- gaston::lmm.diago.profile.likelihood(tau=res00$tau, s2=res00$sigma2, Y=pheno, X=X0, eigenK=eiKs)[1,1]
      p_self <- stats::pchisq(-2*(LL00-LL_self0), 1, lower.tail=FALSE)

      LL_self1 <- gaston::lmm.diago.profile.likelihood(tau=res01$tau, s2=res01$sigma2, Y=pheno, X=X0, eigenK=eiK)[1,1]
      LL_nei <- gaston::lmm.diago.profile.likelihood(tau=res01$tau, s2=res01$sigma2, Y=pheno, X=X1, eigenK=eiK)[1,1]
      p_nei <- stats::pchisq(-2*(LL_self1-LL_nei), 1, lower.tail=FALSE)

      lmm_nei <- gaston::lmm.diago(Y=pheno, X=X1, eigenK=eiK, verbose=FALSE)
      beta <- lmm_nei$BLUP_beta[-1:0+length(lmm_nei$BLUP_beta)]

      resList <- c(beta, p_self, p_nei)
      return(resList)
    } else { # if(response=="binary") {
      res <- gaston::logistic.mm.aireml(Y=pheno, X=cbind(1,addcovar,geno[,i],g_nei[,i]), K=list(K_self,K_nei), verbose=FALSE)
      beta <- res$BLUP_beta[-1:0+length(res$BLUP_beta)]
      z_val <- res$BLUP_beta/sqrt(diag(res$varbeta))
      p_val <- stats::pchisq(z_val^2, 1, lower.tail=FALSE)[c((length(z_val)-1), length(z_val))]

      resList <- c(beta, p_val)
      return(resList)
    }
  }
  results <- do.call(rbind, parallel::mcmapply(test_marker_i, 1:q, mc.cores=getOption("mc.cores",n_core), SIMPLIFY=FALSE, USE.NAMES=FALSE))

  colnames(results) <- c("beta_self", "beta_nei", "p_self", "p_nei")
  results <- as.data.frame(results)

  return(results)
}

