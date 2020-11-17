#' Standard linear models for testing self and neighbor effects
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
#'  \item{\code{p_self}} {p-value for self effects by a likelihood ratio test between a null and standard linear model}
#'  \item{\code{p_nei}} {p-value for neighbor effects by a likelihood ratio test between models with or without neighbor effects}
#' }
#' @details This function is a subset of \code{neiGWAS()}. \code{nei_lm()} gives detailed results when the option \code{model="lm"} is selected in \code{neiGWAS()}.
#' @author Yasuhiro Sato (\email{sato.yasuhiro.36c@kyoto-u.jp})
#' @import Matrix gaston parallel
#' @seealso \code{\link{neiGWAS}}
#' @export
nei_lm = function(geno, g_nei, pheno, addcovar=NULL, response=c("quantitative","binary"), n_core=1L) {
  response <- match.arg(response)
  
  RcppParallel::setThreadOptions(numThreads=n_core)
  
  geno[is.na(geno)] <- 0L
  q <- ncol(geno)
  
  if(is.null(addcovar)) {
    X <- matrix(1, nrow=length(pheno))
  } else {
    X <- cbind(addcovar)
  }
  
  if(response=="quantitative") {
    res00 <- stats::lm(pheno~X)
    LL00 <- stats::logLik(res00)[1]
  } else {
    res00 <- stats::glm(pheno~X, family="binomial")
    LL00 <- stats::logLik(res00)[1]
  }
  
  test_marker_i = function(i) {
    if(response=="quantitative") {
      X0 <- cbind(addcovar, geno[,i])
      X1 <- cbind(X0, g_nei[,i])
      
      LL_self0 <- stats::logLik(stats::lm(pheno~X0))[1]
      p_self <- stats::pchisq(-2*(LL00-LL_self0), 1, lower.tail=FALSE)
      
      LL_nei <- stats::logLik(stats::lm(pheno~X1))[1]
      p_nei <- stats::pchisq(-2*(LL_self0-LL_nei), 1, lower.tail=FALSE)
      
      lmm_nei <- summary(stats::lm(pheno~X1))
      beta <- lmm_nei$coef[-1:0+nrow(lmm_nei$coef),1]
      
      resList <- c(beta, p_self, p_nei)
      return(resList)
    } else { # if(response=="binary") {
      X0 <- cbind(addcovar, geno[,i])
      X1 <- cbind(X0, g_nei[,i])
      
      LL_self0 <- stats::logLik(stats::glm(pheno~X0, family="binomial"))[1]
      p_self <- stats::pchisq(-2*(LL00-LL_self0), 1, lower.tail=FALSE)
      
      LL_nei <- stats::logLik(stats::glm(pheno~X1, family="binomial"))[1]
      p_nei <- stats::pchisq(-2*(LL_self0-LL_nei), 1, lower.tail=FALSE)
      
      lmm_nei <- summary(stats::glm(pheno~X1, family="binomial"))
      beta <- lmm_nei$coef[-1:0+nrow(lmm_nei$coef),1]
      
      resList <- c(beta, p_self, p_nei)
      return(resList)
    }
  }
  results <- do.call(rbind, parallel::mcmapply(test_marker_i, 1:q, mc.cores=getOption("mc.cores",n_core), SIMPLIFY=FALSE, USE.NAMES=FALSE))
  
  colnames(results) <- c("beta_self", "beta_nei", "p_self", "p_nei")
  results <- as.data.frame(results)
  
  return(results)
}
