#' Genome-wide association mapping of neighbor effects
#'
#' A function to test neighbor effects for each marker and to calculate p-values at a given reference scale.
#' @param geno An individual x marker matrix. Bialleles (i.e., A or a) must be converted into -1 or 1 digit.
#' @param pheno A numeric vector including phenotypes for individuals
#' @param gmap A matrix or data.frame including chromosome numbers in the first column, and SNP positions in the second column.
#' @param smap A matrix showing a spatial map for individuals. The first and second column include spatial points along an x-axis and y-axis, respectively.
#' @param scale A numeric scalar indicating the maximum spatial distance between a focal individual and neighbors to define neighbor effects.
#' @param addcovar An optional matrix including additional non-genetic covariates. It contains no. of individuals x no. of covariates.
#' @param grouping A positive integer vector assigning each individual to a group. This argument can be useful when a "smap" contains different experimental replicates. Default setting means that all individuals are belong to a single group.
#' @param response An option to select if the phenotype is a \code{"quantitative"} trait subject to linear models, or a \code{"binary"} trait subject to logistic models.
#' @param model An option to select linear mixed model \code{"lmm"} or linear model \code{"lm"}. Default setting is to use a mixed model.
#' @param n_core No. of cores for a multi-core computation. This does not work for Windows OS. Default is a single-core computation.
#' @return A data.frame including the chromosome number, marker position, and p-values.
#' \itemize{
#'  \item{\code{chr}} {Chromosome number}
#'  \item{\code{pos}} {Marker position}
#'  \item{\code{p}} {p-value by a likelihood ratio test between models with or without neighbor effects}
#' }
#' @details
#' This function calls a mixed model via the \code{gaston} package. If \code{"lmm"} with \code{"binary"} is selected, p-values are based on Wald tests.
#' This is because the logistic mixed model is based on a pseudo-likelihood and thus likelihood ratio tests are not applicable. See Chen et al. (2016) for the theory.
#' @author Yasuhiro Sato (\email{sato.yasuhiro.36c@kyoto-u.jp})
#' @references
#' Chen H, Wang C, Conomos M. et al. (2016) Control for population structure and relatedness for binary traits in genetic association studies via logistic mixed models. The American Journal of Human Genetics 98: 653-666.
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
#' scale <- 43
#' gwas_out <- neiGWAS(geno=fake_nei$geno, pheno=fake_nei$pheno[,1],
#'                     gmap=fake_nei$gmap, smap=fake_nei$smap,
#'                     scale=scale, addcovar=as.matrix(fake_nei$pheno$grouping),
#'                     grouping=fake_nei$pheno$grouping)
#'
#' gaston::manhattan(gwas_out)
#' gaston::qqplot.pvalues(gwas_out$p)
#' @export
neiGWAS = function(geno, pheno, gmap, smap, scale, addcovar=NULL, grouping, response=c("quantitative","binary"), model=c("lmm","lm"), n_core=1L) {
  response <- match.arg(response)
  model <- match.arg(model)

  if(ncol(geno)!=nrow(gmap)) {
    warning("error: no. of SNPs does not match between 'geno' and 'gmap'")
    return(NULL)
    }

  RcppParallel::setThreadOptions(numThreads=n_core)
  g_nei <- nei_coval(geno=geno, smap=smap, scale=scale, alpha=Inf, kernel="exp", grouping=grouping, n_core=n_core)

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

  if(response=="quantitative") {
    res0 <- gaston::lmm.aireml(Y=pheno, X=X, K=list(K_self,K_nei), verbose=FALSE)
    K <- res0$tau[1]*K_self + res0$tau[2]*K_nei
    eiK <- eigen(K)
    res0 <- gaston::lmm.diago(Y=pheno, X=X, eigenK=eiK, verbose=FALSE)
  }

  test_marker_i = function(i) {
    if(model=="lmm") {
      if(response=="quantitative") {
        LL_self <- gaston::lmm.diago.profile.likelihood(tau=res0$tau, s2=res0$sigma2, Y=pheno, X=cbind(1,addcovar,geno[,i]), eigenK=eiK)[1,1]
        LL_nei <- gaston::lmm.diago.profile.likelihood(tau=res0$tau, s2=res0$sigma2, Y=pheno, X=cbind(1,addcovar,geno[,i],g_nei[,i]), eigenK=eiK)[1,1]
        p_nei <- stats::pchisq(-2*(LL_self-LL_nei), 1, lower.tail=FALSE)
        return(p_nei)
      } else { # if(response=="binary") {
        res <- gaston::logistic.mm.aireml(Y=pheno, X=cbind(1,addcovar,geno[,i],g_nei[,i]), K=list(K_self,K_nei), verbose=FALSE)
        z_val <- res$BLUP_beta[length(res$BLUP_beta)]/sqrt(res$varbeta[length(res$BLUP_beta), length(res$BLUP_beta)])
        p_nei <- stats::pchisq(z_val^2, 1, lower.tail=FALSE)
        return(p_nei)
      }
    } else {
      if(response=="quantitative") {
        X <- cbind(1, addcovar, geno[,i])
        LL_self_lm <- stats::logLik(stats::lm(pheno~X-1))
        X <- cbind(1, addcovar, geno[,i], g_nei[,i])
        LL_nei_lm <- stats::logLik(stats::lm(pheno~X-1))
        p_nei <- stats::pchisq(-2*(LL_self_lm-LL_nei_lm), 1, lower.tail=FALSE)
        return(p_nei)
      } else { # if(response=="binary") {
        X <- cbind(1, addcovar, geno[,i])
        LL_self_lm <- stats::logLik(stats::glm(pheno~X-1,family="binomial"))
        X <- cbind(1, addcovar, geno[,i], g_nei[,i])
        LL_nei_lm <- stats::logLik(stats::glm(pheno~X-1,family="binomial"))
        p_nei <- stats::pchisq(-2*(LL_self_lm-LL_nei_lm), 1, lower.tail=FALSE)
        return(p_nei)
      }
    }
  }
  results <- unlist(parallel::mcmapply(test_marker_i, 1:q, mc.cores=getOption("mc.cores",n_core), SIMPLIFY=FALSE, USE.NAMES=FALSE))
  results <- data.frame(gmap, results)
  colnames(results) <- c("chr", "pos", "p")
  return(results)
}
