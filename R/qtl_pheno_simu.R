#' Simulating phenotype values with neighbor effects.
#'
#' A function to simulate phenotype values with multiple sources of variation controlled
#' @param b_self A n x 1 genotype vector to design major additive genetic effect.
#' @param b_nei A vector of explanatory variable for beta_2 (neighbor effects)
#' @param eigenK_self Products of function 'eigen' with self covariance matrices that were used as explanatory variables for the phenotype.
#' @param eigenK_nei Products of function 'eigen' with neighbor covariance matrices that were used as explanatory variables for the phenotype.
#' @param b_ratio Ratio for contribution of engenK_self,  eigenK_nei and eigenK_sxn to the phenotype.
#' @param pveB Proportion of variance explained by genetic effect designed by b_.. vector.
#' @param pve Proportion of variance explained by all genetic effects (i.e., b_.. and eigenK_..)
#' @return A list of simulated phenotypes
#' \itemize{
#' \item{y}{Simulated phenotype values}
#' \item{beta_self}{major self-genetic effects}
#' \item{beta_nei}{major neighbor effects}
#' \item{sigma_self}{self polygenic effects}
#' \item{sigma_nei}{neighbor polygenic effects}
#' \item{epsilon}{residuals}
#' \item{res_pveB}{realized proportion of variation explained by major-effect genes}
#' \item{res_pve}{realized proportion of variation explained by major-effect genes and polygenic effects}
#' }
#' @author Eiji Yamamoto, and Yasuhiro Sato
qtl_pheno_simu = function(b_self, b_nei, eigenK_self, eigenK_nei, b_ratio = c(1,1), pveB, pve)
{
  beta_self <- b_self / stats::sd(b_self)
  beta_nei <- b_nei / stats::sd(b_nei)
  Maj_eff <- b_ratio[1]*beta_self + b_ratio[2]*beta_nei

  sigma_self <- eigenK_self$vectors %*% stats::rnorm(nrow(eigenK_self$vectors), sd = sqrt(b_ratio[1]*eigenK_self$values))
  sigma_nei <- eigenK_nei$vectors %*% stats::rnorm(nrow(eigenK_nei$vectors), sd = sqrt(b_ratio[2]*eigenK_nei$values))
  sigma_all <- sigma_self + sigma_nei
  sigma_all <- sigma_all / stats::sd(sigma_all)
  adj.sigma = function(adj, Maj_eff, sigma_all, pveB, pve)	{
    (pveB/pve - (stats::var(Maj_eff) / (stats::var(Maj_eff + (sigma_all / adj)))))^2
  }#adj.omega()
  adj <- stats::optimize(adj.sigma , interval=c(0.0000001, 10000000), Maj_eff, sigma_all, pveB, pve)$minimum
  adj.sigma_all <- sigma_all / adj
  Model_effect <- Maj_eff + adj.sigma_all

  epsilon <- stats::rnorm(nrow(b_self), sd = 1)
  pve.adj <- function(adj, Model_effect, epsilon, pve)	{
    (pve - (stats::var(Model_effect) / (stats::var(Model_effect + (epsilon / adj)))))^2
  }#pveM.adj()
  adj <- stats::optimize(pve.adj, interval=c(0.0000001, 10000000),  Model_effect, epsilon, pve)$minimum
  adj.epsilon <- epsilon / adj

  y <- Model_effect + adj.epsilon

  realized_pveB <- stats::var(Maj_eff) / stats::var(y)
  realized_pve <- stats::var(Model_effect) / stats::var(y)

  return(list(y=y, beta_self=beta_self, beta_nei=beta_nei, sigma_self=sigma_self, sigma_nei=sigma_nei, epsilon=epsilon, res_pveB = rep(realized_pveB,length(y)), res_pve = rep(realized_pve,length(y))))
}
