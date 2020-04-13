#' Calculating neighbor genotypic identity
#'
#' A function to calculate neighbor genotypic identity, with a given reference scale and a degree of distance decay.
#' @param geno An individual x marker matrix. Bialleles (i.e., A or a) must be converted into -1 or 1 digit.
#' @param smap A matrix showing a spatial map for individuals. The first and second column include spatial points along a x-axis and y-axis, respectively.
#' @param scale A numeric scalar indicating the maximum spatial distance between a focal individual and neighbors to define neighbor effects.
#' @param alpha Distance decay coefficient \eqn{\alpha} in a dispersal kernel. Default is set at Inf, meaning no distance decay.
#' @param kernel Type of dispersal kernel in the distance decay. Select \code{"exp"} or \code{"gaussian"} for a negative exponential kernel (fat-tailed) or Gaussian kernel (thin-tailed), respectively.
#' @param grouping A integer vector assigning each individual to a group. This argument can be useful when a "smap" contains different experimental replicates. Default setting means that all individuals are belong to a single group.
#' @param n_core No. of cores for a multi-core computation. This does not work for Windows OS. Default is a single-core computation.
#' @return A numeric matrix for neighbor covariates, with no. of individuals x markers.
#' @author Yasuhiro Sato (\email{sato.yasuhiro.36c@kyoto-u.jp})
#' @references
#' Nathan R, Klein E, Robledo-Arnuncio JJ, Revilla E, (2012) Dispersal kernels: review. In: Clobert J, Baguette M, Benton TG, Bullock JM (Eds.), *Dispersal Ecology and Evolution*. Oxford University Press, pp.186-210.
#' @import parallel
#' @examples
#' set.seed(1)
#' g <- matrix(sample(c(-1,1),100*1000,replace = TRUE),100,1000)
#' gmap <- cbind(c(rep(1,nrow(g)/2),rep(2,nrow(g)/2)),c(1:ncol(g)))
#' x <- runif(nrow(g),1,100)
#' y <- runif(nrow(g),1,100)
#' smap <- cbind(x,y)
#' grouping <- c(rep(1,nrow(g)/2), rep(2,nrow(g)/2))
#'
#' g_nei <- nei_coval(g,smap,44,grouping = grouping)
#' @export
nei_coval = function(geno, smap, scale, alpha=Inf, kernel="exp", grouping=rep(1,nrow(smap)), n_core=1L) {
  p <- nrow(smap)

  coval_i = function(i) {
    id <- c(1:p)[grouping == grouping[i]]

    geno_i <- as.numeric(geno[i,])
    d_i <- mapply(function(x) { return(sqrt((smap[x,1]-smap[i,1])^2 + (smap[x,2]-smap[i,2])^2)) },id)
    j_id <- id[(d_i!=0)&(d_i<=scale)]
    d_i <- d_i[(d_i!=0)&(d_i<=scale)]

    if(length(j_id)==0) {
      return(geno_i*0)
    }

    if(alpha==Inf){
      res <- t(geno[j_id,])
    } else {
      res <- mapply(function(x) { return(w(d_i[x], a=alpha, kernel=kernel)*geno[j_id[x],]) },1:length(j_id))
    }

    if(length(j_id)==1) {
      return(geno_i*res)
    } else {
      return(geno_i*apply(res,1,sum)/length(j_id))
    }
  }
  g_nei <- parallel::mcmapply(coval_i, 1:p, mc.cores=getOption("mc.cores", n_core), mc.preschedule=FALSE)
  g_nei <- t(g_nei)

  colnames(g_nei) <- colnames(geno)
  rownames(g_nei) <- rownames(geno)
  return(g_nei)
}
