#' Calculating neighbor genotypic identity
#'
#' A function to calculate neighbor genotypic identity, with a given reference scale and a degree of distance decay.
#' @param geno An individual x marker matrix. Bialleles (i.e., A or a) must be converted into -1 or 1 digit.
#' @param smap A matrix showing a spatial map for individuals. The first and second column include spatial points along an x-axis and y-axis, respectively.
#' @param scale A numeric scalar indicating the maximum spatial distance between a focal individual and neighbors to define neighbor effects.
#' @param alpha An option to set a distance decay coefficient \eqn{\alpha} in a dispersal kernel. Default is set at Inf, meaning no distance decay.
#' @param kernel An option to select either \code{"exp"} or \code{"gaussian"} for a negative exponential kernel or Gaussian kernel, respectively.
#' @param grouping A positive integer vector assigning each individual to a group. This argument can be useful when a "smap" contains different experimental replicates. Default setting means that all individuals are belong to a single group.
#' @param n_core No. of cores for a multi-core computation. This does not work for Windows OS. Default is a single-core computation.
#' @return A numeric matrix for neighbor covariates, with no. of individuals x markers.
#' @author Yasuhiro Sato (\email{sato.yasuhiro.36c@kyoto-u.jp})
#' @details 
#' Default setting is recommended for \code{alpha} and \code{kernel} arguments unless spatial distance decay of neighbor effects needs to be modeled.
#' If \code{alpha} is not \code{Inf}, output variables are weighted by a distance decay from a focal individual to \code{scale}.
#' For the type of dispersal kernel in the distance decay, we can choose a negative exponential or Gaussian kernel as a fat-tailed or thin-tailed distribution, respectively.
#' See Nathan et al. (2012) for detailed characteristics of the two dispersal kernels.
#' @references 
#' Nathan R, Klein E, Robledo-Arnuncio JJ, Revilla E. (2012) Dispersal kernels: review. In: Clobert J, Baguette M, Benton TG, Bullock JM (Eds.), *Dispersal Ecology and Evolution*. Oxford University Press, pp.186-210.
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
nei_coval = function(geno, smap, scale, alpha=Inf, kernel=c("exp","gaussian"), grouping=rep(1,nrow(smap)), n_core=1L) {
  kernel <- match.arg(kernel)

  p <- nrow(smap)

  g.d2 <- lapply(unique(grouping), function(gi, s2) {
    id <- which(grouping == gi)
    d2 <- outer(smap[id,1], smap[id,1], "-")^2 + outer(smap[id,2], smap[id,2], "-")^2
    d2[d2>s2] <- 0
    if (alpha==Inf) {
      wd <- NULL
    } else {
      wd <- w(sqrt(d2), a=alpha, kernel=kernel)
    }
    list(id=id,d2=d2,wd=wd)
  }, scale^2)

  n_div <- ceiling(ncol(geno)/3000/n_core)*n_core
  div.i <- div.seq(ncol(geno), max(n_core, n_div))

  coval_i = function(i1, i2) {
    res <- matrix(0, p, i2-i1+1L)
    for (g.d2_i in g.d2) {
      for (i in seq_along(g.d2_i$id)) {
        id <- g.d2_i$id[i]
        j_id <- g.d2_i$id[g.d2_i$d2[,i] != 0]
        if (length(j_id) > 0) {
          if (alpha==Inf) {
            res_i12 = geno[j_id,i1:i2,drop=F]
          } else {
            w_i <- g.d2_i$wd[g.d2_i$d2[,i] != 0, i]
            res_i12 <- w_i * geno[j_id,i1:i2,drop=F]
          }
          if(length(j_id)==1) {
            res[id,] <- res_i12
          } else {
            res[id,] <- colSums(res_i12, na.rm = T)/colSums(!is.na(res_i12))
          }
        }
      }
    }
    res
  }
  g_nei <- geno*do.call(cbind,
                        parallel::mcmapply(coval_i, div.i$i1, div.i$i2,
                                           mc.cores=getOption("mc.cores", n_core),
                                           USE.NAMES = FALSE, SIMPLIFY = FALSE))
  g_nei[is.na(g_nei)] <- 0

  colnames(g_nei) <- colnames(geno)
  rownames(g_nei) <- rownames(geno)
  return(g_nei)
}

div.seq <- function(n, n.div) {
  if (n.div >= n) return(list(i1 = 1L:n, i2 = 1L:n))
  c.div <- rep(n%/%n.div, n.div)
  if ((r <- n%%n.div) > 0) {
    c.div[1:r] <- c.div[1:r] + 1L
  }
  i2 <- cumsum(c.div)
  list(i1 = c(1L, i2[-length(i2)]+1L), i2 = i2)
}
