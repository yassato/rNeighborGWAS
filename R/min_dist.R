#' Calculating the minimum distance
#'
#' A function to calculate a Euclidian distance including at least one neighbor for all individuals.
#' @param smap A matrix showing a spatial map for individuals. The first and second column include spatial points along a x-axis and y-axis, respectively.
#' @param grouping A integer vector assigning each individual to a group. This argument can be useful when a "smap" contains different experimental replicates. Default setting means that all individuals are belong to a single group.
#' @return Return a scalar of the minimum Euclidian distance that allows all individuals to have at least one neighbor.
#' @author Yasuhiro Sato (\email{sato.yasuhiro.36c@kyoto-u.jp})
#' @examples
#' set.seed(1)
#' g <- matrix(sample(c(-1,1),100*1000,replace = TRUE),100,1000)
#' gmap = cbind(c(rep(1,nrow(g)/2),rep(2,nrow(g)/2)),c(1:ncol(g)))
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
#'
#' min_s <- min_dist(fake_nei$smap, fake_nei$pheno$grouping)
#' @export
min_dist = function(smap, grouping=rep(1,nrow(smap))) {
  min_d_i <- c()
  for(i in 1:nrow(smap)) {
    id <- c(1:nrow(smap))[grouping==grouping[i]]
    d_i <- mapply(function(x) { return(sqrt((smap[x,1]-smap[i,1])^2 + (smap[x,2]-smap[i,2])^2)) },id)
    min_d_i <- c(min_d_i, min(d_i[d_i!=0]))
  }
  return(max(min_d_i)[1])
}
