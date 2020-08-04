#' Estimating the effective scale of neighbor effects
#'
#' A function to calculate \eqn{\Delta}PVE that estimates the effective scale of neighbor effects.
#' @param res Output results of \code{calc_PVEnei()}.
#' @param fig TRUE/FALSE to plot the results (or not). Default is TRUE.
#' @param ... Arguments to be passed to \code{plot()}.
#' @return Estimated effective scale and proportion of phenotypic variation explained by neighbor effects at that scale.
#' @author Yasuhiro Sato (\email{sato.yasuhiro.36c@kyoto-u.jp})
#' @seealso \code{\link{calc_PVEnei}}
#' @export
delta_PVE = function(res, fig=TRUE, ...) {
  res <- res[order(res$scale),]
  pve <- res$PVEnei[-1]
  s <- res$scale[-1]
  delta_pve <- pve - c(0,pve[1:(length(pve)-1)])

  est_scale <- s[delta_pve==max(delta_pve)][1]
  est_pve <- pve[delta_pve==max(delta_pve)][1]

  if(fig==TRUE) {
    args <- list(...)
    args$x <- s
    args$y <- pve
    args$type <- "l"
    args$xlab <- "scale"
    args$ylab <- "PVE_nei"
    do.call(graphics::plot,args)
    graphics::points(s,pve)
    graphics::points(est_scale,est_pve,pch=16)
  }

  est <- c(est_scale,est_pve)
  names(est) <- c("scale","PVEnei")
  return(est)
}
