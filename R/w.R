#' Calculating a distance decay weight
#'
#' A function to calculate, with a negative exponential or Gaussian dispersal kernel.
#'
#' @param s A numeric scalar indicating spatial distance at which the distance decay is referred
#' @param a A numeric scalar indicating a decay coefficient
#' @param kernel An option to select a negative exponential kernel \code{"exp"} or Gaussian kernel \code{"gaussian"}.
#' @return A numeric scalar for a distance decay weight.
#' @author Yasuhiro Sato (\email{sato.yasuhiro.36c@kyoto-u.jp})
w = function(s,a,kernel="exp") {
  if(kernel=="exp") {
    return(exp(-(s/a)))
  } else if(kernel=="gaussian") {
    return(exp(-(s^2/a^2)))
  } else {
    warning("error: Select 'exp' or 'gaussian' for the dispersal kernel")
    return(NULL)
  }
}
