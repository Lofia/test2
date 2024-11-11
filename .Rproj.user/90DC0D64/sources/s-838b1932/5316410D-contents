#' @title pvt.isoMean
#'
#' @param y y values
#' @param w weights
#'
#' @return estimated gcm
#' @export
#' @import fdrtool
#'
#' @examples
pvt.isoMean = function(y, w)
{
  # Input:	y: measured values in a regression setting
  #		w: weights
  # Output: 	vector containing estimated (isotonic) values

  n = length(y)

  if(n == 1) return(y)
  else{
    ghat = .C("C_isomean",
              as.double(y),
              as.double(w),
              as.integer(n),
              ghat=double(n), PACKAGE="fdrtool")$ghat

    return(ghat)
  }
}
