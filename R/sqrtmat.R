#' Square root of nonsymmetric matrix
#'
#' \code{sqrtmat} is used to calculate the square root of \eqn{E_i - H_{ii}},
#' which is an adjustment factor in Kauermann and Carroll-type method.
#'
#' @param M Square matrix whose square root is to be computed.
#'
#' @return The square root of \code{M}
#'
#' @references Kauermann, G. and Carroll, R. J. (2001). A note on the efficiency
#'         of sandwich covariance matrix estimation.
#'         \emph{Journal of the American Statistical Association},
#'         96, 1387–1396,
#'         \doi{10.1198/016214501753382309}.
#'
#' @importFrom MASS ginv
#'
#' @export
sqrtmat <- function(M){
  eig <- eigen(M)
  v <- eig$values
  if (!is.complex(v)) v <- complex(real = v, imaginary = 0)

  return(Re(eig$vectors %*% diag(sqrt(v)) %*% ginv(eig$vectors)))
}
