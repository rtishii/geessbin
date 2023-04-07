#' Square root of nonsymmetric matrix
#'
#' \code{sqrtmat} is used to calculate the square root of \eqn{E_i - H_{ii}},
#' which is an adjustment factor in Kauermann and Carroll-type method.
#'
#' @param M Square matrix whose square root is to be computed.
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
