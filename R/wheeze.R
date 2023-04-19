#' Wheeze dataset
#'
#' The data studied the effect of air pollution on the health of 16 children.
#' The outcome variable was the wheezing status measured consistently four
#' times yearly at ages of 9, 10, 11, and 12 years.
#'
#' @name wheeze
#'
#' @docType data
#'
#' @format A data frame with 64 observations on the following 6 variables:
#'   \describe{
#'     \item{\code{ID}}{child identifier.}
#'     \item{\code{Wheeze}}{binary indicator of wheezing presence.}
#'     \item{\code{City}}{binary indicator of whether the child lives in
#'                        Kingston.}
#'     \item{\code{Age}}{age of child in years ranging from 9 to 12.}
#'     \item{\code{Smoke}}{measure of smoking habits of child's mother.}
#'   }
#'
#' @examples data(wheeze)
#'
#' @references \itemize{
#'   \item Hardin, J. and Hilbe, J. (2013).
#'         \emph{Generalized Estimating Equations, 2nd edition}.
#'         Chapman and Hall, London.
#'   \item Fitzmaurice, G.M., Laird, N.M., and Ware, J.H. (2011).
#'         \emph{Applied Longitudinal Analysis, 2nd edition}. Wiley, New York.
#'   }
NULL
