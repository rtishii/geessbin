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
#'                        Kingston (0 = Portage; 1 = Kingston).}
#'     \item{\code{Age}}{age of child in years ranging from 9 to 12.}
#'     \item{\code{Smoke}}{measure of smoking habits (cigarettes per day) of
#'                         child's mother.}
#'   }
#'
#' @examples data(wheeze)
#'
#' @references \itemize{
#'   \item Hardin, J. and Hilbe, J. (2013).
#'         \emph{Generalized Estimating Equations, 2nd edition}.
#'         Chapman and Hall, London.\cr
#'   \item Lipsitz, S. R., Fitzmaurice, G. M., Orav, E. J., and Laird, N. M.
#'         (1994). Performance of Generalized Estimating Equations in Practical
#'         Situations.
#'         \emph{Biometrics}, 50, 270–278,
#'         \doi{10.2307/2533218}.
#'   }
NULL
