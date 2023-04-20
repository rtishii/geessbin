#' Modified Generalized Estimating Equations for Binary Outcome
#'
#' \code{geessbin} analyzes small-sample clustered or longitudinal data using
#' modified generalized estimating equations (GEE) with bias-adjusted covariance
#' estimator. This function assumes binary outcome and uses the logit link
#' function. This function provides any combination of three GEE methods
#' (conventional and two modified GEE methods) and 12 covariance estimators
#' (unadjusted and 11 bias-adjusted estimators).
#'
#' Details of \code{beta.method} are as follows:
#' \itemize{
#' \item "GEE" is the conventional GEE method (Liang and Zeger, 1986)
#' \item "BCGEE" is the bias-corrected GEE method
#'       (Paul and Zhang, 2014; Lunardon and Scharfstein, 2017)
#' \item "PGEE" is the bias reduction of the GEE method obtained by adding a
#'       Firth-type penalty term to the estimating equation
#'       (Mondol and Rahman, 2019)
#' }
#'
#' Details of \code{SE.method} are as follows:
#' \itemize{
#' \item "SA" is the unadjusted sandwich variance estimator
#'       (Liang and Zeger, 1986)
#' \item "MK" is the MacKinnon and White estimator (MacKinnon and White, 1985)
#' \item "KC" is the Kauermann and Carroll estimator
#'        (Kauermann and Carroll, 2001)
#' \item "MD" is the Mancl and DeRouen estimator (Mancl and DeRouen, 2001)
#' \item "FG" is the Fay and Graubard estimator (Fay and Graubard, 2001)
#' \item "PA" is the Pan estimator (Pan, 2001)
#' \item "GS" is the Gosho et al. estimator (Gosho et al., 2014)
#' \item "MB" is the Morel et al. estimator (Morel et al., 2003)
#' \item "WL" is the Wang and Long estimator (Wang and Long, 2011)
#' \item "WB" is the Westgate and Burchett estimator
#'       (Westgate and Burchett, 2016)
#' \item "FW" is the Ford and Wastgate estimator (Ford and Wastgate, 2017)
#' \item "FZ" is the Fan et al. estimator (Fan et al., 2013)
#' }
#'
#' Descriptions and performances of some of the above methods can be found in
#' Gosho et al. (2023).
#'
#' @param formula Object of class formula: symbolic description of model to be
#'        fitted (see documentation of \code{lm} and
#'        \code{formula} for details).
#' @param data  Data frame.
#' @param id  Vector that identifies the subjects or clusters (\code{NULL} by
#'        default).
#' @param repeated Vector that identifies repeatedly measured variable within
#'        each subject or cluster. If \code{repeated = NULL}, as is the case in
#'        function \code{gee}, data are assumed to be sorted so that
#'        observations on a cluster are contiguous rows for all entities
#'        in the formula.
#' @param corstr Working correlation structure. The following are permitted:
#'        "\code{independence}", "\code{exchangeable}", "\code{ar1}", and
#'        "\code{unstructured}" ("\code{independence}" by default).
#' @param beta.method Method for estimating regression parameters (see Details
#'        section). The following are permitted: "\code{GEE}", "\code{PGEE}",
#'        and "\code{BCGEE}" ("\code{PGEE}" by default).
#' @param SE.method Method for estimating standard errors (see Details section).
#'        The following are permitted: "\code{SA}", "\code{MK}", "\code{KC}",
#'        "\code{MD}", "\code{FG}", "\code{PA}", "\code{GS}", "\code{MB}",
#'        "\code{WL}", "\code{WB}", "\code{FW}", and "\code{FZ}"
#'        ("\code{MB}" by default).
#' @param b Numeric vector specifying initial values of regression coefficients.
#'        If \code{b = NULL} (default value), the initial values are calculated
#'        using the ordinary or Firth logistic regression assuming that all the
#'        observations are independent.
#' @param maxitr Maximum number of iterations (50 by default).
#' @param tol Tolerance used in fitting algorithm (\code{1e-5} by default).
#' @param scale.fix Logical variable; if \code{TRUE}, the scale parameter is
#'        fixed at 1 (\code{FALSE} by default).
#' @param conf.level Numeric value of confidence level for confidence intervals
#'        (0.95 by default).
#'
#' @return The object of class "\code{geessbin}" representing the results of
#' modified generalized estimating equations with bias-adjusted covariance
#' estimators. Generic function \code{summary} provides details of the results.
#'
#' @references \itemize{
#'   \item Fan, C., and Zhang, D., and Zhang, C. H. (2013). A comparison of
#'         bias-corrected covariance estimators for generalized estimating
#'         equations.
#'         \emph{Journal of Biopharmaceutical Statistics}, 23, 1172–1187,
#'         \doi{10.1080/10543406.2013.813521}.\cr
#'   \item Fay, M. P. and Graubard, B. I. (2001). Small-sample adjustments for
#'         Wald-type tests using sandwich estimators.
#'         \emph{Biometrics}, 57, 1198–1206,
#'         \doi{10.1111/j.0006-341X.2001.01198.x}.\cr
#'   \item Ford, W. P. and Westgate, P. M. (2017). Improved standard error
#'         estimator for maintaining the validity of inference in cluster
#'         randomized trials with a small number of clusters.
#'         \emph{Biometrical Journal}, 59, 478–495,
#'         \doi{10.1002/bimj.201600182}.\cr
#'   \item Gosho, M., Ishii, R., Noma, H., and Maruo, K. (2023).
#'         A comparison of bias-adjusted generalized estimating equations for
#'         sparse binary data in small-sample longitudinal studies.
#'         \emph{Statistics in Medicine}, \doi{10.1002/sim.9744}.\cr
#'   \item Gosho, M., Sato, T., and Takeuchi, H. (2014). Robust covariance
#'         estimator for small-sample adjustment in the generalized estimating
#'         equations: A simulation study.
#'         \emph{Science Journal of Applied Mathematics and Statistics},
#'         2, 20–25,
#'         \doi{10.11648/j.sjams.20140201.13}.\cr
#'   \item Kauermann, G. and Carroll, R. J. (2001). A note on the efficiency of
#'         sandwich covariance matrix estimation.
#'         \emph{Journal of the American Statistical Association},
#'         96, 1387–1396,
#'         \doi{10.1198/016214501753382309}.\cr
#'   \item Liang, K. and Zeger, S. (1986). Longitudinal data analysis using
#'         generalized linear models.
#'         \emph{Biometrika}, 73, 13–22,
#'         \doi{10.1093/biomet/73.1.13}.\cr
#'   \item Lunardon, N. and Scharfstein, D. (2017). Comment on ‘Small sample GEE
#'         estimation of regression parameters for longitudinal data’.
#'         \emph{Statistics in Medicine}, 36, 3596–3600,
#'         \doi{10.1002/sim.7366}.\cr
#'   \item MacKinnon, J. G. and White, H. (1985). Some
#'         heteroskedasticity-consistent covariance matrix estimators with
#'         improved finite sample properties.
#'         \emph{Journal of Econometrics}, 29, 305–325,
#'         \doi{10.1016/0304-4076(85)90158-7}.\cr
#'   \item Mancl, L. A. and DeRouen, T. A. (2001). A covariance estimator for
#'         GEE with improved small-sample properties.
#'         \emph{Biometrics}, 57, 126–134,
#'         \doi{10.1111/j.0006-341X.2001.00126.x}.\cr
#'   \item Mondol, M. H. and Rahman, M. S. (2019). Bias-reduced and
#'         separation-proof GEE with small or sparse longitudinal binary data.
#'         \emph{Statistics in Medicine}, 38, 2544–2560,
#'         \doi{10.1002/sim.8126}.\cr
#'   \item Morel, J. G., Bokossa, M. C., and Neerchal, N. K. (2003). Small
#'         sample correlation for the variance of GEE estimators.
#'         \emph{Biometrical Journal}, 45, 395–409,
#'         \doi{10.1002/bimj.200390021}.\cr
#'   \item Pan, W. (2001). On the robust variance estimator in generalised
#'         estimating equations.
#'         \emph{Biometrika}, 88, 901–906,
#'         \doi{10.1093/biomet/88.3.901}.\cr
#'   \item Paul, S. and Zhang, X. (2014). Small sample GEE estimation of
#'         regression parameters for longitudinal data.
#'         \emph{Statistics in Medicine}, 33, 3869–3881,
#'         \doi{10.1002/sim.6198}.\cr
#'   \item Wang, M. and Long, Q. (2011). Modified robust variance estimator for
#'         generalized estimating equations with improved small-sample
#'         performance.
#'         \emph{Statistics in Medicine}, 30, 1278–1291,
#'         \doi{10.1002/sim.4150}.\cr
#'   \item Westgate, P. M. and Burchett, W. W. (2016). Improving power in
#'         small-sample longitudinal studies when using generalized estimating
#'         equations.
#'         \emph{Statistics in Medicine}, 35, 3733–3744,
#'         \doi{10.1002/sim.6967}.
#' }
#'
#' @examples
#' data(wheeze)
#'
#' # analysis of PGEE method with Morel et al. covariance estimator
#' res <- geessbin(formula = Wheeze ~ City + factor(Age), data = wheeze, id = ID,
#'                 corstr = "ar1", repeated = Age, beta.method = "PGEE",
#'                 SE.method = "MB")
#'
#' # hypothesis tests for regression coefficients
#' summary(res)
#'
#' @importFrom MASS ginv
#' @importFrom stats model.matrix model.response model.frame model.extract
#' @importFrom stats glm cov pnorm qnorm
#'
#' @export
geessbin <- function (formula, data = parent.frame(), id = NULL,
                      corstr = "independence", repeated = NULL,
                      beta.method = "PGEE", SE.method = "MB", b = NULL,
                      maxitr = 50, tol = 1e-5, scale.fix = FALSE,
                      conf.level = 0.95)
{

  Call <- match.call()
  mc <- match.call(expand.dots = FALSE)
  mc$corstr <- mc$beta.method <- mc$SE.method <-
    mc$b <- mc$maxitr <- mc$tol <- mc$scale.fix <- mc$conf.level <- NULL
  mc[[1]] <- as.name("model.frame")
  dat <- eval(mc, parent.frame())

  id <- model.extract(dat, "id")
  repeated <- model.extract(dat, "repeated")
  names(id) <- names(repeated) <- NULL

  if (is.null(id) & !is.null(repeated)) {
    stop("'id' must be specified when 'repeated' is not NULL")
  }

  if (is.null(id) & is.null(repeated)) {
    message(paste("'id' and 'repeated' are not specified", "\n",
                  "all observations are assumed to be independent"))

    idseq <- 1:nrow(dat)
    repval <- NULL
    repseq <- rep(1, nrow(dat))
  }

  if (!is.null(id) & is.null(repeated)) {
    id <- deparse(substitute(id))
    idval <- dat[, "(id)"]

    chg <- (1:length(idval))[c(TRUE, idval[-length(idval)] != idval[-1])]
    nidat <- c(chg[-1], length(idval) + 1) - chg
    idseq <- rep(1:length(nidat), time = nidat)
    repseq <- unlist(tapply(nidat, unique(idseq), function(x) 1:x))
    names(repseq) <- repval <- NULL
  }

  if (!is.null(id) & !is.null(repeated)) {
    dat <- dat[order(id, repeated), ]

    idseq <- as.numeric(factor(dat[, "(id)"]))
    repval <- dat[, "(repeated)"]
    repseq <- as.numeric(factor(dat[, "(repeated)"]))
  }

  n <- length(unique(repseq))
  K <- length(unique(idseq))
  ndat <- as.numeric(table(idseq))
  replst <- split(repseq, idseq)

  Terms <- attr(dat, "terms")
  y <- as.matrix(model.extract(dat, "response"))
  X <- model.matrix(Terms, dat)
  p <- ncol(X)

  if (!is.numeric(y) | !setequal(unique(y), 0:1)) {
    stop("outcome vector must be numeric and take values in {0, 1}")
  }

  for(v in c("corstr", "beta.method", "SE.method")){
    if (eval(parse(text = paste0("length(", v, ")"))) > 1) {
      stop(paste0("'", v, "'", " has length > 1"))
    }
  }

  corstrs <- c("independence", "exchangeable", "ar1", "unstructured")
  if (is.na(match(corstr, corstrs))) {
    stop(
      paste(c("invalid correlation structure", "\n",
              "'corstr' must be specified from the following list:", "\n",
              paste(paste0("\"", corstrs, "\""), collapse = ", ")))
    )
  }

  beta.methods <- c("GEE", "BCGEE", "PGEE")
  if (is.na(match(beta.method, beta.methods))) {
    stop(
      paste(c("invalid estimation method", "\n",
              "'beta.method' must be specified from the following list:", "\n",
              paste(paste0("\"", beta.methods, "\""), collapse = ", ")))
    )
  }

  SE.methods <- c("SA", "MK", "KC", "MD", "FG", "PA",
                  "GS", "MB", "WL", "WB", "FW", "FZ")
  if (is.na(match(SE.method, SE.methods))) {
    stop(
      paste(c("invalid SE estimator", "\n",
              "'SE.method' must be specified from the following list:", "\n",
              paste(paste0("\"", SE.methods, "\""), collapse = ", ")))
    )
  }

  if (!is.null(b) & length(b) != p) {
    stop(paste("length of 'b' must be ncol(X) =", p))
  }

  comp <- unlist(lapply(replst, function(x) identical(x, unique(repseq))))
  if (sum(!comp) > 0) {
    if (!is.na(match(SE.method, c("PA", "GS", "WL", "WB")))) {
      stop(paste0("\"", SE.method, "\"",
                  " method cannot be used for incomplete data"))
    }
  }

  if (conf.level <= 0 | conf.level >= 1) {
    stop("'conf.level' must be in interval (0,1)")
  }

  if (is.null(b)) {
    if (beta.method == "PGEE") {
      b <- numeric(p)
      del <- 100
      nitr <- 0
      while (del > 1e-5) {
        mu <- 1 / (1 + exp(-X %*% b))
        I <- t(X) %*% diag(c(mu * (1 - mu))) %*% X
        U <- t(X) %*%
          (y - mu + diag(X %*% ginv(I) %*% t(X)) * mu * (1 - mu) * (0.5 - mu))
        del <- max(abs(U))
        if (del > 1e-5) b <- b + ginv(I) %*% U

        nitr <- nitr + 1
        if (nitr == 50) break
      }
    } else {
      res_glm <- glm(formula = formula, family = "binomial", data = dat)
      b <- res_glm$coefficients
    }
  } else {
    if (anyNA(b) | sum(is.infinite(b)) > 0) stop("'b' contains Na/NaN/Inf")
    names(b) <- colnames(X)
  }

  conv <- "converged"
  nitr <- 0
  del <- 100
  while (del > tol) {
    mu <- 1 / (1 + exp(- X %*% b))
    r <- (y - mu) / sqrt(mu * (1 - mu))

    if (min(mu) < 0.0001 | max(mu) > 0.9999) {
      conv <- "fitted probabilities numerically 0 or 1 occurred."
      warning(conv)
      break
    }

    if (scale.fix == TRUE) phi <- 1
    if (scale.fix == FALSE) phi <- sum(r ^ 2) / (sum(ndat) - p)
    if (is.infinite(phi)) {
      conv <- "infinite scale parameter"
      warning(conv)
      break
    }

    if (corstr == "independence") R <- diag(n)

    if (corstr == "exchangeable") {
      a0 <- 0
      for (i in 1:K) {
        ri <- r[idseq == i]
        pmat <- ri %*% t(ri)
        a0 <- a0 + sum(pmat[upper.tri(pmat)])
      }
      alpha <- a0 / ((0.5 * sum(ndat * (ndat - 1)) - p) * phi)
      R <- matrix(alpha, n, n) + diag(1 - alpha, n, n)
    }

    if (corstr == "ar1") {
      a0 <- d0 <- 0
      for (i in 1:K) {
        ti <- replst[[i]]
        ri <- numeric(n)
        ri[replst[[i]]] <- r[idseq == i]
        a0 <- a0 + sum(ri[-1] * ri[-n])
        d0 <- d0 + sum((ti[-1] - ti[-length(ti)]) == 1)
      }
      alpha <- a0 / ((d0 - p) * phi)
      R <- alpha ^ abs(matrix(0:(n - 1), nrow = n, ncol = n, byrow = TRUE)
                       - 0:(n - 1))
    }

    if (corstr == "unstructured") {
      m <- count <- matrix(0, n, n)
      for (i in 1:K) {
        ri <- ci <- numeric(n)
        ri[replst[[i]]] <- r[idseq == i]
        ci[replst[[i]]] <- 1
        m <- m + t(t(ri)) %*% t(ri)
        count <- count + t(t(ci)) %*% t(ci)
      }
      R <- m / (phi * (count - p))
      diag(R) <- 1
    }

    U <- numeric(p)
    I <- matrix(0, p, p)
    dI <- array(0, c(p, p, p))
    for (i in 1:K) {
      mat <- calc_mat(X[idseq == i, , drop = FALSE], y[idseq == i], b,
                      R[replst[[i]], replst[[i]]], phi)

      U <- U + t(mat$VD) %*% mat$e
      I <- I + t(mat$D) %*% mat$VD

      if (beta.method == "PGEE") {
        dI <- dI + array(apply(X[idseq == i, , drop = FALSE], 2, function (x) {
          t(mat$D) %*% diag(c(1 - 2 * mat$mu) * x, ndat[i], ndat[i]) %*%
            mat$VD
        }), c(p, p, p))
      }
    }

    Iinv <- ginv(I)

    if (beta.method == "PGEE") {
        U <- U + 0.5 * apply(dI, 3, function (x) sum(diag(Iinv %*% x)))
    }

    del <- max(abs(U))

    if (del > tol) b <- b + Iinv %*% U

    nitr <- nitr + 1
    if (nitr == maxitr) {
      if (del > tol) conv <- "maximum number of iterations consumed"
      warning(conv)
      break
    }
  }

  if (conv == "converged" & del > tol) {
    conv <- "convergence failure"
    warning(conv)
  }

  if (conv == "converged") {
    if (beta.method == "BCGEE") {
      k11 <- array(0, c(p, p))
      k21 <- k3 <- array(0, c(p, p, p))
      for (i in 1:K) {
        mat <- calc_mat(X[idseq == i, , drop = FALSE], y[idseq == i], b,
                        R[replst[[i]], replst[[i]]], phi)

        k11 <- k11 + t(mat$VD) %*% mat$emat %*% mat$VD
        dDi <- array(0, c(ndat[i], p, p))
        dDVi <- array(0, c(p, ndat[i], p))
        for (u in 1:p) {
          XuPi <- diag(X[idseq == i, u] * c(1 - 2 * mat$mu), ndat[i], ndat[i])
          dDi[, , u] <- XuPi %*% mat$D
          dDVi[, , u] <- t(mat$D) %*% XuPi %*% mat$Vinv -
            0.5 * t(mat$D) %*% (mat$Vinv %*% XuPi + XuPi %*% mat$Vinv)
        }
        dDVimat <- matrix(dDVi, c(p, ndat[i] * p))

        for (u in 1:p) {
          k21[, , u] <- k21[, , u] +
            dDVimat %*% (diag(p) %x% (mat$emat %*% mat$VD[, u]))
          k3[, , u] <- k3[, , u] -
            (dDVimat %*% (diag(p) %x% mat$D[, u]) +
               dDVi[, , u] %*% mat$D + t(mat$VD) %*% dDi[, , u])
        }
      }

      bhat0 <- numeric(p)
      for (u in 1:p) {
        bhat0 <- bhat0 +
          (k21[, , u] + 0.5 * k3[, , u] %*% Iinv %*% k11) %*% Iinv[, u]
      }

      b <- b - Iinv %*% bhat0

      I <- matrix(0, p, p)
      for (i in 1:K) {
        mat <- calc_mat(X[idseq == i, , drop = FALSE], y[idseq == i], b,
                        R[replst[[i]], replst[[i]]], phi)

        I <- I + t(mat$D) %*% mat$VD
      }

      Iinv <- ginv(I)
    }


    J <- matrix(0, p, p)

    if (SE.method == "SA" | SE.method == "MK") {
      for (i in 1:K) {
        mat <- calc_mat(X[idseq == i, , drop = FALSE], y[idseq == i], b,
                        R[replst[[i]], replst[[i]]], phi)

        J <- J + t(mat$VD) %*% mat$emat %*% mat$VD
      }

      if (SE.method == "MK") J <- J * K / (K - p)
    }

    if (SE.method == "KC") {
      for (i in 1:K) {
        mat <- calc_mat(X[idseq == i, , drop = FALSE], y[idseq == i], b,
                        R[replst[[i]], replst[[i]]], phi)

        HKC <- sqrtmat(ginv(diag(ndat[i]) - mat$D %*% Iinv %*% t(mat$VD)))
        J <- J + t(mat$VD) %*% HKC %*% mat$emat %*% t(HKC) %*% mat$VD
      }
    }

    if (SE.method == "MD") {
      for (i in 1:K) {
        mat <- calc_mat(X[idseq == i, , drop = FALSE], y[idseq == i], b,
                        R[replst[[i]], replst[[i]]], phi)

        HMD <- ginv(diag(ndat[i]) - mat$D %*% Iinv %*% t(mat$VD))
        J <- J + t(mat$VD) %*% HMD %*% mat$emat %*% t(HMD) %*% mat$VD
      }
    }

    if (SE.method == "FG") {
      for (i in 1:K) {
        mat <- calc_mat(X[idseq == i, , drop = FALSE], y[idseq == i], b,
                        R[replst[[i]], replst[[i]]], phi)

        Fi <- diag(
          (1 - pmin(0.75, diag(t(mat$VD) %*% mat$D %*% Iinv))) ^ (-0.5), p, p)
        J <- J + Fi %*% t(mat$VD) %*% mat$emat %*% mat$VD %*% Fi
      }
    }

    if (SE.method == "PA" | SE.method == "GS") {
      M <- matrix(0, n, n)
      for (i in 1:K) {
        mat <- calc_mat(X[idseq == i, , drop = FALSE], y[idseq == i], b,
                        R[replst[[i]], replst[[i]]], phi)

        Aiinv05 <- diag(c(sqrt(1 / mat$nu)), ndat[i], ndat[i])
        M <- M + Aiinv05 %*% mat$emat %*% Aiinv05
      }

      for (i in 1:K) {
        mat <- calc_mat(X[idseq == i, , drop = FALSE], y[idseq == i], b,
                        R[replst[[i]], replst[[i]]], phi)
        Ai05 <- diag(c(sqrt(mat$nu)), ndat[i], ndat[i])

        J <- J + t(mat$VD) %*% Ai05 %*% M %*% Ai05 %*% mat$VD
      }

      if (SE.method == "PA") J <- J / K
      if (SE.method == "GS") J <- J / (K - p)
    }

    if (SE.method == "MB") {
      d <- matrix(0, K, p)
      for (i in 1:K) {
        mat <- calc_mat(X[idseq == i, , drop = FALSE], y[idseq == i], b,
                        R[replst[[i]], replst[[i]]], phi)

        d[i, ] <- t(mat$VD) %*% mat$e
      }

      I1 <- (sum(ndat) - 1) * K * cov(d) / (sum(ndat) - p)
      q <- min(0.5, p / (K - p)) * max(1, sum(diag(Iinv %*% I1)) / p)
      J <- I1 + q * I
    }

    if (SE.method == "WL") {
      M <- matrix(0, n, n)
      for (i in 1:K) {
        mat <- calc_mat(X[idseq == i, , drop = FALSE], y[idseq == i], b,
                        R[replst[[i]], replst[[i]]], phi)
        Aiinv05 <- diag(c(sqrt(1 / mat$nu)), ndat[i], ndat[i])

        HMD <- ginv(diag(ndat[i]) - mat$D %*% Iinv %*% t(mat$VD))
        M <- M + Aiinv05 %*% HMD %*% mat$emat %*% t(HMD) %*% Aiinv05
      }

      for (i in 1:K) {
        mat <- calc_mat(X[idseq == i, , drop = FALSE], y[idseq == i], b,
                        R[replst[[i]], replst[[i]]], phi)
        Ai05 <- diag(c(sqrt(mat$nu)), ndat[i], ndat[i])

        J <- J + t(mat$VD) %*% Ai05 %*% M %*% Ai05 %*% mat$VD / K
      }
    }

    if (SE.method == "WB") {
      M <- matrix(0, n, n)
      for (i in 1:K) {
        mat <- calc_mat(X[idseq == i, , drop = FALSE], y[idseq == i], b,
                        R[replst[[i]], replst[[i]]], phi)
        Aiinv05 <- diag(c(sqrt(1 / mat$nu)), ndat[i], ndat[i])

        HKC <- sqrtmat(ginv(diag(ndat[i]) - mat$D %*% Iinv %*% t(mat$VD)))
        M <- M + Aiinv05 %*% HKC %*% mat$emat %*% t(HKC) %*% Aiinv05
      }

      for (i in 1:K) {
        mat <- calc_mat(X[idseq == i, , drop = FALSE], y[idseq == i], b,
                        R[replst[[i]], replst[[i]]], phi)
        Ai05 <- diag(c(sqrt(mat$nu)), ndat[i], ndat[i])

        J <- J + t(mat$VD) %*% Ai05 %*% M %*% Ai05 %*% mat$VD / K
      }
    }

    if (SE.method == "FW") {
      for (i in 1:K) {
        mat <- calc_mat(X[idseq == i, , drop = FALSE], y[idseq == i], b,
                        R[replst[[i]], replst[[i]]], phi)

        Hi <- mat$D %*% Iinv %*% t(mat$VD)
        HKC <- sqrtmat(ginv(diag(ndat[i]) - Hi))
        HMD <- ginv(diag(ndat[i]) - Hi)
        J <- J + 0.5 * t(mat$VD) %*%
          (HKC %*% mat$emat %*% t(HKC) + HMD %*% mat$emat %*% t(HMD)) %*% mat$VD
      }
    }

    if (SE.method == "FZ") {
      for (i in 1:K) {
        mati <- calc_mat(X[idseq == i, , drop = FALSE], y[idseq == i], b,
                         R[replst[[i]], replst[[i]]], phi)

        HMD <- ginv(diag(ndat[i]) - mati$D %*% Iinv %*% t(mati$VD))

        M <- matrix(0, p, p)
        for (j in setdiff(1:K, i)) {
          matj <- calc_mat(X[idseq == j, , drop = FALSE], y[idseq == j], b,
                           R[replst[[j]], replst[[j]]], phi)

          M <- M + t(matj$VD) %*% matj$emat %*% matj$VD
        }

        J <- J + t(mati$VD) %*% HMD %*%
          (mati$emat - mati$D %*% Iinv %*% M %*% Iinv %*% t(mati$D)) %*%
          t(HMD) %*% mati$VD
      }
    }

    covb <- Iinv %*% J %*% Iinv

  } else {
    covb <- matrix(NA, p, p)
  }

  if (!exists("phi")) phi <- NA

  if (!exists("b")) b <- rep(NA, p)

  if (!exists("R")) R <- matrix(NA, n, n)

  b <- as.vector(b)
  names(b) <- colnames(X)

  if (is.null(repval)) {
    colnames(R) <- rownames(R) <- NULL
  } else {
    rep_unique <- unique(data.frame(repseq = repseq, repval = repval))
    rep_unique <- rep_unique[order(rep_unique$repseq), ]
    colnames(R) <- rownames(R) <- rep_unique$repval
  }

  structure(class = "geessbin",
            list(call = Call,
                 coefficients = b,
                 scale = phi,
                 covb = covb,
                 wcorr = R,
                 iterations = nitr,
                 beta.method = beta.method,
                 SE.method = SE.method,
                 K = K,
                 max.ni = n,
                 corstr = corstr,
                 convergence = conv,
                 conf.level = conf.level,
                 data = data))
}

calc_mat <- function(X, y, b, R, phi) {
  mu <- 1 / (1 + exp( - X %*% b))
  nu <- mu * (1 - mu)
  e <- y - mu
  Ainv05 <- diag(c(sqrt(1 / nu)), nrow(X), nrow(X))
  D <- diag(c(nu), nrow(X), nrow(X)) %*% X
  Vinv <- Ainv05 %*% ginv(R) %*% Ainv05 / phi

  list(mu = mu, nu = nu, D = D,
       Vinv = Vinv, VD = Vinv %*% D,
       e = e, emat = e %*% t(e))
}

#' @export
print.geessbin <- function(x, digits = 3, ...) {
  if(is.null(digits)) digits <- options()$digits else options(digits = digits)

  cat("Call:\n")
  dput(x$call)

  cat("\nCorrelation Structure: ", x$corstr, "\n")
  cat("Estimation Method for Regression Coefficients: ", x$beta.method, "\n")
  cat("Estimation Method for Standard Errors: ", x$SE.method, "\n")

  cat("\nNumber of observations: ", nrow(x$data), "\n")
  cat("Number of clusters: ", x$K, "\n")
  cat("Maximum cluster size: ", x$max.ni, "\n")

  cat("\nCoefficients:\n")
  print(x$coefficients, digits = digits, ...)

  cat("\nEstimated Scale Parameter: ", format(round(x$scale, digits)))
  cat("\nNumber of Iterations: ", x$iterations, "\n")

  cat("\nWorking Correlation:\n")
  print(x$wcorr, digits = digits, ...)

  if (x$convergence == "converged") {
    cat("\nConvergence status: ", "Converged")
  } else {
    cat("\nConvergence status: ", "Failed (", x$convergence, ")\n")
  }

  invisible(x)
}

#' @export
summary.geessbin <- function(object, ...){
  b <- object$coefficients
  se <- sqrt(diag(object$covb))
  coef <- matrix(c(b, se, b/se, pnorm(-abs(b/se), 0, 1) * 2), ncol = 4)
  colnames(coef) <- c("Estimate", "Std.err", "Z", "P.value")
  rownames(coef) <- names(b)

  b0 <- b[toupper(names(b)) != "(INTERCEPT)"]
  if (length(b0) > 0) {
    se0 <- se[toupper(names(b)) != "(INTERCEPT)"]

    OR <- exp(b0)
    q <- qnorm(1 - (1 - object$conf.level) / 2)
    OR <- cbind(matrix(exp(b0), ncol = 1),
                matrix(exp(b0 - q * se0), ncol = 1),
                matrix(exp(b0 + q * se0), ncol = 1))
    colnames(OR) <- c("Odds Ratio", "Lower Limit", "Upper Limit")
    rownames(OR) <- names(b0)
  } else {
    OR <- NULL
  }

  structure(class = "summary.geessbin",
            list(call = object$call,
                 coefficients = coef,
                 OR = OR,
                 scale = object$scale,
                 wcorr = object$wcorr,
                 iterations = object$iterations,
                 beta.method = object$beta.method,
                 SE.method = object$SE.method,
                 corstr = object$corstr,
                 conf.level = object$conf.level))
}

#' @export
print.summary.geessbin <- function(x, digits = 3, ...) {
  if(is.null(digits)) digits <- options()$digits else options(digits =
                                                                digits)
  cat("Call:\n")
  dput(x$call)

  cat("\nCorrelation Structure: ", x$corstr, "\n")
  cat("Estimation Method for Regression Coefficients: ", x$beta.method, "\n")
  cat("Estimation Method for Standard Errors: ", x$SE.method, "\n")

  cat("\nCoefficients:\n")
  print(x$coefficients, digits = digits)

  cat("\nOdds Ratios with", paste0(x$conf.level * 100, "%"),
      "Confidence Intervals", ":\n")
  print(x$OR, digits = digits)

  cat("\nEstimated Scale Parameter: ", format(round(x$scale, digits)))
  cat("\nNumber of Iterations: ", x$iterations, "\n")

  cat("\nWorking Correlation:\n")
  print(x$wcorr, digits = digits)

  invisible(x)
}
