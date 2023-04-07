#' Modified Generalized Estimating Equations for Binary Outcome
#'
#' \code{geessbin} provides results of modified generalized estimating equations
#' with bias-adjusted covariance estimators for longitudinal or clustered
#' data with binary outcomes
#'
#' @param formula Object of class formula: symbolic description of model to be
#'        fitted (see documentation of \code{lm} and
#'        \code{formula} for details).
#' @param data  Data frame.
#' @param id  Vector that identifies the subjects or clusters (NULL by default).
#' @param repeated Vector that identifies repeatedly measured variable within
#'        each subject or cluster. If \code{repeated = NULL}, as is the case in
#'        function \code{gee}, data are assumed to be sorted so that
#'        observations on a cluster are contiguous rows for all entities
#'        in the formula.
#' @param corstr Working correlation structure. The following are permitted:
#'        "\code{independence}", "\code{exchangeable}", "\code{ar1}", and
#'        "\code{unstructured}" ("\code{independence}" by default).
#' @param beta.method Method for estimating regression parameters.
#'        The following are permitted: "\code{GEE}", "\code{PGEE}", and
#'        "\code{BCGEE}" ("\code{PGEE}" by default).
#' @param SE.method Method for estimating standard errors. The following are
#'        permitted: "\code{SA}", "\code{MK}", "\code{KC}", "\code{MD}",
#'        "\code{FG}", "\code{PA}", "\code{GS}", "\code{MB}", "\code{WL}",
#'        "\code{WB}", "\code{FW}", and "\code{FZ}" ("\code{MB}" by default).
#' @param b Numeric vector specifying initial values of regression coefficients.
#'        If \code{b = NULL} (i.e., default value), the initial values are
#'        calculated using the ordinary or Firth logistic regression assuming
#'        that all the observations are independent.
#' @param maxitr Maximum number of iterations (50 by default).
#' @param tol Tolerance used in fitting algorithm (\code{1e-5} by default).
#' @param scale.fix Logical variable; if \code{TRUE}, the scale parameter is
#'        fixed at 1 (\code{FALSE} by default).
#' @param conf.level Numeric value of confidence level for confidence intervals
#'        (0.95 by default).
#'
#' @return an object of class "\code{geessbin}" representing the results of
#' modified generalized estimating equations with bias-adjusted covariance
#' estimators. Generic function \code{summary} provides details of the results.
#'
#' @references Liang, K. and Zeger, S. (1986). Longitudinal data analysis using
#'         generalized linear models.
#'         \emph{Biometrika}, 73, 13-22.
#'
#' @examples
#' data(wheeze)
#'
#' # analysis of PGEE method with Morel et al. covariance estimator
#' res <- geessbin(formula = Wheeze ~ City + Time, id = ID,
#'                 repeated = Time, corstr = "ar1", data = wheeze,
#'                 beta.method = "PGEE", SE.method = "MB")
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

  if (length(unique(y)) != 2) stop("outcome vector is not binary")

  if (!is.numeric(y)) stop("outcome vector must be numeric")

  if (!setequal(unique(y), 0:1)) {
    stop("values of outcome vector must be 0 or 1")
  }

  if (length(corstr) > 1) stop("'corstr' has length > 1")

  if (length(beta.method) > 1) stop("'beta.method' has length > 1")

  if (length(SE.method) > 1) stop("'SE.method' has length > 1")

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
    message("dataset has missing value")

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

      fitted.p <- 1 / (1 + exp(-X %*% b))
    } else {
      res_glm <- glm(formula = formula, family = "binomial", data = dat)
      fitted.p <- res_glm$fitted.values
      b <- res_glm$coefficients
    }

    if(min(fitted.p) < 0.0001 | max(fitted.p) > 0.9999) {
      stop(paste0("error in calculating initial value: ",
                  "fitted probabilities numerically 0 or 1 occurred."))
    }
  } else {
    if (anyNA(b) | sum(is.infinite(b)) > 0) stop("'b' contains Na/NaN/Inf")
    names(b) <- colnames(X)
  }

  conv <- 1
  nitr <- 0
  del <- 100
  while (del > tol) {
    mu <- 1 / (1 + exp(- X %*% b))
    r <- (y - mu) / sqrt(mu * (1 - mu))

    if(min(mu) < 0.0001 | max(mu) > 0.9999) {
      stop(paste0("error in iteration = ", nitr, " :",
                  "fitted probabilities numerically 0 or 1 occurred."))
    }

    phi <- 1
    if (scale.fix == F) phi <- sum(r ^ 2) / (sum(ndat) - p)
    if (is.infinite(phi)) {
      stop(paste0("error in iteration = ", nitr, " :",
                  "infinite scale parameter"))
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
    if (nitr == maxitr) break
  }

  if (del > tol) conv <- 0

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

      Fi <- diag((1 - pmin(0.75, diag(t(mat$VD) %*% mat$D %*% Iinv))) ^ (-0.5),
                 p, p)
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
  if(is.null(digits)) digits <- options()$digits else options(digits =
                                                                digits)
  cat("Call:\n")
  dput(x$call)

  if (x$conv == 1) cstat <- "Success"
  if (x$conv == 0) cstat <- "Failure"

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

  cat("\nConvergence: ", x$convergence, "(", cstat, ")")

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
