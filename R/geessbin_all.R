#' Function for analysis using all combinations of GEE methods and covariance
#' estimators
#'
#' \code{geessbin_all} provides analysis results using all combinations of three
#' GEE methods and 12 covariance estimators.
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
#' @return The list containing two data frames.
#' The first is a table of estimates of regression coefficients, standard
#' errors, z-values, and p-values.
#' The second is a table of odds ratios and confidence intervals.
#'
#' @importFrom MASS ginv
#' @importFrom stats model.matrix model.response model.frame model.extract
#' @importFrom stats glm.fit cov pnorm qnorm binomial
#'
#' @export
geessbin_all <- function (formula, data = parent.frame(), id = NULL,
                          corstr = "independence", repeated = NULL, b = NULL,
                          maxitr = 50, tol = 1e-5, scale.fix = FALSE,
                          conf.level = 0.95)
{
  Call <- match.call()
  mc <- match.call(expand.dots = FALSE)
  mc$corstr <- mc$b <- mc$maxitr <- mc$tol <-
    mc$scale.fix <- mc$conf.level <- NULL
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

  if (length(corstr) > 1) {
    stop("'corstr' has length > 1")
  }

  corstrs <- c("independence", "exchangeable", "ar1", "unstructured")
  if (is.na(match(corstr, corstrs))) {
    stop(
      paste(c("invalid correlation structure", "\n",
              "'corstr' must be specified from the following list:", "\n",
              paste(paste0("\"", corstrs, "\""), collapse = ", ")))
    )
  }

  if (!is.null(b) & length(b) != p) {
    stop(paste("length of 'b' must be ncol(X) =", p))
  }

  if (conf.level <= 0 | conf.level >= 1) {
    stop("'conf.level' must be in interval (0,1)")
  }

  comp <- unlist(lapply(replst, function(x) length(x) == n))
  if (sum(!comp) > 0) {
    stop(paste0("\"PA\", \"GS\", \"WL\", and \"WB\"",
                " methods cannot be used for incomplete data"))
  }

  id <- Call$id
  repeated <- Call$repeated
  beta.methods <- c("GEE", "BCGEE", "PGEE")
  SE.methods <- c("SA", "MK", "KC", "MD", "FG", "PA",
                  "GS", "MB", "WL", "WB", "FW", "FZ")

  out1 <- out2 <- c()
  for (beta in beta.methods) {
    for (se in SE.methods) {
      c <- paste0("try(geessbin(formula = formula, data = data, id = ",
                  deparse(substitute(id)),
                  ", corstr = corstr, repeated = ",
                  deparse(substitute(repeated)),
                  ", beta.method = \"",
                  beta, "\", SE.method = \"",
                  se, "\", b = ", b, ", maxitr = ", maxitr,
                  ", tol = ", tol, ", scale.fix = ", scale.fix,
                  ", conf.level = ", conf.level, "), silent = TRUE)")

      suppressMessages(
        suppressWarnings(
          res1 <- eval(parse(text = c))
        )
      )

      if (!inherits(res1, "try-error")) {
        coef <- summary(res1)$coefficients
        OR <- summary(res1)$OR

        m_coef <- data.frame(beta.method = beta,
                             SE.method = se,
                             Variable = rownames(coef))

        m_OR <- data.frame(beta.method = beta,
                           SE.method = se,
                           Variable = rownames(OR))

        out1 <- rbind(out1, cbind(m_coef, coef))
        out2 <- rbind(out2, cbind(m_OR, OR))
      } else {
        m_coef <- m_OR <- data.frame(beta.method = beta,
                                     SE.method = se,
                                     Variable = NA)
        m_coef$Estimate <- m_coef$Std.err <-
          m_coef$Z <- m_coef$`Pr(>|Z|)` <- NA

        m_OR$`Odds Ratio` <- m_OR$`Lower Limit` <- m_OR$`Upper Limit` <- NA

        out1 <- rbind(out1, m_coef)
        out2 <- rbind(out2, m_OR)
      }
    }
  }

  rownames(out1) <- rownames(out2) <- NULL
  out <- list(coefficients = out1, OR = out2)
  out
}

