#' Estimate mutual dependency and related measures
#'
#' @description Estimate generalized correlation measure (`gcor`),
#' generalized dissimilarity between variables (`gdis`), and
#' predictability score (`pscore`) based on mutual dependency.
#'
#' @param x a vector, matrix, data frame or formula. If formula, `data` should be specified.
#' `gdis` requires a matrix or data frame.
#' @param y `NULL` (default) or a vector, matrix or data frame with compatible dimensions to `x`.
#' @param k `NULL` (default) or an integer specifying the number of groups for discretization.
#' Numerical data are divided into `k` groups using `k`-quantiles.
#' If `NULL`, it is determined automatically.
#' @param data `NULL` (default) or a data frame. Required if `x` is a formula.
#' @param drop a logical. If `TRUE`, the returned value is coerced to
#' a vector when one of its dimensions is one.
#' @param ... additional arguments (`diag` and `upper`) passed to `as.dist` function.
#' See \link{as.dist} for details.
#'
#' @return For `gcor` and `pscore`, a numeric matrix is returned (or a vector if `drop = TRUE`).
#' For `gdis`, an object of class `"dist"` is returned.
#'
#' @references
#' Suzuki, R. (2025). *A generalization of correlation coefficient*. preprint.
#' \url{https://r-suzuki.github.io/preprints}
#'
#' @examples
#' # Generalized correlation measure
#' gcor(iris)
#'
#' # Clustering
#' gd <- gdis(iris)
#' hc <- hclust(gd, method = "ward.D2")
#' plot(hc)
#'
#' # Multidimensional scaling
#' mds <- cmdscale(gd, k = 2)
#' plot(mds, type = "n", xlab = "", ylab = "", asp = 1, axes = FALSE,
#'      main = "cmdscale with gdis(iris)")
#' text(mds[,1], mds[,2], rownames(mds))
#'
#' # Predictability of Species from other variables
#' ps <- pscore(Species ~ ., data = iris)
#' dotchart(sort(ps), xlim = c(0, 1), main = "Predictability of Species")
#' @name mdep-package
#' @docType _PACKAGE
#' @aliases mdep
NULL

# @param measure a character specifying the type of measure, one of `"cor"`, `"dist"`, `"pred"`.
# `gcor` is a wrapper for `mdep` with `measure = "cor"`.
# Similarly, `gdis` wraps `measure = "dist"`, and `pscore` wraps `measure = "pred"`.
# @param xname a character to be used as the name of `x`, when x is an atomic vector.
# @param yname a character used as the name of `y` (same as `xname` for `x`).
mdep <- function(x, y = NULL, k = NULL, data = NULL, drop = FALSE, measure,
                 xname = deparse1(substitute(x)), yname = deparse1(substitute(y)),
                 ...
                 ) {
  IS_XY_SYNMETRIC <- FALSE
  MEASURES <- c("cor", "dist", "pred")
  xx <- yy <- kk <- ret <- NULL

  if(is.na(match(measure, MEASURES))) {
    stop(.gen_msg("measure", measure, MEASURES))
  }

  if(is.null(y) & is.atomic(x) && is.null(dim(x))) {
    stop("Supply non-NULL y if x is a non-matrix atomic vector.")
  }

  if(inherits(x, "formula")) {
    if(!is.null(y)) {
      stop("y should be NULL if x is a formula")
    }

    if(is.null(data)) {
      stop("Supply non-null data if x is a formula")
    }

    mf <- model.frame(formula = x, data = data)
    xx <- mf[, -1, drop = FALSE]
    yy <- as.data.frame(model.response(mf))

    if(ncol(yy) == 1) {
      names(yy) <- all.vars(x[[2]])
    }
  }

  if(is.null(xx)) {
    if(is.atomic(x) && is.null(dim(x))) {
      xx <- .vec2df(x, xname)
    } else {
      xx <- as.data.frame(x)
    }
  }

  if(is.null(yy)) {
    if(is.null(y)) {
      yy <- xx
      IS_XY_SYNMETRIC <- TRUE
    } else if(is.atomic(y) && is.null(dim(y))) {
      yy <- .vec2df(y, yname)
    } else {
      yy <- as.data.frame(y)
    }
  }

  stopifnot(nrow(xx) == nrow(yy))

  if(is.null(k)) k <- pmax(2, floor(sqrt(nrow(xx) / 10)))
  stopifnot(length(k) == 1)

  ret <- matrix(rep(NA_real_, ncol(xx) * ncol(yy)),
                nrow = ncol(xx), ncol = ncol(yy),
                dimnames = list(names(xx), names(yy)))

  for(i in seq_len(ncol(xx))) {
    for(j in seq_len(ncol(yy))) {
      if(IS_XY_SYNMETRIC && i > j) {
        # DO NOTHING
      } else if(IS_XY_SYNMETRIC && i == j) {
        ret[i, j] <- if(measure == "dist") 0.0 else 1.0
      } else {
        m_ij <- .mdep_quantile_grid(xx[,i], yy[,j], k)

        kk <- sqrt(m_ij$kx) * sqrt(m_ij$ky)
        r2 <- 1 - 1/m_ij$estimate
        if(measure == "cor") {
          ret[i, j] <- sqrt(r2 / (1 - 1/kk))
          if(IS_XY_SYNMETRIC) ret[j, i]  <- ret[i, j]
        } else if(measure == "dist") {
          ret[i, j] <- sqrt(1 - r2 / (1 - 1/kk))
          if(IS_XY_SYNMETRIC) ret[j, i] <- ret[i, j]
        } else if(measure == "pred") {
          ret[i, j] <- sqrt(r2 / (1 - 1/m_ij$ky))
          if(IS_XY_SYNMETRIC) ret[j, i] <- sqrt(r2 / (1 - 1/m_ij$kx))
        }
      }
    }
  }

  if(measure == "dist") {
    ret <- as.dist(ret, ...)
  } else if(drop) {
    if(nrow(ret) == 1 && ncol(ret) == 1) {
      ret <- as.vector(ret)
    } else if(nrow(ret) == 1) {
      ret <- setNames(as.vector(ret), colnames(ret))
    } else if(ncol(ret) == 1) {
      ret <- setNames(as.vector(ret), rownames(ret))
    }
  }

  return(ret)
}

#' @rdname mdep-package
#' @export
gcor <- function(x, y = NULL, k = NULL, data = NULL, drop = TRUE) {
  mdep(x = x, y = y, k = k, data = data, drop = drop, measure = "cor",
       xname = deparse1(substitute(x)), yname = deparse1(substitute(y)))
}

#' @rdname mdep-package
#' @export
pscore <- function(x, y = NULL, k = NULL, data = NULL, drop = TRUE) {
  mdep(x = x, y = y, k = k, data = data, drop = drop, measure = "pred",
       xname = deparse1(substitute(x)), yname = deparse1(substitute(y)))
}

#' @rdname mdep-package
#' @export
gdis <- function(x, k = NULL, ...) {
  if(!is.matrix(x) && !is.data.frame(x)) {
    stop("x should be a matrix or data frame.")
  }

  mdep(x = x, y = NULL, k = k, data = data, measure = "dist",
       xname = deparse1(substitute(x)), yname = deparse1(substitute(y)), ...)
}
