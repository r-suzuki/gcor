# discretize (numeric) values
# TODO: Consider treating NA and NaN differently.
#       It may not be natural in R since is.na(NaN) is TRUE
#       and NaNs are omitted by na.omit()
# TODO: Add tests for Date/POSIXt, difftime
#' @importFrom stats quantile na.omit
#' @importFrom utils head tail
.div <- function(x, k, max_levels) {
  is_date_dt <- inherits(x, c("Date", "POSIXt"))

  if(inherits(x, "difftime")) x <- as.double(x)
  
  if(is.numeric(x) || is_date_dt) {
    # default selection of k (10-by-2 rule)
    if(is.null(k)) {
      k <- pmax(2, floor(length(na.omit(x))^log10(2)/2))
    }
    
    if(length(unique(x)) > k) {
      # Set type = 1 to use the empirical distribution function
      qt <- quantile(x, probs = seq(0, 1, length.out = k + 1),
                     na.rm = TRUE, type = 1)
      bk <- unique(qt)
      # Set custom labels for Date/POSIXt
      lab <- if(is_date_dt) {
        paste0(
              c("[", rep("(", length(bk) - 2)),
              head(bk, -1),
              ",",
              tail(bk, -1),
              "]")
      } else {
        NULL
      }
      
      ret <- cut(x, breaks = bk,
                 labels = lab, include.lowest = TRUE,
                 right = TRUE) # required for Date/POSIXt
    } else {
      ret <- as.factor(x)
    }
  } else {
    ret <- if(length(unique(x)) > max_levels) x else as.factor(x)
  }

  return(ret)
}

# returns a list containing the estimated values, and other quantities
# for further computations.
#' @importFrom stats xtabs
.mdep_quantile_grid <- function(x, y, k, max_levels, useNA = TRUE) {
  stopifnot(length(x) == length(y))
  n <- length(x)

  if(n == 0) {
    phi <- NA_real_
    kx <- ky <- 0L
  } else if(!is.factor(x) || !is.factor(y)) {
    # .div returns non-factor vector is length(unique(.)) > max_levels.
    # In such a case we return NA unless x and y are identical,
    # with length(unique(.)) as k where NAs are counted as one level.
    kx <- length(unique(x))
    ky <- length(unique(y))

    phi <- if(identical(x, y)) kx else NA_real_
  } else {
    # drop.unused.levels = TRUE is required to avoid marginal probability 0
    nn <- xtabs(~ x + y, addNA = useNA, drop.unused.levels = TRUE)
    nx  <- apply(nn, 1, sum)
    ny  <- apply(nn, 2, sum)

    phi <- sum(nn^2 / outer(nx, ny))
    kx <- length(nx)
    ky <- length(ny)
  }

  ret <- list(estimate = phi, kx = kx, ky = ky)

  return(ret)
}

.vec2df <- function(x, xname) {
  ret <- data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
  names(ret) <- xname
  return(ret)
}

.gen_msg <- function(par, val, lst) {
  if(is.character(val)) val <- paste0('"', val, '"')
  if(is.character(lst)) lst <- paste0('"', lst, '"')

  paste0("Unsupported argument: ", par, " = ", val, "\n",
         "  Should be one of ", paste(lst, collapse = ", "))
}
