.div <- function(x, k) {
  if(length(unique(x)) <= k) {
    ret <- as.factor(x)
  } else {
    qt <- quantile(x, probs = seq(0, 1, length.out = k + 1), na.rm = TRUE)
    ret <- cut(x, breaks = unique(qt), include.lowest = TRUE)
  }
  return(ret)
}

# returns a list containing the estimated values, and other quantities for further computations.
.mudep_quantile_grid <- function(x, y, k, useNA = TRUE) {
  stopifnot(length(x) == length(y))

  if(length(x) == 0) {
    phi <- rep(NA_real_, length(k))
    kx <- ky <- rep(0, length(k))
  } else {
    phi <- kx <- ky <- numeric(length(k))

    for(i in seq_along(k)) {
      xx <- if(is.numeric(x)) .div(x, k[i]) else x
      yy <- if(is.numeric(y)) .div(y, k[i]) else y

      # drop.unused.levels = TRUE is required to avoid marginal probability 0
      nn <- xtabs(~ xx + yy, addNA = useNA, drop.unused.levels = TRUE)
      nx  <- apply(nn, 1, sum)
      ny  <- apply(nn, 2, sum)

      phi[i] <- sum(nn^2 / outer(nx, ny))
      kx[i] <- length(nx)
      ky[i] <- length(ny)
    }
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
