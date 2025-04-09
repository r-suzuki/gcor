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
.mdep_quantile_grid <- function(x, y, k, na.rm = FALSE) {
  stopifnot(length(x) == length(y))

  if(length(x) == 0) {
    phi <- rep(NA_real_, length(k))
    kx <- ky <- rep(0, length(k))
  } else {
    phi <- kx <- ky <- numeric(length(k))

    for(i in seq_along(k)) {
      xx <- if(is.numeric(x)) .div(x, k[i]) else x
      yy <- if(is.numeric(y)) .div(y, k[i]) else y

      nn <- xtabs(~ xx + yy, addNA = !na.rm)
      nx  <- apply(nn, 1, sum)
      ny  <- apply(nn, 2, sum)

      # remove rows with p(x) = 0, and cols with p(y) = 0
      nn <- nn[nx > 0, ny > 0, drop =FALSE]
      nx <- nx[nx > 0]
      ny <- ny[ny > 0]

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

.gen_msg <- function(par, val, lst, pmatch = FALSE) {
  paste0('Unsupported argument: ', par, ' = "', val, '"\n',
         'Should be ',
         (if(pmatch) '(an abbreviation of) ' else ''),
         'one of "', paste(lst, collapse = '", "'), '"')
}
