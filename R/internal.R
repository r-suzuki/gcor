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
.mdep_quantile_grid <- function(x, y, k, includeNA = TRUE) {
  stopifnot(length(x) == length(y))
  stopifnot(length(x) > 0)

  phi <- kx <- ky <- numeric(length(k))
  n <- length(x)

  for(i in seq_along(k)) {
    xx <- if(is.numeric(x)) .div(x, k[i]) else x
    yy <- if(is.numeric(y)) .div(y, k[i]) else y

    nn <- xtabs(~ xx + yy, addNA = includeNA)

    nx  <- apply(nn, 1, sum)
    ny  <- apply(nn, 2, sum)

    phi[i] <- sum(nn^2 / outer(nx, ny))
    kx[i] <- length(nx)
    ky[i] <- length(ny)
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
  paste0("Unsupported ", par, ": ", val, "\n",
         "Should be one of ", paste(lst, collapse = ", "))
}
