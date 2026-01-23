# discretize (numeric) values
# TODO: Consider treating NA and NaN differently.
#       It may not be natural in R since is.na(NaN) is TRUE
#       and NaNs are omitted by na.omit()
# TODO: Discretize date/datetime (Date/POSIXct/POSIXlt/difftime).
#       Fixes will be needed in some documents
#' @importFrom stats quantile
#' @importFrom stats na.omit
.div <- function(x, k, max_levels) {

  if(is.numeric(x)) {
    # default selection of k (10-by-2 rule)
    if(is.null(k)) {
      k <- pmax(2, floor(length(na.omit(x))^log10(2)/2))
    }
    
    if(length(unique(x)) > k) {
      # Set type = 1 to use the empirical distribution function
      qt <- quantile(x, probs = seq(0, 1, length.out = k + 1),
                     na.rm = TRUE, type = 1)
      ret <- cut(x, breaks = unique(qt), include.lowest = TRUE)
    } else {
      ret <- as.factor(x)
    }
  } else {
    ret <- if(length(unique(x)) > max_levels) x else as.factor(x)
  }

  return(ret)
}

# compute entropy from x = table(.)
.entropy <- function(x) {
  stopifnot(all(x >= 0))

  y <- x[x > 0]

  if(length(y) == 0) {
    return(NA_real_)
  }
  
  p <- y / sum(y)
  ent <- sum(-p * log(p))

  return(ent)
}

# returns a list containing the estimated values, and other quantities
# for further computations.
#' @importFrom stats xtabs
.quantile_grid_aprox <- function(x, y, k, type, max_levels, useNA = TRUE) {
  stopifnot(length(x) == length(y))
  n <- length(x)

  # Initial values used when n == 0
  phi <- NA_real_
  kx <- ky <- 0L
  Ixy <- Hx <- Hy <- NA_real_

  if(n > 0) {
    if(!is.factor(x) || !is.factor(y)) {
      kx <- length(unique(x))
      ky <- length(unique(y))
      
      if(type == "chisq") {
        # .div returns non-factor vector is length(unique(.)) > max_levels.
        # In such a case we return NA unless x and y are identical,
        # with length(unique(.)) as k where NAs are counted as one level.
        phi <- if(identical(x, y)) kx else NA_real_
      } else if(type == "KL") {
        # TODO: Consider using NaN when useNA is TRUE.
        #       See "exclude" argument of ?table
        useNA_ent <- if(useNA) "ifany" else "no"
        Hx <- .entropy(table(x, useNA = useNA_ent))
        Hy <- .entropy(table(y, useNA = useNA_ent))
        Ixy <- if(identical(x, y)) Hx else NA_real_
      }
    } else {
      # drop.unused.levels = TRUE is required to avoid marginal probability 0
      nn <- xtabs(~ x + y, addNA = useNA, drop.unused.levels = TRUE)
      nx  <- apply(nn, 1, sum)
      ny  <- apply(nn, 2, sum)

      kx <- length(nx)
      ky <- length(ny)

      if(type == "chisq") {
        phi <- sum(nn^2 / outer(nx, ny))
      } else if(type == "KL") {
        pp <- nn / n
        pq <- nn / outer(nx, ny) * n

        # - 0 * log(0) = NaN and removed with na.rm = TRUE
        Ixy <- sum(pp * log(pq), na.rm = TRUE)
        Hx <- .entropy(nx)
        Hy <- .entropy(ny)
      }
    }
  }
  
  if(type == "chisq") {
    ret <- list(phi = phi, kx = kx, ky = ky)
  } else if(type == "KL") {
    ret <- list(Ixy = Ixy, Hx = Hx, Hy = Hy, kx = kx, ky = ky)
  }

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
