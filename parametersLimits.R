PCMParamLowerLimit.BM <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))
  
  if(is.Global(o$X0)) {
    o$X0[1:6] <- 0.1
  }
  
  if(is.Global(o$Sigma_x)) {
    diag(o$Sigma_x) <- 0.0001
    if(!is.Diagonal(o$Sigma_x)) {
      diag(o$Sigma_x) <- 0.0001
      o$Sigma_x[upper.tri(o$Sigma_x)] <- -1.0
    }
  } else {
    for(r in seq_len(R)) {
      diag(o$Sigma_x[, , r]) <- 0.0001
      if(!is.Diagonal(o$Sigma_x)) {
        diag(o$Sigma_x[, , r]) <- 0.0001
        o$Sigma_x[, , r][upper.tri(o$Sigma_x[, , r])] <- -1.0
      }
    }
  }
  
  o
}

PCMParamUpperLimit.BM <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))
  
  if(is.Global(o$X0)) {
    o$X0[1:6] <- 4
  }
  
  if(is.Global(o$Sigma_x)) {
    diag(o$Sigma_x) <- 2
    if(!is.Diagonal(o$Sigma_x)) {
      o$Sigma_x[upper.tri(o$Sigma_x,diag = TRUE)] <- 2
    }
  } else {
    for(r in seq_len(R)) {
      diag(o$Sigma_x[, , r]) <- 2
      if(!is.Diagonal(o$Sigma_x)) {
        o$Sigma_x[, , r][upper.tri(o$Sigma_x[, , r], diag = TRUE)] <- 2
      }
    }
  }
  o
}

PCMParamLowerLimit.BMkappa <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))
  
  if(is.Global(o$X0)) {
    o$X0[1:6] <- 0.1
  }
  
  if(is.ScalarParameter(o$kappa)) {
    o$kappa[1] <- 0.1
  }
  o
}

PCMParamUpperLimit.BMkappa <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))
  
  if(is.Global(o$X0)) {
    o$X0[1:6] <- 4
  }
  
  if(is.ScalarParameter(o$kappa)) {
    o$kappa[1] <- 2
  }
  o
}
