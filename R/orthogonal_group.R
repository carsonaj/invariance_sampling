#' Calculate value of truncated normal density on nonnegative reals up to
#' multiplicative constant
#'
#' @param n The dimension of original space
#' @param s The element of Rn / O(n) at which the value of the truncated normal
#' density is calculated up to a multiplicative constant (the proposed state)
#' @param t The element of Rn / O(n) at which the mode of the truncated normal
#' density occurs (the current state)
#' @param sig The scale parameter of the truncated normal
#' @return value of trunctated normal density at t up to multiplicative constant
#' @examples
#' sample <- g_bar(s, t, sig);
#' @export
g_bar <- function(n, s, t, sig) {
  if (s > 0) {
    y <- t^(n-1) * exp(-1/(2 * sig^2) * (s - t)^2)
    return(y)
  } else {
    return(0)
  }
}












#' Sample value of truncated normal density on nonnegative reals
#'
#' @param s The element of Rn / O(n) at which the value of the truncated normal
#' density is calculated up to a multiplicative constant (the proposed state)
#' @param t The element of Rn / O(n) at which the mode of the truncated normal
#' density occurs (the current state)
#' @param sig The scale parameter of the truncated normal
#' @return value of trunctated normal density at t up to multiplicative constant
#' @examples
#' sample <- sample_g_bar(s, t, sig);
#' @export
sample_g_bar <- function(t, sig) {
  alpha <- -t / sig
  u <- runif(1)
  x <- qnorm(pnorm(alpha) + u * (1 - pnorm(alpha))) * sig + t

  return(x)
}















#' Sample from induced target distribution using Quotient Metropolis-Hastings
#' on Rn / O(n)
#' @param f_bar (function) The induced density on the quotient
#' @param t_0 (positive number) The starting location for the chain
#' @param n_samples (positive integer) The number of consecutive samples desired
#' @param sig (positive number) The standard dev of the proposal
#' @return Samples from target distribution
#' @examples
#' samples <- qmh_orthogonal_grp(f_bar, t_0, n_samp);
qmh_orthogonal_grp <- function(f_bar, t_0, n_samp, sig) {
  samples <- rep(0, n_samp)
  t <- t_0
  n <- 1
  while (n <= n_samp) {
    s <- sample_g_bar(t, sig)
    ratio <- (f_bar(s) * g_bar(n, t, s, sig)) / (f_bar(t) * g_bar(n, s, t, sig))
    print(ratio)
    print(" and ")
    print(s)
    alpha <- min(1, ratio)

    u <- runif(1)
    if (u <= alpha) {
      t <- s
    }
    samples[n] <- t
    n <- n + 1
  }

  return(samples)
}








#' Calculate value of truncated normal density on nonnegative reals up to
#' multiplicative constant
#' @param t1 (positive number) The mean of the first normal
#' distribution
#' @param t2 (positive number) The mean of the second normal
#' distribution
#' @param sig1 (positive number) The standard deviation of the first normal
#' distribution
#' @param sig2 (positive number) The standard deviation of the second normal
#' distribution
#' @return (function) The function that maps t to its value of the truncated normal mixture density
#' @examples
#' s <- f_bar_conc_bimodal(t);
f_bar_conc_bimodal <- function(t1, t2, sig1, sig2) {
  f_bar_conc_bimodal_specific <- function(s) {
    val <- dnorm(s, mean = t1, sd = sig1) + dnorm(s, mean = t2, sd = sig2)

    return(val)
  }

  return(f_bar_conc_bimodal_specific)
}











#' Sample uniformly on surphace of n-sphere
#'
#' @param n (positive integer) The dimension of the sphere
#' @param r (positive number) The radius of the sphere
#' @return (col matrix) A uniform sample on the radius of the sphere
#' @examples
#' sample <- runif_sphere(n, r);
#' @export
runif_sphere <- function(n, r) {
  samp_cube <- rnorm(n)
  samp_cube <- matrix(samp_cube, n, 1)
  samp_sphere <- r * samp_cube / norm(samp_cube, "F")

  return(samp_sphere)
}



