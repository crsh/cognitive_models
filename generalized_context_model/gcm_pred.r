# param  A vector of starting parameters: c(w1, c, b)
#             w = Attentional weight for dimension 1 of the psychological similarity space (assuming two dimensions)
#             c = Similarity sensitivity
#             b = Bias towards category 1
# mem   A matrix of exemplars in memory with one column for each dimension in psychological space
# obs   A matrix of observed exemplars with one column for each dimension in psychological space
# rho   An integer determining the distance metric in psychological space (1 = City block distance; 2 = Eucledian distance)

gcm_pred <- function(param, mem, obs, rho = 2) {
  w <- param[1]
  w[2] <- 1-w
  c <- param[2]
  b <- param[3]

  # Prepare objects
  m <- 1
  n_obs <- nrow(obs)
  mem <- as.matrix(mem)
  obs <- as.matrix(obs)
  all_resp <- matrix(rep(NA, n_obs), nrow = n_obs)

  # Model computations
  for(i in 1:n_obs) {
    iobs <- as.vector(obs[i, 1:ncol(obs)])

    ## Determine similarities & activation
    d <- w*abs(iobs - t(mem[, 1:2]))^rho
    d <- colSums(d)^(1/rho)           # Eq. 1, Nosofsky, Little, Donkin, & Fific (2011)
    s <- exp(-c*d^2)                  # Eq. 2, Nosofsky, Little, Donkin, & Fific (2011)
    a <- m*s
    s_ab <- b*sum(a[mem[, 3] == 1]) + (1-b)*sum(a[mem[, 3] == 2])

    ## Compute response probability for category 1
    p <- b*sum(a[mem[, 3] == 1])/s_ab # Eq. 13b, Nosofsky & Palmeri (1997)
    all_resp[i,] <- p
  }
  return(all_resp)
}
