# param  A vector of starting parameters: c(w1, c, b)
#             w = Attentional weight for dimension 1 of the psychological similarity space (assuming two dimensions)
#             c = Similarity sensitivity
#             b = Bias towards category 1
# mem   A matrix of exemplars in memory with one column for each dimension in psychological space
# obs   A matrix of observed exemplars with one column for each dimension in psychological space
# rho   An integer determining the distance metric in psychological space (1 = City block distance; 2 = Eucledian distance)
# p     An integer determining the form of the similarity function (1 = Exponential; 2 = Gaussian)

gcm_pred <- function(param, mem, obs, rho = 2, p = 2) {
  w <- param[1]
  w[2] <- 1-w
  c <- param[2]
  b <- param[3]

  # Prepare objects
  n_obs <- nrow(obs)
  mem <- as.matrix(mem)
  obs <- as.matrix(obs)
  all_resp <- matrix(rep(NA, n_obs), nrow = n_obs)

  # Model computations
  for(i in 1:n_obs) {
    iobs <- as.vector(obs[i, 1:ncol(obs)])

    ## Determine similarities & activation
    d <- w*abs(iobs - t(mem[, 1:2]))^rho
    d <- colSums(d)^(1/rho)           # Eq. 3, Nosofsky (1988)
    s <- exp(-c*d^p)                  # Eq. 4, Nosofsky (1989)
    s_ab <- b*sum(s[mem[, 3] == 1]) + (1-b)*sum(s[mem[, 3] == 2])

    ## Compute response probability for category 1
    p_a <- b*sum(s[mem[, 3] == 1])/s_ab # Eq. 2, Nosofsky (1989)
    all_resp[i,] <- p_a
  }
  return(all_resp)
}
