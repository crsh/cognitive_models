# param  A vector of starting parameters: c(w, c, k)
#             w = Vector of two attentional weights for dimension 1 and 2 of the psychological similarity space (assuming three dimensions)
#             c = Similarity sensitivity
#             k = Response criterion parameter
# mem   A matrix of exemplars in memory with one column for each dimension in psychological space
# obs   A matrix of observed exemplars with one column for each dimension in psychological space
# rho   An integer determining the distance metric in psychological space (1 = City block distance; 2 = Eucledian distance)
# p     An integer determining the form of the similarity function (1 = Exponential; 2 = Gaussian)

gcm_rec_pred <- function(param, mem, obs, rho = 2, p = 1, pred = "prop") {
  w <- param["w1"]
  w[2] <- param["w2"]
  w[3] <- param["w3"]
  w[4] <- param["w4"]
  w[5] <- param["w5"]
  w[6] <- 1-sum(w)
  c <- param["c"]
  k <- param["k"]

  # Prepare objects
  n_obs <- nrow(obs)
  mem <- as.matrix(mem)
  ndim <- ncol(mem)
  obs <- as.matrix(obs)
  all_resp <- c()

  # Model computations
  for(i in 1:n_obs) {
    iobs <- as.vector(obs[i, ])

    ## Determine similarities & activation
    d <- w*abs(iobs - t(mem))^rho
    d <- colSums(d)^(1/rho)               # Eq. 3, Nosofsky (1988)
    s <- exp(-c*d^p)

    f <- sum(s)                           # Eq. 8, Shin & Nosofsky (1992)

    if(pred == "single") {
      iresp <- ifelse(f > k, 1, 0)        # adapted from Eq. 6, Nosofksy (1991)
    } else if(pred == "prop") {
      iresp <- f/(f+k)                    # Eq. 9, Shin & Nosofsky (1992)
    }

    all_resp[i] <- iresp
  }
  return(all_resp)
}
