# param  A vector of starting parameters: c(w, c, k)
#             w = Vector of two attentional weights for dimension 1 and 2 of the psychological similarity space (assuming three dimensions)
#             c = Similarity sensitivity
#             k = Response criterion parameter
# mem   A matrix of exemplars in memory with one column for each dimension in psychological space
# obs   A matrix of observed exemplars with one column for each dimension in psychological space
# rho   An integer determining the distance metric in psychological space (1 = City block distance; 2 = Eucledian distance)
# p     An integer determining the form of the similarity function (1 = Exponential; 2 = Gaussian)

gcm_rec_pred <- function(param, mem, obs, rho = 2, p = 1) {
  w <- param["w1"]
  w[2] <- param["w2"]
  w[3] <- 1-sum(w)
  c <- param["c"]
  k <- param["k"]

  if(is.na(w[1])) {             # Assume one psychological dimension
    w <- 1
  } else if(is.na(w[2])) {      # Assume two psychological dimensions
    w[2] <- 1 - w[1]
    w <- w[-3]
  }
  if(is.na(param["m"])) m <- 1  # Memory strength

  # Prepare objects
  n_obs <- nrow(obs)
  mem_cat <- mem[, "category"]
  obs_cat <- obs[, "category"]
  mem <- as.matrix(subset(mem, select = -category))
  ndim <- ncol(mem)
  obs <- as.matrix(subset(obs, select = -category))
  all_resp <- matrix(rep(NA, n_obs), nrow = n_obs)

  # Model computations
  for(i in 1:n_obs) {
    iobs <- as.vector(obs[i, ])

    ## Determine similarities & activation
    if(sum(mem_cat %in% obs_cat[i]) > 0) { # Only compare studied items to memory
      d <- w*abs(iobs - t(mem[mem_cat %in% obs_cat[i], ]))^rho        # Only compare to items from the same category, where similarity information is available
      d <- colSums(d)^(1/rho)             # Eq. 3, Nosofsky (1988)
      if(!is.na(drest)) d <- c(d, rep(drest, nrow(mem) - length(d)))  # Add assumed mean distance to each items from other categories
      s <- exp(-c*d^p)
      if(!is.na(srest)) s <- c(s, srest*(nrow(mem) - length(s)))      # Add assumed similarity for each items from other categories
    } else {
      if(!is.na(srest)) {
        s <- srest*nrow(mem)
      } else if(!is.na(drest)) {
        d <- rep(drest, nrow(mem))
        s <- exp(-c*d^p)
      } else {                            # Assume no inter-category similarity
        s <- 0
      }
    }

    a <- m*s
    f <- sum(a)                           # Eq. 8, Shin & Nosofsky (1992)

    if(pred == "single") {
      iresp <- ifelse(f > k, 1, 0)        # adapted from Eq. 6, Nosofksy (1991)
    } else if(pred == "prop") {
      iresp <- f/(f+k)                    # Eq. 9, Shin & Nosofsky (1992)
    }

    all_resp[i,] <- iresp
  }
  return(all_resp)
}
