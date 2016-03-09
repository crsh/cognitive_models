# param  A vector of starting parameters: c(w, c, A, alpha, k, mu)
#             w = A vector of attention weights
#             c = Similarity sensitivity
#             K = Threshold for 'Category A' response (equal thresholds assumed)
#             alpha = Retrieval constant
#             k = Response time scaling constant
#             mu = Reponse time constant
# mem   A matrix of exemplars in memory with one column for each dimension in psychological space
# obs   A matrix of observed exemplars with one column for each dimension in psychological space
# rho   An integer determining the distance metric in psychological space (2 = Eucledian distance; 1 = City block distance)

ebrw_pred <- function(param, mem, obs, rho = 2) {

  # Define parameters
  ndim <- ncol(mem) - 1 # -0, if no category information
  w <- param[1:(ndim - 1)]
  w[ndim] <- 1 - sum(w)
  c <- param[ndim]
  A <- param[ndim + 1]
  B <- A
  alpha <- param[ndim + 2]
  k <- param[ndim + 3]
  mu <- param[ndim + 4]


  # Prepare objects
  n_obs <- nrow(obs)
  mem <- as.matrix(mem)
  obs <- as.matrix(obs)
  mean_rt <- rep(NA, n_obs)
  rt_A <- rep(NA, n_obs)
  rt_B <- rep(NA, n_obs)
  accuracy <- rep(NA, n_obs)

  # Model computations
  for(i in 1:n_obs) {
    iobs <- as.vector(obs[i, 1:(ncol(obs)-1)])

    ## Determine similarities & activation
    d <- w*abs(iobs - t(mem[, 1:2]))^rho
    d <- colSums(d)^(1/rho)           # Eq. 3, Nosofsky (1988)
    s <- exp(-c*d)                  # Eq. 4, Nosofsky (1989)
    s_ab <- sum(s[mem[, 3] == 1]) + sum(s[mem[, 3] == 2])

    ## Compute response probability for category 1
    p <- sum(s[mem[, 3] == 1])/s_ab   # Eq. 2, Nosofsky (1989)
    q <- 1-p
    t_step <- alpha + 1/s_ab          # Eq. 10, Nosofsky & Palmeri (1997)

    if(p != 0.5) {
      p_A <- (1-(q/p)^B) / (1-(q/p)^(A+B))              # Eq. 16a, Nosofsky & Palmeri (1997)

      theta1 <- ((p/q)^(A+B) + 1) / ((p/q)^(A+B) - 1)     # Eq. 19, Nosofsky & Palmeri (1997)
      theta2 <- ((p/q)^B + 1) / ((p/q)^B - 1)             # Eq. 19, Nosofsky & Palmeri (1997)
      n_step_A <- 1/(p-q) * (theta1*(A+B) - theta2*B)   # Eq. 18a, Nosofsky & Palmeri (1997)

      theta1 <- ((p/q)^-(A+B) + 1) / ((p/q)^-(A+B) - 1)   # Eq. 21, Nosofsky & Palmeri (1997)
      theta2 <- ((p/q)^-A + 1) / ((p/q)^-A - 1)           # Eq. 21, Nosofsky & Palmeri (1997)
      n_step_B <- 1/(q-p) * (theta1*(A+B) - theta2*A)   # Eq. 20a, Nosofsky & Palmeri (1997)

      n_steps <- B/(q-p) - (A+B)/(q-p) * ((1-(q/p)^B)/(1-(q/p)^(A+B))) # Eq. 14a, Nosofsky & Palmeri (1997)

    } else {
      p_A <- B/(A+B)                                    # Eq. 16b, Nosofsky & Palmeri (1997)

      n_step_A <- A/3*(2*B + A)                         # Eq. 18b, Nosofsky & Palmeri (1997)
      n_step_B <- B/3*(2*A + B)                         # Eq. 20b, Nosofsky & Palmeri (1997)

      n_steps <- A*B                                      # Eq. 14b, Nosofsky & Palmeri (1997)
    }

    rt_A[i] <- (n_step_A * t_step) * k + mu
    rt_B[i] <- (n_step_B * t_step) * k + mu
    mean_rt[i] <- (n_steps * t_step) * k + mu
    accuracy[i] <- if(obs[i, 3] == 1) p_A else 1 - p_A
  }

  pred <- cbind(accuracy, rt_A, rt_B, mean_rt)
  return(pred)
}
