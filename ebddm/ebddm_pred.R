ebddm_pred <- function(param, mem, obs, rho = 2, n_trials = 100) {
  require("RWiener")

  # Define parameters
  # ndim <- ncol(mem) - 1 # -0, if no category information
  w <- param["w"]
  w[2] <- 1 - sum(w)
  c <- param["c"]
  alpha <- param["alpha"]
  tau <- param["tau"]
  beta <- param["beta"]

  # Prepare objects
  n_obs <- nrow(obs)
  mem <- as.matrix(mem)
  obs <- as.matrix(obs)
  results <- data.frame()

  # Model computations
  for(i in 1:n_obs) {
    iobs <- as.vector(obs[i, 1:ncol(obs)])

    ## Determine similarities & activation
    d <- w*abs(iobs - t(mem[, 1:2]))^rho
    d <- colSums(d)^(1/rho)             # Eq. 3, Nosofsky (1988)
    s <- exp(-c*d)                      # Eq. 4, Nosofsky (1989)
    s_ab <- sum(s[mem[, 3] == 1]) + sum(s[mem[, 3] == 2])

    ## Compute response probability for category 1
    p <- sum(s[mem[, 3] == 1]) / s_ab   # Eq. 2, Nosofsky (1989)
    delta <- log(p / (1 - p))

    results <- rbind(results, cbind(stimulus = i, rwiener(n_trials, alpha, tau, beta, delta), delta = delta))
  }

  results
}
