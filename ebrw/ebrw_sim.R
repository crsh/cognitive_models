ebrw_sim <- function(param, mem, m, obs, rho = 2, n_trial = 1000) {
  
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
  mem_sim <- as.matrix(mem[, 1:ndim])
  obs_sim <- as.matrix(obs[, 1:ndim])
  results <- expand.grid(stimulus = 1:n_obs, trial = 1:n_trial, response = NA, rt = NA)
  results <- results[order(results$stimulus), ]
  
  # Model computations
  for(i in 1:n_obs) {
    iobs <- as.vector(obs_sim[i, 1:ndim])
    
    for(j in 1:n_trial) {
      ## Determine similarities & activation
      d <- w*abs(iobs - t(mem[, 1:2]))^rho
      d <- colSums(d)^(1/rho)           # Eq. 3, Nosofsky (1988)
      s <- exp(-c*d)                    # Eq. 4, Nosofsky (1989)
      s_ab <- sum(s[mem[, 3] == 1]) + sum(s[mem[, 3] == 2])
      
      ## Compute response probability for category 1
      p <- sum(s[mem[, 3] == 1])/s_ab   # Eq. 2, Nosofsky (1989)
      q <- 1-p
      t_step <- alpha + 1/s_ab          # Eq. 10, Nosofsky & Palmeri (1997)
      
      ## Simulate responses
      rw_count <- 0
      n_steps <- 0
      old <- 0
      
      while(rw_count < A && rw_count > -B) {
        rw_count <- rw_count + sample(c(1, -1), 1, prob = c(p, q))
        n_steps <- n_steps + 1
      }
    
      if(rw_count >= A) {
        results[results$stimulus == i & results$trial == j, c("response", "rt")] <- c(1, (n_steps * t_step) * k + mu)
      } else {
        results[results$stimulus == i & results$trial == j, c("response", "rt")] <- c(0, (n_steps * t_step) * k + mu)
      }
    }
  }

  merge(results, cbind(obs, stimulus = 1:nrow(obs)))
}