# Simplified conjoint recognition model for aggregated responses

data {
  dimx <- dim(x)
  n_subject <- dimx[1]
  n_c1 <- dimx[2] # Number of conditions in factor 1
  n_c2 <- dimx[3] # Number of conditions in factor 2
}

model {
  
  # Data generating model ---------------------------------------------------

  for(i in 1:n_subject) {
    for(c1 in 1:n_c1) {
      for(c2 in 1:n_c2) {
        ## Targets
        x[i, c1, c2, 1] ~ dbin(V[i, c1, c2] + (1 - V[i, c1, c2]) * (G[i, c1, c2] + (1 - G[i, c1, c2]) * b[i]), n_items[1])

        ## Lure distractors
        x[i, c1, c2, 2] ~ dbin(G[i, c1, c2] + (1 - G[i, c1, c2]) * b[i], n_items[1])
      }
    }
    y[i] ~ dbin(b[i], n_items[2]) # New distractors
  }

  # Prior ------------------------------------------------------------------

  for(i in 1:n_subject) {
    for(c1 in 1:n_c1) {
      for(c2 in 1:n_c2) {
        V[i, c1, c2] ~ dbeta(1, 1)
        G[i, c1, c2] ~ dbeta(1, 1)
      }
    }
    b[i] ~ dbeta(1, 1)
  }
}
