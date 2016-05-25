# Simplified conjoint recognition model for aggregated responses

model {
  
  # Data generating model ---------------------------------------------------

  for(i in 1:n_subject) {

    ## Targets
    # x[i, 1] ~ dbin(V[i] + (1 - V[i]) * G[i] + (1 - V[i]) * (1 - G[i]) * b[i], n_items[1])
    x[i, 1] ~ dbin(V[i] + (1 - V[i]) * (G[i] + (1 - G[i]) * b[i]), n_items[1]) # This form speeds up computations

    ## Lure distractors
    x[i, 2] ~ dbin(G[i] + (1 - G[i]) * b[i], n_items[1])

    ## New distractors
    y[i] ~ dbin(b[i], n_items[2])
  }

  # Prior ------------------------------------------------------------------

  for(i in 1:n_subject) {
    V[i] ~ dbeta(1, 1)
    G[i] ~ dbeta(1, 1)
    b[i] ~ dbeta(1, 1)
  }
}
