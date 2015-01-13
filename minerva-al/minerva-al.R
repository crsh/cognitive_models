calc_echo <- function (probe, memory, cue_features) {
  if(is.null(memory)) { # Empty memory
    echo <- runif(length(probe), -0.001, 0.001) # First trial is noise (p. 65, Jamieson, Crump, & Hannah, 2012)
    normalized_echo <- echo / max(echo) # Eq. 4, Jamieson, Crump, & Hannah (2012)
    return(normalized_echo)
  } else {
    # Compare only features coded in probe (p. 64, Jamieson, Crump, & Hannah, 2012)
    relevant_features <- rep(FALSE, length(probe))
    relevant_features[cue_features] <- TRUE
    probe <- probe[relevant_features]
    relevant_memory <- matrix(memory[, relevant_features], ncol = sum(relevant_features))
    
    # Calculate echo
#     similarity <- colSums(probe * t(relevant_memory)) / (sqrt(sum(probe^2)) * sqrt(rowSums(relevant_memory^2))) # Eq. 7, Jamieson, Crump, & Hannah (2012)
    similarity <- colSums(probe * t(relevant_memory)) / (sqrt(sum(probe^2) * rowSums(relevant_memory^2))) # simplified Eq. 7, Jamieson, Crump, & Hannah (2012)
    activation <- similarity^3 # Eq. 2, Jamieson, Crump, & Hannah (2012)
    echo <- colSums(activation * memory) # Eq. 3, Jamieson, Crump, & Hannah (2012)
    echo <- echo + runif(length(echo), -0.001, 0.001) # Add noise (p. 64, Jamieson, Crump, & Hannah, 2012)
    normalized_echo <- echo / max(echo) # Eq. 4, Jamieson, Crump, & Hannah (2012)
    return(normalized_echo)
  }
}

calc_expectancy <- function (event, echo) {
#     expectancy <- sum(event * echo) / sum(event != 0) # Eq. 5, Jamieson, Crump, & Hannah (2010)
  expectancy <- sum(event * echo) / sum(event != 0 & echo != 0) # Eq. 4, Jamieson, Crump, & Hannah (2012)
  expectancy
}

learn <- function (echo, event, p_encode, memory) {
  # Probability a feature is encoded in memory
  if (p_encode < 1) { # Speeds up simulation
    encoding_error <- rbinom(length(event), 1, p_encode)
  } else {
    encoding_error <- rep(1, length(event))
  }
  
  # Discrepency encoding
  memory <- rbind(memory, (event - echo) * encoding_error) # Eq. 6, Jamieson, Crump, & Hannah, 2012
  memory
}
