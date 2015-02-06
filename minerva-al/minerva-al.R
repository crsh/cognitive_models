# probe           Feature vector for probe event (A).
# memory          Memory matrix with columns representin features and rows traces in memory (Mij).
# cue_features    A vector giving the indeces of features that are associated with cues.

probe_memory <- function (probe, memory, cue_features) {
  if(is.null(memory)) { # Empty memory
    echo <- runif(length(probe), -0.001, 0.001) # First trial is noise (p. 65, Jamieson, Crump, & Hannah, 2012)
    normalized_echo <- echo / max(echo) # Eq. 4, Jamieson, Crump, & Hannah (2012)

    return(normalized_echo)
  } else {
    # Compare only features associated with cues (p. 64, Jamieson, Crump, & Hannah, 2012)
    relevant_features <- rep(FALSE, length(probe))
    relevant_features[cue_features] <- TRUE
    probe <- probe[relevant_features]
    relevant_memory <- matrix(memory[, relevant_features], ncol = sum(relevant_features))

    # Calculate echo
#     similarity <- colSums(probe * t(relevant_memory)) / (sqrt(sum(probe^2)) * sqrt(rowSums(relevant_memory^2))) # Eq. 7, Jamieson, Crump, & Hannah (2012)
    similarity <- colSums(probe * t(relevant_memory)) / (sqrt(sum(probe^2) * rowSums(relevant_memory^2))) # simplified Eq. 7, Jamieson, Crump, & Hannah (2012)
    activation <- similarity^3                          # Eq. 2, Jamieson, Crump, & Hannah (2012)
    echo <- colSums(activation * memory)                # Eq. 3, Jamieson, Crump, & Hannah (2012)
    echo <- echo + runif(length(echo), -0.001, 0.001)   # Add noise (p. 64, Jamieson, Crump, & Hannah, 2012)
    normalized_echo <- echo / max(echo)                 # Eq. 4, Jamieson, Crump, & Hannah (2012)

    return(normalized_echo)
  }
}


# outcome             Feature vector for outcome event (X).
# normalized_echo     Normalized echo (C'j) produced by probe event (A).

expect_event <- function (outcome, normalized_echo) {
#     expectancy <- sum(outcome * normalized_echo) / sum(outcome != 0) # Eq. 5, Jamieson, Crump, & Hannah (2010)
  expectancy <- sum(outcome * normalized_echo) / sum(outcome != 0 & normalized_echo != 0)   # Eq. 4, Jamieson, Crump, & Hannah (2012)

  expectancy
}


# normalized_echo     Normalized echo (C'j) produced by probe event (A).
# event               Feature vector for the encountered event (e.g., E = A + X).
# p_encode            Probability with which a feature is encoded in memory (L).
# memory              Memory matrix with columns representin features and rows traces in memory (Mij).

learn <- function (normalized_echo, event, p_encode, memory) {
  # Probability a feature is encoded in memory
  if (p_encode < 1) { # Speeds up simulation
    encoding_error <- rbinom(length(event), 1, p_encode)
  } else {
    encoding_error <- rep(1, length(event))
  }

  # Discrepency encoding
  memory <- rbind(memory, (event - normalized_echo) * encoding_error)    # Eq. 6, Jamieson, Crump, & Hannah, (2012)

  memory
}
