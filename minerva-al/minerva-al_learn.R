calc_echo <- function (memory, event) {
  similarity <- colSums(event * t(memory)) / (sqrt(sum(event^2)) * sqrt(rowSums(memory^2))) # Eq. 7, Jamieson, Crump, & Hannah (2012)
  activation <- similarity^3 # Eq. 2, Jamieson, Crump, & Hannah (2012)
  echo <- colSums(activation * memory) # Eq. 3, Jamieson, Crump, & Hannah (2012)
  echo <- echo + runif(length(event), -0.001, 0.001) # Add noise (p. 64, Jamieson, Crump, & Hannah, 2012)
  normalized_echo <- echo / max(echo) # Eq. 4, Jamieson, Crump, & Hannah (2012)
#   normalized_echo <- normalized_echo + runif(length(normalized_echo), -0.001, 0.001) # Add noise (p. 64, Jamieson, Crump, & Hannah, 2012)
  normalized_echo
}

calc_expectancy <- function (event, probe, memory) {
  if(is.null(memory)) { # Empty memory
    normalized_echo <- runif(length(probe), -0.001, 0.001)
    normalized_echo <- normalized_echo / max(normalized_echo)
  } else {
    normalized_echo <- calc_echo(memory, probe)
  }
  expectancy <- sum(event * normalized_echo) / (sqrt(sum(event^2)) * sqrt(sum(normalized_echo^2))) # Eq. 4, Jamieson, Crump, & Hannah (2012)
#   expectancy <- sum(event * normalized_echo) / sum(event != 0 | normalized_echo != 0) # Eq. 4, Jamieson, Crump, & Hannah (2012)
  expectancy
}

learn <- function (memory, event, p_encode) {
  if (p_encode < 1) {
    encoding_error <- rbinom(length(event), 1, p_encode) # Probability a feature is encoded in memory
  } else {
    encoding_error <- rep(1, length(event))
  }
  
  # Discrepency encoding
  if(is.null(memory)) { # Empty memory
    normalized_echo <- runif(length(probe), -0.001, 0.001)
    normalized_echo <- normalized_echo / max(normalized_echo)
    memory <- matrix((event - normalized_echo) * encoding_error, ncol = length(event))
    rownames(memory) <- deparse(substitute(event))
  } else {
    normalized_echo <- calc_echo(memory, event)
    memory <- rbind(memory, (event - normalized_echo) * encoding_error)
    rownames(memory)[nrow(memory)] <- deparse(substitute(event))
  }

  memory
}



