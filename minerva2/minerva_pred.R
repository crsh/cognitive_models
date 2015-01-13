# Co-occurance of two events
a <- c(0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0)
b <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0)
e <- a + b

# Storing in memory
n_features <- 12
memory <- matrix(ncol = n_features)
encode_error <- function(p, n) rbinom(n, 1, p)

p_store <- runif(1) # Probability a feature is encoded in memory
memory[1, ] <- e * encode_error(p_store, n_features)

# Probing memory
memory <- matrix( # Create memory with 5 traces á 12 features
  sample(c(-1, 0, 1), size = n_features*5, replace = TRUE)
  , ncol = n_features
  , dimnames = list(
    rows = paste0("trace", 1:5)
    , cols = paste0("feature", 1:n_features)
  )
)
probe <- sample(c(1, -1, 0), size = n_features, replace = TRUE) # Create probe
memory <- rbind(memory, probe) # Add probe to memory

## Echo
similarity <- colSums(probe * t(memory)) / colSums((probe != 0 | t(memory) != 0)) # Eq. 1, Jamieson, Crump, & Hannah (2012)
activation <- similarity^3 # Eq. 2, Jamieson, Crump, & Hannah (2012)
echo <- colSums(activation * memory) # Eq. 3, Jamieson, Crump, & Hannah (2012)
echo_intensity <- sum(activation) # Eq. 3, Hintzman (1984)

## Cued recall
normalized_echo <- echo / max(abs(echo)) # Eq. 4, Jamieson, Crump, & Hannah (2012)
retrieval <- sum(probe * normalized_echo) / sum(event != 0 | normalized_echo != 0) # Eq. 5, Jamieson, Crump, & Hannah (2012)






probe <- c(0, 1, -1, 1)
memory <- t(as.matrix(data.frame(c(0, 1, 0, 1), c(1, 1, -1, 0), c(-1, 0, 0, -1))))

## Echo
similarity <- colSums(probe * t(memory)) / colSums((probe != 0 | t(memory) != 0)) # Eq. 1, Jamieson, Crump, & Hannah (2012)
activation <- similarity^3 # Eq. 2, Jamieson, Crump, & Hannah (2012)
echo <- colSums(activation * memory) # Eq. 3, Jamieson, Crump, & Hannah (2012)
echo_intensity <- sum(activation) # Eq. 3, Hintzman (1984)

## Cued recall
normalized_echo <- echo / max(abs(echo)) # Eq. 4, Jamieson, Crump, & Hannah (2012)
retrieval <- sum(probe * normalized_echo) / sum(((probe != 0) + (normalized_echo != 0)) >= 1) # Eq. 5, Jamieson, Crump, & Hannah (2012)
