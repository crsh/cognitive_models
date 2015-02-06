probe_memory <- function(probe, memory, normalize = FALSE) {
  similarity <- colSums(probe * t(memory)) / colSums((probe != 0 | t(memory) != 0)) # Eq. 1, Hintzman (1984)
  activation <- similarity^3  # Eq. 2, Hintzman (1984)
  echo_intensity <- sum(activation)  # Eq. 3, Hintzman (1984)
  echo_content <- colSums(activation * memory) # Eq. 4, Hintzman (1984)
  if(normalize) echo_content <- echo_content / max(abs(echo_content))

  list(content = echo_content, intensity = echo_intensity)
}

encode <- function(episode, memory, p_encode) {
  encoding_error <- rbinom(length(episode), 1, p_encode)
  new_memory <- rbind(memory, episode * encoding_error)

  new_memory
}

forget <- function(memory, p_forget) {
  forgetting <- rbinom(length(memory), 1, p_forget)
  forgetting <- matrix(forgetting, ncol = ncol(memory))
  new_memory <- memory * forgetting

  new_memory
}
