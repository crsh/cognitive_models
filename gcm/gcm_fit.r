# data  A matrix or data.frame with frequencies of category 1 responses, category 2 responses and total number of responses as columns
# ...   Arguments to be passed to gcm_pred()

gcm_fit <- function(data, ...) {
  pred <- gcm_pred(...)
  data$pred <- pred
  dev <- -sum(apply(data, 1, function(x) dbinom(x[1], x[3], x[4], log = TRUE)))
  return(dev)
}
