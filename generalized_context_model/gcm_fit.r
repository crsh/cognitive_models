gcm_fit <- function(par, mem, obs, rho = 2, data) {
  pred <- gcm_pred(par, mem, obs, rho)
  data$pred <- pred
  dev <- -sum(apply(data, 1, function(x) dbinom(x[1], x[3], x[4], log = TRUE)))
  return(dev)
}
