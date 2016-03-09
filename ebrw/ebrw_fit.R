ebrw_fit <- function(data, ...) {
  predictions <- ebrw_pred(...)
  dev <- sum(((predictions[, 1] - data[, 1])^2)/(sd(data[, 1])/sqrt(length(data[, 1]))) + ((predictions[, 4] - data[, 2])^2)/(sd(data[, 2])/length(data[, 2]))) # Nosofsky & Stanton, 2005
  return(dev)
}
