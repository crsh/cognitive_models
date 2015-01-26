gcm_rec_fit <- function(
  par
  , data
  , mem
  , n
  , minimize = "individual"
  , design = NULL
  , ...
) {
  pred <- gcm_rec_pred(param = par, mem = mem, ...)
  pred <- ifelse(pred == 1, 0.99999, pred)
  pred <- ifelse(pred == 0, 0.00001, pred)

  if(minimize == "individual") {
    dev <- -2*sum(dbinom(data[, "response"], n, pred, log = TRUE))
  } else if(minimize == "condition") {
    cond_data <- aggregate(as.vector(data[, "response"]), by = design, FUN = sum)
    cond_n <- aggregate(as.vector(data[, "response"]), by = design, FUN = function(x) length(x)*n)
    cond_pred <- aggregate(as.vector(data[, "response"]), by = design, FUN = mean)

    dev <- -2*sum(dbinom(cond_data$x, cond_n$x, cond_pred$x, log = TRUE))
  }
  return(dev)
}
