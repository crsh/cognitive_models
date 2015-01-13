gcm_rec_fit <- function(
  par
  , data
  , mem
  , dimensions = 1
  , rho = 2
  , n
  , minimize = "individual"
  , design = NULL
  , ...
) {
  switch(dimensions
    , dims <- "x1"
    , dims <- c("x2", "y2")
    , dims <- c("x3", "y3", "z3")
  )
  param <- c(par, ...) # Add not-to-be optimized parameters
  fix_par <- param[!grepl("(_\\D+|_\\d+)", names(param))] # Parameters fixed across all conditions
  vary_par <- param[grepl("(_\\D+|_\\d+)", names(param))] # Parameters varying across some conditions
  vary_par_names <- unique(gsub("_.+", "", names(vary_par))) # Get parameter types (e.g., c, k)

  dev <- c()
  n_obs <- c()
  if(!is.null(design) & length(fix_par) != length(param)) { # Fit with parameters varying across conditions
    for(i in 1:nrow(design)) {
      factor_levels <- as.character(design[i, ])
      ddata <- subset(data, eval(parse(text = paste(names(design[i, ]), factor_levels, sep = " == ", collapse = " & "))))

      par_regex <- paste0(vary_par_names, "_", factor_levels[1], "_", factor_levels[2], "|", vary_par_names, "_", factor_levels[1], "_x|", vary_par_names, "_x_", factor_levels[2], collapse = "|")
      dpar <- c(fix_par, vary_par[grep(par_regex, names(vary_par))])
      names(dpar) <- gsub("_.+", "", names(dpar))

      pred <- gcm_rec_pred(dpar, mem[, c("category", dims)], ddata[, c("category", dims)], rho)
      pred <- ifelse(pred == 1, 0.99999, pred)
      pred <- ifelse(pred == 0, 0.00001, pred)

      if(minimize == "individual") {
        ddev <- -2*sum(dbinom(ddata[, "response"], n, pred, log = TRUE))
      } else if(minimize == "condition") {
        cond_data <- aggregate(as.vector(ddata[, "response"]), by = design, FUN = sum)
        cond_n <- aggregate(as.vector(ddata[, "response"]), by = design, FUN = function(x) length(x)*n)
        cond_pred <- aggregate(as.vector(ddata[, "response"]), by = design, FUN = mean)

        ddev <- -2*sum(dbinom(cond_data$x, cond_n$x, cond_pred$x, log = TRUE))
      }

      dev <- sum(dev, ddev)
    }

    # Fit new probe condition

  } else { # Fit parameters fixed across all conditions
    pred <- gcm_rec_pred(param, mem[, c("category", dims)], data[, c("category", dims)], rho)
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
  }
  return(dev)
}
