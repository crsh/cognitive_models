
# Libraries ---------------------------------------------------------------

library("runjags")
library("dplyr")
library("ggplot2")
library("abind")

runjags.options(mode.continuous = TRUE)

# Prediction functions ----------------------------------------------------

cr_G <- function(fa, b) return((fa - b) / (1 - b))
cr_V <- function(hits, G, b) return((hits - G - (1 - G) * b) / (1 - G - (1 - G) * b))

cr_fa <- function(G, b) return(G + (1-G) * b)
cr_hits <- function(V, G, b) return(V + (1-V) * G + (1-V) * (1-G) * b)

cr_pred <- function(V, V_delta, G, G_delta, b, n_items) {
  n_subjects <- nrow(V)
  V2 <- V + V_delta
  V2[V2 > 1] <- 1
  V2[V2 < 0] <- 0
  G2 <- G + G_delta
  G2[G2 > 1] <- 1
  G2[G2 < 0] <- 0
  V <- abind::abind(V, V2, along = 3)
  G <- abind::abind(G, G2, along = 3)
  x <- array(NA, dim = c(n_subjects, ncol(V), 2, 2))
  y <- c()

  for(i in 1:n_subjects) {
    for(c1 in 1:ncol(V)) {
      for(c2 in 1:(length(V_delta) + 1)) {
        # Targets
        x[i, c1, c2, 1] <- V[i, c1, c2] + (1 - V[i, c1, c2]) * G[i, c1, c2] + (1 - V[i, c1, c2]) * (1 - G[i, c1, c2]) * b[i]

        # Lure distractors
        x[i, c1, c2, 2] <- G[i, c1, c2] + (1 - G[i, c1, c2]) * b[i]
      }
    }
    y[i] <- b[i] # New distractors
  }

  x <- round(x * n_items[1])
  y <- round(y * n_items[2])

  return(list(x = x, y = y, param = list(V = V, G = G, b = b, n_items = n_items, n_subjects = n_subjects)))
}


# Init functions ----------------------------------------------------------

cr_init <- function(V, G, b, chains) {

  inits <- list()
  for(i in 1:chains) {
    V_init <- if(!is.null(V)) array(runif(length(V), 0, 1), dim = dim(V)) else NULL
    G_init <- if(!is.null(G)) array(runif(length(G), 0, 1), dim = dim(G)) else NULL
    b_init <- if(!is.null(b)) runif(length(b), 0, 1) else NULL

    i_inits <- list(V = V_init, G = G_init, b = b_init)
    inits[[i]] <- Filter(Negate(function(x) is.null(unlist(x))), i_inits)
  }

  inits
}


# Plotting functions ------------------------------------------------------

plot_deviations <- function(samples, truth) {
  recovered <- as.data.frame(summary(samples))
  recovered$param <- factor(gsub("[^GVb]", "", rownames(recovered)))
  recovered$lie <- recovered$Mode - truth

  recovered %>%
    group_by(param) %>%
    ggplot(aes(x = lie, color = param)) +
    geom_density(aes(fill = param), alpha = 0.5, adjust = 0.5) +
    geom_vline(xintercept = 0, linetype = 2) +
    theme_minimal() +
    xlab("Estimation error")
}


# One condition, aggregated responses -------------------------------------

n_subject <- 25
V <- 0.5
G <- 0.5
b <- 0.5
n_items <- 200
n_items <- c(n_items, 1.5 * n_items)

synthetic <- list(x = matrix(NA, ncol = 2, nrow = n_subject), y = NA)
synthetic$x[, 1] <- rep(cr_hits(V, G, b) * n_items[1], n_subject)
synthetic$x[, 2] <- rep(cr_fa(G, b) * n_items[1], n_subject)
synthetic$y <- rep(b * n_items[2], n_subject)
synthetic$n_items <- n_items
synthetic$n_subject <- n_subject

poi <- c("V", "G", "b")

inits <- list(V = runif(n_subject, 0, 1), G = runif(n_subject, 0, 1), b = runif(n_subject, 0, 1))

cr_samples <- run.jags(
  model = "cr_model.txt"
  , monitor = poi
  , inits = list(inits, inits, inits)
  , data = synthetic
  , n.chains = 3
  , sample = 1e4
  , burnin = 1e3
  , thin = 10
  , method = "rjparallel"
)

plot(cr_samples)

plot_deviations(cr_samples, c(V, G, b))
rm(cr_samples)

# Aggregated responses ----------------------------------------------------

V <- matrix(c(0.5, 0.5), ncol = 2)
V_delta <- 0
G <- matrix(c(0.5, 0.5), ncol = 2)
G_delta <- 0
b <- 0.5
n_items <- 324
n_items <- c(n_items, 1.5 * n_items)

synthetic <- cr_pred(V, V_delta, G, G_delta, b, n_items)

synthetic_jags <- synthetic[c("x", "y")]
synthetic_jags$n_subject <- synthetic$param$n_subjects
synthetic_jags$n_items <- synthetic$param$n_items

poi <- c("V", "G", "b")

cr_samples <- run.jags(
  model = "cr_two_within_conditions.txt"
  , monitor = poi
  , inits = cr_init(synthetic$param[[c("V")]], synthetic$param[[c("G")]], synthetic$param[[c("b")]], 3)
  , data = synthetic_jags
  , n.chains = 3
  , sample = 5e4
  , burnin = 5e4
  , thin = 10
  , method = "rjparallel"
)

plot(cr_samples)

plot_deviations(cr_samples, c(synthetic$param$V, synthetic$param$G, synthetic$param$b))
rm(cr_samples)


# Individual participants -------------------------------------------------

n_subjects <- 25
V <- matrix(rep(c(0.5, 0.5), each = n_subjects), ncol = 2)
V_delta <- 0
G <- matrix(rep(c(0.5, 0.5), each = n_subjects), ncol = 2)
G_delta <- 0
b <- rep(0.5, n_subjects)
n_items <- 324
n_items <- c(n_items, 1.5 * n_items)

synthetic <- cr_pred(V, V_delta, G, G_delta, b, n_items)

synthetic_jags <- synthetic[c("x", "y")]
synthetic_jags$n_subject <- synthetic$param$n_subjects
synthetic_jags$n_items <- synthetic$param$n_items

poi <- c("V", "G", "b")

cr_samples <- run.jags(
  model = "cr_two_within_conditions.txt"
  , monitor = poi
  , data = synthetic_jags
  , inits = cr_init(synthetic$param[[c("V")]], synthetic$param[[c("G")]], synthetic$param[[c("b")]], 3)
  , n.chains = 3
  , sample = 5e4
  , burnin = 5e4
  , thin = 10
  , method = "rjparallel"
)

plot_deviations(cr_samples, c(synthetic$param$V, synthetic$param$G, synthetic$param$b))
rm(cr_samples)

# Hierarchical model with correlated participant effects ------------------

n_subjects <- 10
V <- matrix(pnorm(rnorm(n_subjects * 2, rep(qnorm(c(0.5, 0.5)), each = n_subjects), 0.75)), ncol = 2)
G <- matrix(pnorm(rnorm(n_subjects * 2, rep(qnorm(c(0.5, 0.5)), each = n_subjects), 0.75)), ncol = 2)
b <- pnorm(rnorm(n_subjects, qnorm(0.05), 0.75))
V_delta <- 0
G_delta <- 0
n_items <- 324
n_items <- c(n_items, 1.5 * n_items)

synthetic <- cr_pred(V, V_delta, G, G_delta, b, n_items)

synthetic_jags <- synthetic[c("x", "y")]
synthetic_jags$n_subject <- synthetic$param$n_subjects
synthetic_jags$n_items <- synthetic$param$n_items
synthetic_jags$I_part <- diag((ncol(V) + ncol(G)) * 2 + 1)  # Identiy matrix for participants

poi <- c("V", "G", "b", "xi_part", "sigma", "pd", "dic")

cr_samples <- run.jags(
  model = "hierarchical_cr_two_within_conditions0.R"
  , monitor = poi
  , data = synthetic_jags
  , n.chains = 3
  , sample = 1e4
  , burnin = 5e3
  , thin = 10
  , method = "rjparallel"
)

str(cr_samples$end.state)

param_varcor <- function(x) {
  matches <- runjags:::matchvars(runjags:::checkvalidmonitorname("sigma_inv"), varnames(x))
  sigma_cols <- varnames(x)[matches]
  n_cols <- sqrt(length(sigma_cols))
  sigma_inv <- x[, sigma_cols, drop = FALSE]
  sigma <- apply(
    sigma_inv
    , 1
    , function(y) {
      sigma_inv_matrix <- matrix(y, ncol = n_cols)
      sigma_matrix <- solve(sigma_inv_matrix)
      as.vector(sigma_matrix)
    }
  )
  sigma <- t(sigma)
  colnames(sigma) <- gsub("_inv", "", colnames(sigma_inv))

  sigma_V <- x[, "xi_part[1]"] * sqrt(sigma[, "sigma[1,1]"])
  sigma_G <- x[, "xi_part[2]"] * sqrt(sigma[, "sigma[2,2]"])
  sigma_b <- x[, "xi_part[3]"] * sqrt(sigma[, "sigma[3,3]"])

  rho <- matrix(NA, nrow = nrow(x), ncol = n_cols ^ 2)
  colnames(rho) <- gsub("sigma", "rho", colnames(sigma))
  for (i in 1:n_cols) {
    for (j in 1:n_cols) {
      rho[, paste0("rho[", i, ",", j, "]")] <-
        sigma[, paste0("sigma[", i, ",", j, "]")] / sqrt(sigma[, paste0("sigma[", i, ",", i, "]")] * sigma[, paste0("sigma[", j, ",", j, "]")])
    }
  }

  cbind(cbind(sigma_V, sigma_G, sigma_b), rho)
}

add.summary(cr_samples, mutate = param_cor, vars = c("V", "G", "b", "sigma_V", "sigma_G", "sigma_b", "rho"))

plot_deviations(cr_samples, unlist(list(c(V, V + V_delta), c(G, G + G_delta), b)))
