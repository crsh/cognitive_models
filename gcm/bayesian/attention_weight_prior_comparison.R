library("runjags")
library("rstan")

trials <- 3 * 50

tests <- c()
data <- c()
for(i in 1:3) {
  tests <- rbind(tests, read.csv2(paste0("../data/shin_nosofsky_1992_cat", i, ".csv")))
  data <- rbind(data, read.csv2(paste0("../data/shin_nosofsky_1992_responses_cat", i, ".csv")))
}
data$response <- round(data$Observed * trials)

recognition_data <- list(
  y = data$response
  , tests = as.matrix(tests[, -c(1, 5:7)])
  , memory = as.matrix(subset(tests, Exemplar %in% paste0("O", 1:6))[, -c(1, 5:7)])
  , ntrials = rep(trials, nrow(tests))
  , p = 1 # Shape of relationship between similarity and psychological distance
  , rho = 2 # Power of the Minkowski distance
)

poi <- c("c", "w", "k", "pred_y")


# Fit Dirichlet prior -----------------------------------------------------

model1 <- system.time(
  gcm_samples1 <- run.jags(
    model = "GCM_agg_recognition.txt"
    , monitor = poi
    , data = recognition_data
    , n.chains = 3
    , sample = 5e4
    , burnin = 100
    , thin = 1
    , method = "rjparallel"
  )
)


# Fit Gelman prior --------------------------------------------------------

model2 <- system.time(
  gcm_samples2 <- run.jags(
    model = "GCM_agg_recognition2.txt"
    , monitor = poi
    , data = recognition_data
    , n.chains = 3
    , sample = 5e4
    , burnin = 100
    , thin = 1
    , method = "rjparallel"
  )
)


model1 - model2



# STAN --------------------------------------------------------------------

memory <- as.matrix(subset(tests, Exemplar %in% paste0("O", 1:6))[, -c(1, 5:7)])

recognition_data <- list(
  y = data$response
  , tests = as.matrix(tests[, -c(1, 5:7)])
  , memory = memory
  , ntrials = rep(trials, nrow(tests))
  , p = 1 # Shape of relationship between similarity and psychological distance
  , rho = 2 # Power of the Minkowski distance
  , ntests = dim(tests)[1]
  , nmemory = dim(memory)[1]
  , ndim = dim(memory)[2]
)

# Fit Dirichlet prior -----------------------------------------------------

model3 <- system.time(
  gcm_samples3 <- stan(
    file = "GCM_agg_recognition.stan"
    , pars = poi
    , data = recognition_data
    , chains = 3
    , iter = 5e4 + 1100
    , warmup = 1100
    , thin = 1
    , cores = 3
  )
)

model1 - model3

# Fit Gelman prior --------------------------------------------------------

model4 <- system.time(
  gcm_samples4 <- stan(
    file = "GCM_agg_recognition2.stan"
    , pars = poi
    , data = recognition_data
    , chains = 3
    , iter = 5e4 + 1100
    , warmup = 1100
    , thin = 1
    , cores = 3
  )
)

model2 - model4
