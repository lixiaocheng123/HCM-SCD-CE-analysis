
# run BUGS model script
# cf BCEA book chapter 5
#
# estimate transition probabilities
# from discrete time counts
# using multinomial likelihood



library(R2jags)

load("data/data_obs.RData")


n_S <- 3  # number of states
J <- 5    # number of time points

r.0 <-
  data_obs$`0` %>%
  filter(yr_grp == "[0,5)") %>%
  select(scd_count, non_scd_count, healthy)

r.1 <-
  data_obs$`1` %>%
  filter(yr_grp == "[0,5)") %>%
  select(scd_count, non_scd_count, healthy)

n.0 <- rowSums(r.0)
n.1 <- rowSums(r.1)

scale <- 1                            # level of informativeness for
alpha.0 <- alpha.1 <- rep(scale, n_S) # the Dirichlet prior

dataJags <-
  list("n.0", "n.1",
       "r.0", "r.1",
       "alpha.0", "alpha.1",
       "n_S")

filein <- "model.txt"
params <- c("lambda.0", "lambda.1")

inits <- function() {

  temp.0 <- matrix(
    rgamma(4*n_S, scale, 1) , 4, n_S)

  sum.temp.0 <- apply(temp.0, 1, sum)

  mat.0 <- temp.0/sum.temp.0

  temp.1 <-
    matrix(rgamma(4*n_S, scale, 1), 4, n_S)

  sum.temp.1 <- apply(temp.1, 1, sum)

  mat.1 <- temp.1/sum.temp.1

  list(lambda.0 = rbind(mat.0, rep(NA, n_S)),
       lambda.1 = rbind(mat.1, rep(NA, n_S)))
}

n.iter <- 10000
n.burnin <- 5000
n.thin <- floor((n.iter - n.burnin)/500)

mm1 <-
  jags(dataJags,
       inits,
       params,
       model.file = filein,
       n.chains = 2,
       n.iter,
       n.burnin,
       n.thin,
       DIC = TRUE)

print(mm1, digits = 3, intervals = c(0.025, 0.975))
attach.bugs(mm1$BUGSoutput)

