
# run BUGS model script
# cf BCEA book chapter 5
#
# estimate transition probabilities
# from discrete time counts
# using multinomial likelihood


library(R2jags)

# select:
data("trans_counts_obs")
data("trans_counts_risk6")


n_S <- 3  # number of states
J <- 5    # number of time points

r.0 <- data_obs$`0`[, c("healthy", "scd_count", "non_scd_count")]
r.1 <- data_obs$`1`[, c("healthy", "scd_count", "non_scd_count")]

n.0 <- unname(rowSums(r.0))
n.1 <- unname(rowSums(r.1))

scale <- 1                            # level of informativeness for
alpha.0 <- alpha.1 <- rep(scale, n_S) # the Dirichlet prior

dataJags <-
  list(n.0 = n.0,
       n.1 = n.1,
       r.0 = r.0,
       r.1 = r.1,
       alpha.0 = alpha.0,
       alpha.1 = alpha.1,
       n_S = n_S)

filein <- "BUGS/model.txt"
params <- c("lambda.0", "lambda.1")

#
inits <- function() {

  temp.0 <- rgamma(n_S, scale, 1)
  sum.temp.0 <- sum(temp.0)
  p.0 <- temp.0/sum.temp.0

  temp.1 <- rgamma(n_S, scale, 1)
  sum.temp.1 <- sum(temp.1)
  p.1 <- temp.1/sum.temp.1

  list(lambda.0 = rbind(p.0,
                        rep(NA, n_S),
                        rep(NA, n_S)),
       lambda.1 = rbind(p.1,
                        rep(NA, n_S),
                        rep(NA, n_S)))
}

inits()

n.iter <- 10000
n.burnin <- 5000
n.thin <- floor((n.iter - n.burnin)/500)

mm1 <-
  jags(data = dataJags,
       inits = inits,
       parameters.to.save = params,
       model.file = filein,
       n.chains = 2,
       n.iter,
       n.burnin,
       n.thin,
       DIC = TRUE)

R2WinBUGS::attach.bugs(mm1$BUGSoutput)

#select:
save(mm1, file = "data/jags_obs.RData")
saveRDS(lambda.0, file = "data/lambda0_obs.Rds")
saveRDS(lambda.1, file = "data/lambda1_obs.Rds")

save(mm1, file = "data/jags_risk6.RData")
saveRDS(lambda.0, file = "data/lambda0_risk6.Rds")
saveRDS(lambda.1, file = "data/lambda1_risk6.Rds")

# print(mm1, digits = 3, intervals = c(0.025, 0.975))

