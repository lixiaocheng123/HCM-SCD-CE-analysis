
# run BUGS model script
# cf BCEA book chapter 5
#
# estimate transition probabilities
# from discrete time counts
# using multinomial likelihood


library(R2jags)

data("trans_counts_obs")
data("trans_counts_risk6")


# combine low and high risk
# into single transition matrix
empty_mat <- matrix(rep(0,9), nrow = 3)

r.0 <-
  rbind(
    cbind(data_obs[[1]][,-4], empty_mat),
    cbind(empty_mat, data_obs[[2]][,-4]))

r.1 <-
  rbind(
    cbind(data_risk6[[1]][, -4], empty_mat),
    cbind(empty_mat, data_risk6[[2]][, -4]))

n_S <- 6  # number of states

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

save(mm1, file = "data/jags_output.RData")
saveRDS(lambda.0, file = "data/lambda0.Rds")
saveRDS(lambda.1, file = "data/lambda1.Rds")

# print(mm1, digits = 3, intervals = c(0.025, 0.975))

