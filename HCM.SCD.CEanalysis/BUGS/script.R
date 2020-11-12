
# run BUGS model script
# cf BCEA book chapter 5
#
# estimate transition probabilities
# from discrete time counts
# using multinomial likelihood
#
# transition matrix:
# 1 ICD      x x x 0 0 0
# 2 shock    1 0 0 0 0 0
# 3 death    0 0 1 0 0 0
# 4 low risk 0 0 0 x x x
# 5 scd      0 0 0 0 1 0
# 6 death    0 0 0 0 0 1


library(R2jags)

data("trans_counts_obs")
data("trans_counts_risk6")


n_S <- nrow(data_obs$ICD)  # number of states

# combine low and high risk
# into single transition matrix
# pad with 0s
empty_mat <- matrix(0, nrow = n_S, ncol = n_S)

r.0 <-
  rbind(
    cbind(data_obs$ICD[, 1:n_S], empty_mat),
    cbind(empty_mat, data_obs$low_risk[, 1:n_S]))

r.1 <-
  rbind(
    cbind(data_risk6$ICD[, 1:n_S], empty_mat),
    cbind(empty_mat, data_risk6$low_risk[, 1:n_S]))

n.0 <- unname(rowSums(r.0))
n.1 <- unname(rowSums(r.1))

scale <- 1                            # level of informativeness for
alpha.0 <- alpha.1 <- rep(scale, 2*n_S) # the Dirichlet prior

dataJags <-
  list(n.0 = n.0,
       n.1 = n.1,
       r.0 = r.0,
       r.1 = r.1,
       from_shock = c(1,0,0,0,0,0),
       from_ICD_death = c(0,0,1,0,0,0),
       alpha.0 = alpha.0,
       alpha.1 = alpha.1)

filein <- "BUGS/model.txt"
params <- c("lambda.0", "lambda.1") #probabilities

#
inits <- function() {

  matgam.0 <- matrix(rgamma(2*2*n_S, scale, 1), nrow = 2)
  sum.matgam.0 <- rowSums(matgam.0)
  p.0 <- matgam.0/sum.matgam.0

  matgam.1 <- matrix(rgamma(2*2*n_S, scale, 1), nrow = 2)
  sum.matgam.1 <- rowSums(matgam.1)
  p.1 <- matgam.1/sum.matgam.1

  obs_mat <-
    rbind(p.0,
          matrix(NA, ncol = 2*n_S, nrow = 2*n_S - 2))

  risk6_mat <-
    rbind(p.1,
          matrix(NA, ncol = 2*n_S, nrow = 2*n_S - 2))

  # rearrange rows
  new_order <- c(1,3,4,2,5,6)

  list(lambda.0 = obs_mat[new_order, ],
       lambda.1 = risk6_mat[new_order, ])
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

