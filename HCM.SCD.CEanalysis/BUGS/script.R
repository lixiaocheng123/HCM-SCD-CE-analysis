
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
# pad with 0s
empty_mat <- matrix(0,
                    nrow = nrow(data_obs),
                    ncol = ncol(data_obs))

r.0 <-
  rbind(
    cbind(data_obs[[1]][, -4], empty_mat),
    cbind(empty_mat, data_obs[[2]][, -4]))

r.1 <-
  rbind(
    cbind(data_risk6[[1]][, -4], empty_mat),
    cbind(empty_mat, data_risk6[[2]][, -4]))

n_S <- nrow(r.0)  # number of states

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

  matgam.0 <- matrix(rgamma(2*n_S, scale, 1), nrow = 2)
  sum.matgam.0 <- rowSums(matgam.0)
  p.0 <- matgam.0/sum.matgam.0

  matgam.1 <- matrix(rgamma(2*n_S, scale, 1), nrow = 2)
  sum.matgam.1 <- rowSums(matgam.1)
  p.1 <- matgam.1/sum.matgam.1

  low_risk <-
    rbind(p.0,
          matrix(NA, ncol = n_S, nrow = 4))

  ICD <-
    rbind(p.1,
          matrix(NA, ncol = n_S, nrow = 4))

  # rearrange rows
  list(lambda.0 = low_risk[c(1,3,4,2,5,6), ],
       lambda.1 = ICD[c(1,3,4,2,5,6), ])
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

