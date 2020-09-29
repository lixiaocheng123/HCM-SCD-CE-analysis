
# run BUGS model script
# cf BCEA book chapter 5
#


library(R2jags)

S <- 5  # number of states
J <- 12 # number of time points

# observed cases for t=0
r.0 <- matrix(c(
  66, 32, 0, 0, 2,
  42, 752,0, 5, 20,
  0,  4,  0, 1, 0),
  c(3, S),
  byrow = TRUE)

# observed cases for t=1
r.1 <- matrix(c(
  210, 60,  0, 1, 1,
  88,  641, 0, 4, 13,
  1,   0,   0, 0, 1),
  c(3, S),
  byrow = TRUE)

n.0 <- apply(r.0, 1, sum) # number of patients in
n.1 <- apply(r.1, 1, sum) # each state for t=0,1
scale <- 1                # level of informativeness for
alpha.0 <- alpha.1 <- rep(scale, ) # the Dirichlet prior


dataJags <-
  list("n.0", "n.1", "r.0", "r.1", "alpha.0", "alpha.1", "S")

filein <- "model.txt"
params <- c("lambda.0","lambda.1")

inits <- function(){
  temp.0 <- matrix(rgamma(4*S,scale,1),4,S)
  sum.temp.0 <- apply(temp.0,1,sum)
  mat.0 <- temp.0/sum.temp.0
  temp.1 <- matrix(rgamma(4*S,scale,1),4,S)
  sum.temp.1 <- apply(temp.1,1,sum)
  mat.1 <- temp.1/sum.temp.1
  list(lambda.0 = rbind(mat.0,rep(NA,S)),
       lambda.1 = rbind(mat.1,rep(NA,S))
  )
}

n.iter <- 10000
n.burnin <- 5000
n.thin <- floor((n.iter-n.burnin)/500)
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

print(mm1, digits = 3,intervals = c(0.025, 0.975))
attach.bugs(mm1$BUGSoutput)

