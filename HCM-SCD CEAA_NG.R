
# Bayesian Models for Cost-Effectiveness Analysis

# Loads packages

# s = 1. Healthy state
#     2. Stroke HCM-Related state
#     3. SCD(Sudden Cardiac Death) state
#     4. DAC (Death All Causes) state

# EXISTING METHOD OF SCD-RISK-PREDICTION (EMSRP): Current Method of SCD Risk Prediction
# HCM - SCD RISK PREDICTION MODEL(HSRPM)

S <- 4 # number of health states

tmax <- 10 # years of follow up

# load observed data on transitions for the two treatments

dataBugs <- source("HCMdata.txt")$value



list(n.0 = c(100,80,60,55), # number of patients in each state
     n.1 = c(100,80,50,55),
     eta.0 = structure(.Data = c(NA,NA,  # beta distribution parameters ~Beta(eta, tau)
                                 30,NA,
                                 40,20,
                                 25,20),
                       .Dim = c(2,4)),
     eta.1 = structure(.Data = c(NA,NA,
                                 20,NA,
                                 30,20,
                                 25,20),
                       .Dim = c(2,4)),
     tau.0 = structure(.Data = c(NA,NA,
                                 70,NA,
                                 30,60,
                                 5,40),
                       .Dim = c(2,4)),
     tau.1 = structure(.Data = c(NA,NA,
                                 80,NA,
                                 50,60,
                                 25,40),
                       .Dim = c(2,4)),
     y.0 = structure(.Data = c(2,NA,NA,NA, # number patient changing state at time step
                               NA,5,NA,NA,
                               NA,NA,6,NA,
                               NA,NA,NA,7),
                     .Dim = c(4,4)),
     y.1 = structure(.Data = c(4,NA,NA,NA,
                               NA,3,NA,NA,
                               NA,NA,4,NA,
                               NA,NA,NA,6),
                     .Dim = c(4,4)))

##TODO: I cant see how these number match up with table in report...
#status quo
c(NA,4,8,4,
  NA,NA,NA,NA,
  NA,NA,NA,NA,
  NA,NA,NA,NA)


inits1 <- source("hcminits1.txt")$value
list(lambda.0 =
       structure(.Data = c(NA,NA,NA,NA,
                           0.45,NA,NA,NA,
                           0.75,0.56,NA,NA,
                           0.79,0.65,NA,NA,
                           NA,NA,NA,NA),
                 .Dim = c(4,5)),
     lambda.1 = structure(.Data = c(NA,NA,NA,NA,
                                    0.36,NA,NA,NA,
                                    0.64,0.54,NA,NA,
                                    0.67,0.54,NA,NA,
                                    NA,NA,NA,NA),
                          .Dim = c(4,5)))

inits2 <- source("hcminits2.txt")$value
inits <- list(inits1, inits2)

# run MCMC

library(R2OpenBUGS)

model_file <- "HCMmodel.txt"
params <- c("lambda.0", "lambda.1")
n.iter <- 10000
n.burnin<- 5000
n.thin <- floor((n.iter - n.burnin)/500)
n.chain <- 2
debug <- FALSE

mm1 <- bugs(
  data = dataBugs,
  inits = inits,
  parameters.to.save = params,
  model_file = model_file,
  n.chains = n.chain,
  n.iter = n.iter,
  n.burnin = n.burnin,
  n.thin = n.thin,
  DIC = TRUE)

print(mm1, digits = 3)
attach.bugs(mm1)

# Analyse MCMC Convergence with CODA
# Check if all entries in Rhat component of Bugs output are less than 1.1
# all(mm1$summary[,"Rhat"]< 1.1)

hcm <- bugs(
  data = dataBugs,
  inits = inits,
  parameters.to.save = params,
  model_file = model_file,
  codaPkg = TRUE,
  n.chains = n.chain,
  n.iter = n.iter,
  n.burnin = n.burnin,
  n.thin = n.thin,
  DIC = TRUE)

# Convert bugs output for coda and create an mcmc object
hcm.coda <- read.bugs(hcm)

# run Markov model
start <- c(1000, 0, 0, 0)         # cohort of 1000 individuals

# NB All patients enter the model from the first state "Healthy"

# Markov transitions

## initialise population
m.0 <- m.1 <- array(NA, c(n.sims, S, tmax + 1))

for (s in 1:S) {
  m.0[, s, 1] <- start[s]
  m.1[, s, 1] <- start[s]
}

#     BUGS only outputs matrices for lambda.0 and lambda.1 with simulations for the "random" part
#     ie only the first 2 rows, as the last two are deterministically defined as c(0,0,1,1)
#     because once a patient is in SCD,and DAC, they can't move away.
#
#     So there is the need to 
#     reconstruct a full matrix with S rows and S columns for each MCMC simulations.
#     This is done by 
#     defining new arrays lambda0 and lambda1 and then stacking up the simulated values for the first (S-2)
#     rows saved in lambda.0[i,,] and lambda.1[i,,] for MCMC simulations i with a row vector
#     containing (S-2) 0s and then two 1's, ie c(0,0, 1,1)

lambda0 <- lambda1 <- array(NA, c(n.sims, S, S))

for (i in 1:n.sims) {
  lambda0[i, , ] <- rbind(rbind(lambda.0[i, , ],
                                c(0,0,1,1)),
                          c(0,0,1,1))
  lambda1[i, , ] <- rbind(rbind(lambda.1[i, , ],
                                c(0,0,1,1)),
                          c(0,0,1,1))
  
  for (j in 2:(tmax + 1)) {
    for (s in 1:S) {
      
      # lambda0,and lambda1, for the matrix multiplication
      m.0[i, s, j] <- sum(m.0[i, , j - 1]*lambda0[i, , s])
      m.1[i, s, j] <- sum(m.1[i, , j - 1]*lambda1[i, , s])
    }
  }
}







# economic analysis

utility.0 <- utility.1 <- array(NA, c(n.sims, 2, tmax)) 

dec.rate <- 0.35
utility.0[, 1, ]  <- rep(0.637, tmax)          # Healthy under t = 0
utility.1[, 1, ]  <- rep(0.637, tmax)          # Healthy under t = 1
utility.0[, 2, 1] <- dec.rate*utility.0[1, 1, 1] # Stroke-HCM Related" under t = 0
utility.1[, 2, 1] <- dec.rate*utility.1[1, 1, 1] # Stroke-HCM Related" under t = 1

for (i in 1:n.sims) {
  for (j in 2:tmax) {
    utility.0[i, 2, j] <- utility.0[i, 2, j - 1]*(1 - dec.rate)
    utility.1[i, 2, j] <- utility.1[i, 2, j - 1]*(1 - dec.rate)
  }
}

# compute QALYs accumulated under each treatment for each year of follow up

Qal.0 <- Qal.1 <- matrix(NA, n.sims, tmax)

for (i in 1:n.sims) {
  for (j in 1:tmax) {
    Qal.0[i, j] <- (m.0[i, 1, j] %*% utility.0[i, 1, j] +
                      m.0[i, 2, j] %*% utility.0[i, 2, j])/m.0[1, 1, 1]
    
    Qal.1[i, j] <- (m.1[i, 1, j] %*% utility.1[i, 1, j] +
                      m.1[i, 2, j] %*% utility.1[i, 2, j])/m.1[1, 1, 1]
  }
}

# sum values across all time points, and create matrix effectiveness

eff <- array(NA, c(n.sims, 2, tmax))
eff[, 1, ] <- apply(Qal.0, 1, sum)
eff[, 2, ] <- apply(Qal.1, 1, sum)

# annual cost for each non-fatal health state under each treatment
unit.cost.0 <- c(4792, 22880)       # Healthy, Stroke
unit.cost.1 <- c(4812, 22880) 

# Create a holding cost variable to track yearly (j > 0) accumulated cost under each treatment
cost.0 <- cost.1 <- matrix(NA, n.sims, tmax)

for (i in 1:n.sims) {
  for (j in 2:(tmax + 1)) {
    
    cost.0[i, j - 1] <-
      (m.0[i, S, j] + m.0[i, S - 1, j])*(unit.cost.0 %*% m.0[i, 1:(S-2), j])/sum(m.0[i, 1:(S-2), j]) +
      unit.cost.0 %*% m.0[i,1:(S-2),j]
    
    cost.1[i, j - 1] <-
      (m.0[i, S, j] + m.0[i, S - 1, j])*(unit.cost.1 %*% m.0[i, 1:(S-2), j])/sum(m.0[i, 1:(S-2), j]) +
      unit.cost.1 %*% m.0[i, 1:(S-2), j]
  }
}


# discount to cost and benefits
rate.b <- 0.035
rate.c <- 0.035

disc.b <- vector(mode = "numeric", length = tmax)
disc.c <- vector(mode = "numeric", length = tmax)
disc.b[1] <- 1
disc.c[1] <- 1

for (j in 2:tmax) {
  disc.b[j] <- (1 + rate.b)^(j - 1)
  disc.c[j] <- (1 + rate.c)^(j - 1)
}

disc.cost.0 <- disc.eff.0 <- matrix(NA, n.sims, tmax)
disc.cost.1 <- disc.eff.1 <- matrix(NA, n.sims, tmax)

for (j in 1:tmax) {
  disc.cost.0[, j] <- cost.0[, j]/disc.c[j]
  disc.cost.1[, j] <- cost.1[, j]/disc.c[j]
  disc.eff.0[, j]  <- eff[, 1, j]/disc.b[j]
  disc.eff.1[, j]  <- eff[, 2, j]/disc.b[j]
}

# sum the values across all time points and create matrix of costs
c <- matrix(NA, n.sims, 2)
c[, 1] <- apply(disc.cost.0, 1, sum)
c[, 2] <- apply(disc.cost.1, 1, sum)

# Sum all discounted values of effectiveness and create a matrix of discounted effectiveness
e <- matrix(NA, n.sims, 2)
e[, 1] <- apply(disc.eff.0, 1, sum)
e[, 2] <- apply(disc.eff.1, 1, sum)


# Cost-effectiveness analysis ----

library(BCEA)

ints <- c("EMSRP", "HSRPM")
m <- bcea(e, c,
          ref = 2,
          interventions = ints,
          Kmax = 25000)


#########
# plots #
#########

contour2(m,25000)
plot(m)
summary(m)


par(mfrow = c(1,2))

barplot(
  apply(m.0, c(2,3), sum),
  names.arg = seq(0,10),
  space = .2,
  xlab = "Virtual follow up",
  ylab = "Proportion of patients in each state",
  main = "EMSRP",
  ylim = c(0, 2e8))

barplot(
  apply(m.1, c(2, 3), sum),
  names.arg = seq(0, 10),
  space = .2,
  xlab = "Virtual follow up",
  ylab = "Proportion of patients in each state",
  main = "HSRPM",
  ylim = c(0, 2e8))

# trace and density for all mcmc chains
plot(hcm.coda,
     trace = TRUE,
     density = TRUE,
     smooth = FALSE)
