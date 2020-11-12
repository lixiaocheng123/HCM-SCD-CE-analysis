
# Markov model
# main_ce_analysis.R
# N Green


library(dplyr)
library(purrr)
library(BCEA)
library(reshape2)
library(HCM.SCD.CEanalysis)

data("jags_output")
R2WinBUGS::attach.bugs(mm1$BUGSoutput)

# individual patient data
data("ipd_risk")

# 1. HCM with ICD
# 2. shock
# 3. all-cause mortality
# 4. HCM
# 5. SCD
# 6. all-cause mortality

# start state populations
init_risk6 <- table(ipd_risk$risk_over_6)
init_obs <- table(ipd_risk$icd)

n_init <- list()
n_init$obs <- c(init_obs["1"], 0, 0,         # ICD
                init_obs["0"], 0, 0)         # low risk
n_init$risk6 <- c(init_risk6["TRUE"], 0, 0,
                  init_risk6["FALSE"], 0, 0)

n_sim <- nrow(lambda.0)
S <- length(n_init)
J <- 12  # max time

# utilities and costs ----
# for each intervention

# 0.637: managed with an ICD (Noyes 2007)

e_unit <- list(obs =   c(0.637, 0, 0, 1, 0.5, 0),
               risk6 = c(0.637, 0, 0, 1, 0.5, 0))

# cost of non-fatal HCM related events £22,880 (UK Stroke Association)
# cost of SCD risk management for only day case procedures_new £4792 (Waight 2019)

##TODO: do we need a tunnel state for shock so that its a one-off cost?

c_unit <- list(obs =   c(0, 0, 22880, 0, 0, 0),
               risk6 = c(0, 0, 22880, 0, 0, 0))

# linear proportion decrease each year in state
pdecr <- list(obs =   c(0, 0, 0, 0, 0, 0),
              risk6 = c(0, 0, 0, 0, 0, 0))

# rearrange and add homogeneous time dimension
probs_empty <- array(NA, dim = c(S, S, 1, n_sim))

probs <- list(probs_empty,
              probs_empty)

probs[[1]][,,1,] <- aperm(lambda.0, c(2,3,1))
probs[[2]][,,1,] <- aperm(lambda.1, c(2,3,1))


## run Markov model
res_new <-
  init_pop(n_init, n_sim, J) %>%
  ce_sim(probs,
         c_unit,
         e_unit,
         pdecr)

c <- map_dfc(res_new$cost, rowSums) %>% as.matrix()
e <- map_dfc(res_new$eff, rowSums) %>% as.matrix()


##TODO: add one-off cost
# cost of SCD risk algorithm £20

labels <- c("obs", "risk6")

m <-
  BCEA::bcea(e, c,
             ref = 2,
             interventions = labels,
             Kmax = 300)
plot(m)
