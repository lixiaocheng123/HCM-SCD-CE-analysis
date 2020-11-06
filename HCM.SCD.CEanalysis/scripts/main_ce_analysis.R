
# Markov model
# main_ce_analysis.R
# N Green


library(dplyr)
library(purrr)
library(BCEA)
library(res_newhape2)
library(HCM.SCD.CEanalysis)

data(mm1)
R2WinBUGS::attach.bugs(mm1$BUGSoutput)

# 1. alive with HCM
# 2. SCD
# 3. all-cause mortality

n_init <- c(1000, 0, 0)
n_sim <- nrow(lambda.0)
S <- length(n_init)
J <- 12  # max time

# utilities and costs
# for each intervention

# 0.637: managed with an ICD (Noyes 2007)

e_unit <- list(low_risk = c(1, 0, 0),
               ICD = c(0.637, 0.5, 0))

# cost of SCD risk algorithm £20
# cost of SCD risk management for only day case procedures_new £4792 (Waight 2019)
# cost of non-fatal HCM related events £22,880 (UK Stroke Association)

c_unit <- list(low_risk = c(0, 0, 0),
               ICD = c(0, 22880, 0))

# linear proportion decrease each year in state
pdecr <- list(low_risk = c(0, 0, 0),
              ICD = c(0, 0, 0))

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

# total cost for this risk stratification
total_cost <- res_new$cost[[1]] + res_new$cost[[2]]
total_eff <- res_new$eff[[1]] + res_new$eff[[2]]

# sum across all time points

c <- rowSums(total_cost)
e <- rowSums(total_eff)

# c <- map_dfc(res_new$cost, rowSums) %>% as.matrix()
# e <- map_dfc(res_new$eff, rowSums) %>% as.matrix()

labels <- c("low_risk", "ICD")

m <-
  BCEA::bcea(e, c,
             ref = 2,
             interventions = labels,
             Kmax = 300)
plot(m)
