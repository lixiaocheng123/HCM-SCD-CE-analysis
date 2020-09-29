
# Markov model
# main_analysis.R
# N Green


library(dplyr)
library(purrr)
library(BCEA)
library(reshape2)
library(HCM.SCD.CEanalysis)


# 1. alive with HCM
# 2. SCD
# 3. all-cause mortality

n_init <- c(1000, 0, 0)
n_sim <- 2
S <- length(n_init)
J <- 12  # max time

# utilities and costs
e_unit <- list(c(0.637, 0, 0),
               c(0.637, 0, 0))

c_unit <- list(c(4792, 0, 0),
               c(4812, 0, 0))

# linear proportion decrease each year in state
pdecr <- list(c(0, 0, 0),
              c(0, 0, 0))

probs <- list(array(c(0.90, 0.075, 0.025,
                      0, 1, 0,
                      0, 0, 1),               # status-quo
                    dim = c(S, S, 1, n_sim)),
              array(c(0.90, 0.075, 0.025,
                      0, 1, 0,
                      0, 0, 1),               # ICD
                    dim = c(S, S, 1, n_sim)))

# rearrange because array() fills by column
probs <- map(probs, aperm, perm = c(2, 1, 3, 4))

## run Markov model
res <-
  init_pop(n_init, n_sim, J) %>%
  ce_sim(probs,
         c_unit,
         e_unit,
         pdecr)

##TODO:
# covariate-dependent transition probs?
# c_unit, e_unit samples?


# sum across all time points

c <- map_dfc(res$cost, rowSums) %>% as.matrix()
e <- map_dfc(res$eff, rowSums) %>% as.matrix()

labels <- c("status-quo", "ICD")

m <- bcea(e, c,
          ref = 2,
          interventions = labels,
          Kmax = 300)

