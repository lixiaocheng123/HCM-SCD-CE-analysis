
# main_analysis.R
# N Green


library(dplyr)
library(purrr)
library(BCEA)
library(reshape2)
library(HCM.SCD.CEanalysis)


# 1. healthy
# 2. stroke
# 3. SCD
# 4. all-cause death

n_init <- c(1000, 0, 0, 0)
n_sim <- 2
S <- length(n_init)
J <- 12

# initial utility
e_unit <- list(c(0.637, 0.637, 0, 0),
               c(0.637, 0.637, 0, 0))

# linear proportion decrease each year in state
pdecr <- list(c(0, 0.35, 0, 0),
              c(0, 0.35, 0, 0))

lambda <- list(array(c(0.90, 0.05, 0.025, 0.025,
                       0.90, 0.05, 0.025, 0.025,
                       0.90, 0.05, 0.025, 0.025,
                       0.90, 0.05, 0.025, 0.025),   # status-quo
                     dim = c(S, S, n_sim)),
               array(c(0.90, 0.05, 0.025, 0.025,
                       0.90, 0.05, 0.025, 0.025,
                       0.90, 0.05, 0.025, 0.025,
                       0.90, 0.05, 0.025, 0.025),   # ICD
                     dim = c(S, S, n_sim)))

# because array() fills by column
lambda <- map(lambda, aperm, perm = c(2,1,3))

c_unit <- list(c(4792, 22880, 0, 0),
               c(4812, 22880, 0, 0))

## run Markov model
res <-
  init_pop(n_init, n_sim, J) %>%
  ce_sim(lambda,
         c_unit,
         e_unit,
         pdecr)


# sum across all time points

c <- map_dfc(res$cost, rowSums) %>% as.matrix()
e <- map_dfc(res$eff, rowSums) %>% as.matrix()

labels <- c("status-quo", "ICD")

m <- bcea(e, c,
          ref = 2,
          interventions = labels,
          Kmax = 300)
