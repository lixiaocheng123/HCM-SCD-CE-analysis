
# Markov model script
# main_ce_analysis.R
# N Green, UCL


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
S <- length(n_init$obs)
J <- 12  # max time

# utilities and costs ----
# for each intervention

# managed with an ICD (Noyes 2007)
u_hcm <- 0.637

##TODO: what is the detriment to having an icd and scd?
u_icd <- -0.05
u_scd <- 0.5

e_unit <- list(obs =   c(u_hcm + u_icd, u_scd, 0,
                         u_hcm,         0,     0),
               risk6 = c(u_hcm + u_icd, u_scd, 0,
                         u_hcm,         0,     0))

# cost of non-fatal HCM related events £22,880 (UK Stroke Association)
c_nfatal <- 22880

##TODO: fatal cost?

# cost of SCD risk management for day case procedures_new £4792 (Waight 2019)

# EY02B NHS tariffs implantation cost
c_icd <- 4666
c_rscore <- 20

c_init <- list(obs =   c(c_icd, 0, 0,
                         0, 0, 0),
               risk6 = c(c_icd + c_rscore, 0, 0,
                         0, 0, 0))

##TODO: for shock so that its a one-off cost: extend code or tunnel state?

# what complexity & comorbidity (CC) score?
# other possible cost from NHS tariffs:
# * attention to icd £1817
# * testing icd £467

c_unit <- list(obs =   c(0, 0, c_nfatal,
                         0, 0, 0),
               risk6 = c(0, 0, c_nfatal,
                         0, 0, 0))

c_entry <- list(obs =   c(0, 0, 0,
                          0, 0, 0),
                risk6 = c(0, 0, 0,
                          0, 0, 0))

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
         c_init,
         pdecr)

# total cost across all time points
c <-
  res_new$cost %>%
  map_dfc(rowSums) %>%
  as.matrix()
e <-
  res_new$eff %>%
  map_dfc(rowSums) %>%
  as.matrix()


##TODO: add one-off cost
# cost of SCD risk algorithm £20

labels <- c("obs", "risk6")

m <-
  BCEA::bcea(e, c,
             ref = 2,
             interventions = labels,
             Kmax = 300)
plot(m)


