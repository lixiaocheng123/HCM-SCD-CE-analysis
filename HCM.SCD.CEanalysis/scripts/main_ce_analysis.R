
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

# structure:
# (1.icd  2.shock 3.mortality)
# (4.~icd 5.scd   6.mortality)


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
J <- 12  # max time (year)


# utilities and costs --------
# for each intervention

u_hcm <- 0.637   # managed with an ICD (Noyes 2007)
u_icd <- -0.05   # ICD annual detriment
u_shock <- 0.5   # annual detriment
u_icd_init <- -6/365 # icd implant one-off (Smith EJH 2013)

e_unit <- list(obs =   c(u_hcm + u_icd, u_shock, 0,
                         u_hcm,         0,       0),
               risk6 = c(u_hcm + u_icd, u_shock, 0,
                         u_hcm,         0,       0))


# non-fatal HCM related events (UK Stroke Association)
c_nfatal <- 22880
c_icd <- 4666         # EY02B NHS tariffs implantation cost
c_rscore <- 20
c_icd_appt <- 10
c_icd_repl <- 45000

# implant complication
p_compl <- 0.047        # Smith EHJ 2013
p_compl_repl <- 0.032   # Smith EHJ 2013
# weighted infection, dislodgement
c_compl <- (0.0227*37116 + 0.00828*6146)/(0.0227 + 0.00828)

c_entry <- list(obs =   c(0, 0, 0,
                          0, 0, 0),
                risk6 = c(0, 0, 0,
                          0, 0, 0))

# time to replace ICD
t_repl <- 10


#' sample costs
#'
rc_unit <- function(J,
                    c_icd_appt,
                    c_icd_repl,
                    c_nfatal) {
  function() {
    rc_appt <- runif(1, c_icd_appt, c_icd_appt)

    c_unit <-
      list(obs =   c(2*rc_appt, c_nfatal, 0,
                     0, 0, 0),
           risk6 = c(2*rc_appt, c_nfatal, 0,
                     0, 0, 0))

    c_unit_t <- purrr::map(1:J, ~c_unit)

    # every 10 years replace ICD
    rc_repl <- runif(1, c_icd_repl, c_icd_repl)

    c_unit_t[[10]]$obs[1] <- c_unit_t[[10]]$obs[1] + rc_repl
    c_unit_t[[10]]$risk6[1] <- c_unit_t[[10]]$risk6[1] + rc_repl

    c_unit_t
  }
}

#'
rc_init <- function(c_icd,
                    c_rscore) {
  function() {
    rc_icd <- runif(1, c_icd, c_icd)

    list(obs =   c(c_icd + p_compl*c_compl, 0, 0,
                   0,     0, 0),
         risk6 = c(c_icd + c_rscore + p_compl*c_compl, 0, 0,
                   0,     0, 0))
  }
}


# linear proportion decrease each year in state
pdecr <- list(obs =   c(0, 0, 0, 0, 0, 0),
              risk6 = c(0, 0, 0, 0, 0, 0))

# rearrange and add homogeneous time dimension
probs_empty <- array(NA, dim = c(S, S, 1, n_sim))

probs <- list(probs_empty,
              probs_empty)

probs[[1]][,,1,] <- aperm(lambda.0, c(2,3,1))
probs[[2]][,,1,] <- aperm(lambda.1, c(2,3,1))


#############
# run model #
#############

res_new <-
  init_pop(n_init, n_sim, J) %>%
  ce_sim(probs,
         rc_unit(J,
                 c_icd_appt,
                 c_icd_repl,
                 c_nfatal),
         e_unit,
         rc_init(c_icd, c_rscore),
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


labels <- c("obs", "risk6")

m <-
  BCEA::bcea(e, c,
             ref = 2,
             interventions = labels,
             Kmax = 300)
plot(m)


