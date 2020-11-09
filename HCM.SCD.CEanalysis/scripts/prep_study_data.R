
# HCM SCD risk prediction:
# prep individual level data to use in BUGS model
# N Green
#
# the individual-level dataset is split into those who
# are given the intervention (ICD) and those who are not
# for different decision rules.
#
# the individual data are then summed into annual counts
# from the health state.
#
# then aggregated across some time horizon and formatted
# as a transition matrix


library(haven)
library(dplyr)
library(purrr)

## STATA data
# mydata <- read_dta("raw data/hcmdata.dta")
# # single imputed data set
# data_set1 <- mydata %>% filter(set == 1)
# save(data_set1, file = "data/data_set1.RData")

data("data_set1")

CYCLE <- 1 #year

# non-deterministic decision rule
FUZZY_RISK <- TRUE

fuzzy_noise <-
  if (FUZZY_RISK) {
    rnorm(nrow(data_set1), 0, 0.001)
  } else {0}


# subset cohort
data_set1 <-
  data_set1 %>%
  mutate(risk_over_6 =
           data_set1$risk_5_years + fuzzy_noise > 0.06,
         risk_over_4 =
           data_set1$risk_5_years + fuzzy_noise > 0.04)

# cox risk score > 6% ----

data_risk6 <-
  data_set1 %>%
  mutate(risk_status = ifelse(risk_over_6,
                              yes = "ICD",
                              no = "low_risk"),
         risk_status = as.factor(risk_status)) %>%
  split(.$risk_status) %>%
  map(.f = annual_trans_counts,
      cycle_length = CYCLE) %>%
  map(.f = obs_aggr_trans_mat)

data_risk6

save(data_risk6, file = "data/trans_counts_risk6.RData")

# cox risk score > 4% ----

data_risk4 <-
  data_set1 %>%
  mutate(risk_status = ifelse(risk_over_4,
                              yes = "ICD",
                              no = "low_risk"),
         risk_status = as.factor(risk_status)) %>%
  split(.$risk_status) %>%
  map(.f = annual_trans_counts,
      cycle_length = CYCLE) %>%
  map(.f = obs_aggr_trans_mat)

data_risk4

save(data_risk4, file = "data/trans_counts_risk4.RData")

# status-quo as risk factor rule ----
#
# mwt30: max wall thickness
# nsvt: NSVT
# fhxscd: family history of SCD
# syncope: unexplained syncope
#
# no noise on decision

rf_threshold <- 1

data_set1 <-
  data_set1 %>%
  mutate(num_rf = mwt30 + nsvt + fhxscd + syncope,
         rule_icd = num_rf > rf_threshold)

data_rule <-
  data_set1 %>%
  mutate(risk_status = ifelse(rule_icd,
                              yes = "ICD",
                              no = "low_risk"),
         risk_status = as.factor(risk_status)) %>%
  split(.$risk_status) %>%
  map(.f = annual_trans_counts,
      cycle_length = CYCLE) %>%
  map(.f = obs_aggr_trans_mat)

data_rule

save(data_rule, file = "data/trans_counts_rule.RData")


# status-quo as what _actually_ happens ----
# i.e. observed ICD implants

data_obs <-
  data_set1 %>%
  mutate(risk_status = ifelse(icd,
                              yes = "ICD",
                              no = "low_risk"),
         risk_status = as.factor(risk_status)) %>%
  split(.$risk_status) %>%
  map(.f = annual_trans_counts,
      cycle_length = CYCLE) %>%
  map(.f = obs_aggr_trans_mat)

data_obs

save(data_obs, file = "data/trans_counts_obs.RData")

