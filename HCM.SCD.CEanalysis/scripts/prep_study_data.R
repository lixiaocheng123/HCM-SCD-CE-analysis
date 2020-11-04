
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

mydata <- read_dta("raw data/hcmdata.dta")

# single imputed data set
data_set1 <- mydata %>% filter(set == 1)

save(data_set1, file = "data/data_set1.RData")

CYCLE <- 1

# non-deterministic decision
fuzzy_risk_noise <- rnorm(nrow(data_set1), 0, 0.001)

# cox risk score > 6% ----

data_set1$risk_over_6 <-
  as.factor((data_set1$risk_5_years + fuzzy_risk_noise) > 0.06)

data_risk6 <-
  group_split(data_set1, risk_over_6) %>%
  setNames(levels(data_set1$risk_over_6)) %>%
  map(.f = annual_trans_counts,
      cycle_length = CYCLE) %>%
  map(.f = obs_aggr_trans_mat)

data_risk6

save(data_risk6, file = "data/data_risk6.RData")

# cox risk score > 4% ----

data_set1$risk_over_4 <-
  as.factor((data_set1$risk_5_years + fuzzy_risk_noise) > 0.04)

data_risk4 <-
  group_split(data_set1, risk_over_4) %>%
  setNames(levels(data_set1$risk_over_4)) %>%
  map(.f = annual_trans_counts,
      cycle_length = CYCLE)

data_risk4

save(data_risk4, file = "data/data_risk4.RData")

# status-quo as risk factor rule ----
#
# mwt30: max wall thickness
# nsvt: NSVT
# fhxscd: family history of SCD
# syncope: unexplained syncope
#
# no noise on decision
#
data_set1 <-
  data_set1 %>%
  mutate(num_rf = mwt30 + nsvt + fhxscd + syncope,
         rule_icd = num_rf > 1)

data_rule <-
  group_split(data_set1, rule_icd) %>%
  setNames(unique(data_set1$rule_icd)) %>%
  map(.f = annual_trans_counts,
      cycle_length = CYCLE)

data_rule

save(data_rule, file = "data/data_rule.RData")


# status-quo as what _actually_ happens ----
# i.e. observed ICD implants
data_obs <-
  group_split(data_set1, icd) %>%
  setNames(unique(data_set1$icd)) %>%
  map(.f = annual_trans_counts,
      cycle_length = CYCLE)

data_obs

save(data_obs, file = "data/data_obs.RData")

