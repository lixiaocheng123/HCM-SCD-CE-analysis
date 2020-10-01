
# prep individual level data
# to use in BUGS model
# N Green


#' aggregate individual data to annual totals
#'
annual_trans_counts <- function(data,
                                cycle_length = 1) {

  n_pop <- nrow(data)

  data %>%
    # create combined non-scd field
    mutate(
      non_scd = as.numeric((cvs_all & !scd) | non_cvs)) %>%
    select(time, non_scd, scd) %>%
    mutate(yr_grp = cut(time,
                        breaks = seq(0, 35, by = cycle_length),
                        right = FALSE)) %>%
    group_by(yr_grp) %>%
    summarise(scd_count = sum(scd),
              non_scd_count = sum(non_scd)) %>%
    mutate(cum_scd = cumsum(scd_count),
           cum_non_scd = cumsum(non_scd_count),
           healthy = n_pop - cum_scd - cum_non_scd,
           at_risk = lag(healthy, default = n_pop))
}

# -------------------------------------------------------------------------

library(haven)
library(dplyr)
library(purrr)

mydata <- read_dta("raw data/hcmdata.dta")

# single imputed data set
data_set1 <- mydata %>% filter(set == 1)

save(data_set1, file = "data/data_set1.RData")


# cox risk score > 6% ----

data_set1$risk_over_6 <- as.factor(data_set1$risk_5_years > 0.06)

data_risk6 <-
  group_split(data_set1, risk_over_6) %>%
  setNames(levels(data_set1$risk_over_6)) %>%
  map(.f = annual_trans_counts,
      cycle_length = 5)

data_risk6

save(data_risk6, file = "data/data_risk6.RData")

# cox risk score > 4% ----

data_set1$risk_over_4 <- as.factor(data_set1$risk_5_years > 0.04)

data_risk4 <-
  group_split(data_set1, risk_over_4) %>%
  setNames(levels(data_set1$risk_over_4)) %>%
  map(.f = annual_trans_counts,
      cycle_length = )

data_risk4

save(data_risk4, file = "data/data_risk4.RData")

# status-quo as risk factor rule ----
#
# mwt30: max wall thickness
# nsvt: NSVT
# fhxscd: family history of SCD
# syncope: unexplained syncope
#
data_set1 <-
  data_set1 %>%
  mutate(num_rf = mwt30 + nsvt + fhxscd + syncope,
         rule_icd = num_rf > 1)

data_rule <-
  group_split(data_set1, rule_icd) %>%
  setNames(unique(data_set1$rule_icd)) %>%
  map(.f = annual_trans_counts,
      cycle_length = 5)

data_rule

save(data_rule, file = "data/data_rule.RData")


# status-quo as what _actually_ happens ----
# i.e. observed ICD implants
data_obs <-
  group_split(data_set1, icd) %>%
  setNames(unique(data_set1$icd)) %>%
  map(.f = annual_trans_counts,
      cycle_length = 5)

data_obs

save(data_obs, file = "data/data_obs.RData")

