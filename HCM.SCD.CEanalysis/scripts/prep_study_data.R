
# prep individual level data
# to use in BUGS model
# N Green

# An ICD costs ~?10k

library(haven)
library(dplyr)


mydata <- read_dta("raw data/hcmdata.dta")

# single imputed data set
data_set1 <- mydata %>% filter(set == 1)

save(data_set1, file = "data/data_set1.RData")


n_pop <- nrow(data_set1)

probs <-
  data_set1 %>%
  mutate(non_scd = as.numeric((cvs_all & !scd) | non_cvs)) %>%
  select(time, non_scd, scd) %>%
  mutate(yr_grp = cut(time,
                      breaks = 0:35,
                      right = FALSE)) %>%
  group_by(yr_grp) %>%
  summarise(scd_count = sum(scd),
            non_scd_count = sum(non_scd)) %>%
  mutate(cum_scd = cumsum(scd_count),
         cum_non_scd = cumsum(non_scd_count),
         healthy = n_pop - cum_scd - cum_non_scd,
         at_risk = lag(healthy, default = n_pop))

probs

write.csv(probs, file = "data/probs.csv")
