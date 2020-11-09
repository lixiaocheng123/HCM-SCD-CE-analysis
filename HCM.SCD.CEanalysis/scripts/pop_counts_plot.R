
# pop counts plot
# stacked percentage bar plot


library(ggplot2)
library(gridExtra)

plot_obs <-
  apply(res_new$pop[[1]], c(1,2), mean) %>%
  as.data.frame()

plot_obs$state <-
  c("ICD", "shock", "all_cause", "low_risk", "scd", "all_cause")

plot_obs <-
  melt(plot_obs,
       id.vars = "state",
       value.name = "pop",
       variable.name = "time") %>%
  group_by(state, time) %>%
  summarise(pop = sum(pop))

plot_obs$time <- as.numeric(gsub("V", "", plot_obs$time))
plot_obs$state <-
  factor(plot_obs$state,
         levels = c("scd", "shock", "all_cause", "ICD", "low_risk"))

p1 <-
  ggplot(plot_obs, aes(fill = state, y = pop, x = time)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_x_continuous(breaks = 1:12)

plot_risk6 <-
  apply(res_new$pop[[2]], c(1,2), mean) %>%
  as.data.frame()

plot_risk6$state <-
  c("ICD", "shock", "all_cause", "low_risk", "scd", "all_cause")

plot_risk6 <-
  melt(plot_risk6,
       id.vars = "state",
       value.name = "pop",
       variable.name = "time") %>%
  group_by(state, time) %>%
  summarise(pop = sum(pop))

plot_risk6$time <- as.numeric(gsub("V", "", plot_risk6$time))
plot_risk6$state <-
  factor(plot_risk6$state,
         levels = c("scd", "shock", "all_cause", "ICD", "low_risk"))

p2 <-
  ggplot(plot_risk6, aes(fill = state, y = pop, x = time)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_x_continuous(breaks = 1:12)

# combine
total_dat <- rbind(cbind(strat = "obs", plot_obs),
                   cbind(strat = "risk6", plot_risk6))

facet_plot <-
  ggplot(total_dat, aes(fill = state, y = pop, x = time)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_x_continuous(breaks = 1:12) +
  facet_grid(.~strat)

ggsave(facet_plot, file = "images/state_pop_over_time.png")


