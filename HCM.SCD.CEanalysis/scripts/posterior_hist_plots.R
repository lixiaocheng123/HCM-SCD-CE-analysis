
# posterior histograms
# transition probabilities


data("jags_output")
R2jags::attach.jags(mm1)

# with ICD
x11()
par(mfrow = c(2,3))
# obs
hist(lambda.0[, 1, 1], breaks = 20, xlab = "Prob", main = "ICD -> ICD")
abline(v = mean(lambda.0[, 1, 1]), col = "red")
hist(lambda.0[, 1, 2], breaks = 20, xlab = "Prob", main = "ICD -> Shock")
abline(v = mean(lambda.0[, 1, 2]), col = "red")
hist(lambda.0[, 1, 3], breaks = 20, xlab = "Prob", main = "ICD -> All-cause death")
abline(v = mean(lambda.0[, 1, 3]), col = "red")
# > 6%
hist(lambda.1[, 1, 1], breaks = 20, xlab = "Prob", main = "ICD -> ICD")
abline(v = mean(lambda.1[, 1, 1]), col = "red")
hist(lambda.1[, 1, 2], breaks = 20, xlab = "Prob", main = "ICD -> Shock")
abline(v = mean(lambda.1[, 1, 2]), col = "red")
hist(lambda.1[, 1, 3], breaks = 20, xlab = "Prob", main = "ICD -> All-cause death")
abline(v = mean(lambda.1[, 1, 3]), col = "red")

# low risk
x11()
# obs
par(mfrow = c(2,3))
hist(lambda.0[, 4, 4], breaks = 20, xlab = "Prob", main = "Low risk -> Low risk")
abline(v = mean(lambda.0[, 4, 4]), col = "red")
hist(lambda.0[, 4, 5], breaks = 20, xlab = "Prob", main = "Low risk -> SCD")
abline(v = mean(lambda.0[, 4, 5]), col = "red")
hist(lambda.0[, 4, 6], breaks = 20, xlab = "Prob", main = "Low risk -> All-cause death")
abline(v = mean(lambda.0[, 4, 6]), col = "red")
# > 6%
hist(lambda.1[, 4, 4], breaks = 20, xlab = "Prob", main = "Low risk -> Low risk")
abline(v = mean(lambda.1[, 4, 4]), col = "red")
hist(lambda.1[, 4, 5], breaks = 20, xlab = "Prob", main = "Low risk -> SCD")
abline(v = mean(lambda.1[, 4, 5]), col = "red")
hist(lambda.1[, 4, 6], breaks = 20, xlab = "Prob", main = "Low risk -> All-cause death")
abline(v = mean(lambda.1[, 4, 6]), col = "red")


## ggplot2 version

library(reshape2)
library(ggplot2)
library(dplyr)

## no ICD

pobs <- melt(lambda.0[, 4, ])
prisk6 <- melt(lambda.1[, 4, ])

pobs <- cbind(pobs, strat = "obs")
prisk6 <- cbind(prisk6, strat = "risk6")

plot_dat <-
  rbind(pobs, prisk6) %>%
  filter(Var2 %in% c(4,5,6)) %>%
  mutate(Var2 =
           ifelse(Var2 == 4, "HCM no ICD",
                  ifelse(Var2 == 5, "SCD",
                         "All-cause death")))

ggplot(plot_dat, aes(x = value, y = ..density..)) +
  facet_grid(strat ~ Var2, scales = "free") +
  geom_histogram() + ggtitle("No ICD")

ggsave(filename = "images/post_hist_noICD.png")

## with ICD

pobs <- melt(lambda.0[, 1, ])
prisk6 <- melt(lambda.1[, 1, ])

pobs <- cbind(pobs, strat = "obs")
prisk6 <- cbind(prisk6, strat = "risk6")

plot_dat <-
  rbind(pobs, prisk6) %>%
  filter(Var2 %in% c(1,2,3)) %>%
  mutate(Var2 =
           ifelse(Var2 == 1, "ICD",
                  ifelse(Var2 == 2, "Shock",
                         "All-cause death")))

ggplot(plot_dat, aes(x = value, y = ..density..)) +
  facet_grid(strat ~ Var2, scales = "free") +
  geom_histogram() + ggtitle("With ICD")

ggsave(filename = "images/post_hist_withICD.png")

