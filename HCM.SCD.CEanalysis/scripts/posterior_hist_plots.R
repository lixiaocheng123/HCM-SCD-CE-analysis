
# posterior histograms
# transition probabilities
#


data("jags_output")
R2WinBUGS::attach.bugs(mm1$BUGSoutput)

x11()
par(mfrow = c(2,3))
hist(lambda.0[, 1, 1], breaks = 20, xlab = "Prob", main = "ICD -> ICD")
abline(v = mean(lambda.0[, 1, 1]), col = "red")
hist(lambda.0[, 1, 2], breaks = 20, xlab = "Prob", main = "ICD -> Shock")
abline(v = mean(lambda.0[, 1, 2]), col = "red")
hist(lambda.0[, 1, 3], breaks = 20, xlab = "Prob", main = "ICD -> All-cause death")
abline(v = mean(lambda.0[, 1, 3]), col = "red")

hist(lambda.0[, 4, 4], breaks = 20, xlab = "Prob", main = "Low risk -> Low risk")
abline(v = mean(lambda.0[, 4, 4]), col = "red")
hist(lambda.0[, 4, 5], breaks = 20, xlab = "Prob", main = "Low risk -> SCD")
abline(v = mean(lambda.0[, 4, 5]), col = "red")
hist(lambda.0[, 4, 6], breaks = 20, xlab = "Prob", main = "Low risk -> All-cause death")
abline(v = mean(lambda.0[, 4, 6]), col = "red")






