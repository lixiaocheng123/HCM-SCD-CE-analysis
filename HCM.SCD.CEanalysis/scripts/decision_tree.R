
#
#
#

c_ICD <- 1
p_ICD_under4
p_ICD_4to6
p_ICD_over6
p_under4
p_4to6
p_over6

# c_assess + c_ICD * (p_under4*p_ICD_under4 + p_4to6*p_ICD_4to6 + p_over6*p_ICD_over6)
#
# and values we need are:
#
# cost of assessment (c_assess)
# cost of ICD (c_ICD)
# probability assessed <4%, 4-6%, >6%
# probability ICD given assessment result
# Alive with HCM: unit cost, unit utility for ICD and non-ICD
# p(alive with HCM -> SCD) for ICD and non-ICD, by age (or other covariate?)
# p(alive with HCM -> all cause death) for ICD and non-ICD, by age (or other covariate?)


