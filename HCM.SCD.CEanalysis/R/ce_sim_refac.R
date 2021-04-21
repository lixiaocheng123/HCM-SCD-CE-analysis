
##TODO...

# ce_sim <- function(pop,
#                    probs,
#                    c_unit,
#                    e_unit,
#                    c_init = NA,
#                    pdecr = NA,
#                    delta = 0.035) {
#
#   stopifnot(delta >= 0, delta <= 1)
#
#   ##TODO: check inputs
#
#   ##TODO: interleave pop, probs, c_unit, c_init. e_unit
#   # so that list over interventions and sim i
#   # inputs <-
#
#   out <-
#     map(inputs,
#         ~ map(.x, ~ce_markov))
#
#   # reorder
#   # purrr::
#
#   list(pop = pop,
#        cost = cost,
#        dcost = dcost,
#        eff = eff,
#        deff = deff)
# }
#
#
# #
# next_pop <- function(pop, j, t) {
#   function(s) {
#     t(pop[, j - 1]) %*% probs[, s, t]
#   }
# }


#
ce_markov <- function(inputs,
                      pdecr = NA,
                      delta = 0.035) {

  pop <- inputs$pop
  probs <- inputs$probs
  c_unit <- inputs$c_unit
  c_init <- inputs$c_init
  e_unit <- inputs$e_unit

  S <- dim(pop)[1]        # number of states
  tmax <- dim(pop)[2]     # time horizon

  if (any(is.na(pdecr))) {
    pdecr <- rep(0, S)
  }

  # sample costs
  rc_unit <- c_unit()

  cost <- vector("numeric", tmax)
  eff <- vector("numeric", tmax)
  dcost <- vector("numeric", tmax)
  deff <- vector("numeric", tmax)

  for (j in seq_len(tmax)) {

    # time-homogeneous dim 1
    t <- min(dim(probs)[3], j - 1)

    if (j > 1) {
      for (s in seq_len(S)) {

        pop[s, j] <-
          t(pop[, j - 1]) %*% probs[, s, t]
      }
    }

    disc <- (1 + delta)^(j-1)

    e_state <- e_unit*(1 - pdecr)^j

    cost[j] <- rc_unit[[j]] %*% pop[, j]

    eff[j] <- e_state %*% pop[, j]

    dcost[j] <- cost[j] / disc
    deff[j] <- eff[j] / disc
  }

  # add one-off starting state costs
  if (is.function(c_init)) {
    cost[1] <-
      cost[1] + c_init() %*% pop[, 1]

    dcost[1] <-
      dcost[1] + c_init() %*% pop[, 1]
  }

  list(pop = pop,
       cost = cost,
       dcost = dcost,
       eff = eff,
       deff = deff)
}

