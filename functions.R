
# m: number of individuals [state, time, n_sim]
# n_sim: number of simulations (from jags)
# tmax: time horizon
# n_interv: number of interventions, including status-quo
#
init_pop <- function(n_init,
                     n_sim,
                     tmax,
                     n_interv = 2) {
  
  S <- length(n_init)
  pop <- list()
  
  for (n in seq_len(n_interv)) {
    pop[[n]] <- array(NA, c(S, tmax, n_sim))
    
    for (i in seq_len(n_sim)) {
      pop[[n]][, 1, i] <- n_init
    }
  }
  
  pop
}


# pop: number of individuals 
# lambda: transition probabilities
# c_unit: unit costs
# e_unit: health unit values
# delta: discount rate
#
##TODO: hard coded for 2 interventions; need to generalise
#
ce_sim <- function(pop,
                   lambda,
                   c_unit,
                   e_unit,
                   delta = 0.035) {
  
  S <- dim(pop[[1]])[1]        #number of states
  tmax <- dim(pop[[1]])[2]
  n_sim <- dim(pop[[1]])[3]
  n_interv <- length(pop)
  
  cost <- map(1:n_interv, ~matrix(NA, nrow = n_sim, ncol = tmax))
  dcost <- map(1:n_interv, ~matrix(NA, nrow = n_sim, ncol = tmax))
  eff <- map(1:n_interv, ~matrix(NA, nrow = n_sim, ncol = tmax))
  deff <- map(1:n_interv, ~matrix(NA, nrow = n_sim, ncol = tmax))
  
  for (i in seq_len(n_sim)) {
    
    for (j in seq_len(tmax)) {
      
      if (j > 1) {
        for (s in seq_len(S)) {
          pop[[1]][s, j, i] <- t(pop[[1]][, j - 1, i]) %*% lambda[[1]][, s, i]
          pop[[2]][s, j, i] <- t(pop[[2]][, j - 1, i]) %*% lambda[[2]][, s, i]
        }
      }
      
      disc <- (1 + delta)^(j-1)
      
      cost[[1]][i, j] <- c_unit[[1]] %*% pop[[1]][, j, i]
      cost[[2]][i, j] <- c_unit[[2]] %*% pop[[2]][, j, i]
      
      dcost[[1]][i, j] <- cost[[1]][i, j] / disc
      dcost[[2]][i, j] <- cost[[2]][i, j] / disc
      
      eff[[1]][i, j] <- e_unit[[1]] %*% pop[[1]][, j, i]
      eff[[2]][i, j] <- e_unit[[2]] %*% pop[[2]][, j, i]
      
      deff[[1]][i, j] <- eff[[1]][i, j] / disc
      deff[[2]][i, j] <- eff[[2]][i, j] / disc
    }
  }
  
  list(pop = pop,
       cost = cost,
       dcost = dcost,
       eff = eff,
       deff = deff)
}
