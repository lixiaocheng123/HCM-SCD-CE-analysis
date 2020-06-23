
# Create initial population
#
# pop: number of individuals [state, time, n_sim]
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


# Markov model cost-effectiveness simulation
#
# pop: number of individuals 
# lambda: transition probabilities
# c_unit: unit costs per state
# e_unit: health unit values per state
# pdecr: linear utility decrease per state
# delta: discount rate
#
ce_sim <- function(pop,
                   lambda,
                   c_unit,
                   e_unit,
                   pdecr = NA,
                   delta = 0.035) {
  
  stopifnot(delta >= 0, delta <= 1)
  
  S <- dim(pop[[1]])[1]        # number of states
  tmax <- dim(pop[[1]])[2]     # time horizon
  n_sim <- dim(pop[[1]])[3]
  n_interv <- length(pop)
  
  if (is.na(pdecr)) {
    pdecr <- map(1:n_interv, ~rep(0, S))}
  
  # initialise empty output matrices
  out_mat <- map(1:n_interv, ~matrix(NA,
                                     nrow = n_sim,
                                     ncol = tmax))
  cost <- out_mat
  eff <- out_mat
  dcost <- out_mat
  deff <- out_mat
  
  for (k in seq_len(n_interv)) {
    for (i in seq_len(n_sim)) {
      for (j in seq_len(tmax)) {
        
        if (j > 1) {
          for (s in seq_len(S)) {
            
            pop[[k]][s, j, i] <-
              t(pop[[k]][, j - 1, i]) %*% lambda[[k]][, s, i]
          }
        }
        
        disc <- (1 + delta)^(j-1)
        
        e_state <- e_unit[[k]]*(1 - pdecr[[k]])^j
        
        cost[[k]][i, j] <- c_unit[[k]] %*% pop[[k]][, j, i]
        eff[[k]][i, j] <- e_state %*% pop[[k]][, j, i]
        
        dcost[[k]][i, j] <- cost[[k]][i, j] / disc
        deff[[k]][i, j] <- eff[[k]][i, j] / disc
      }
    }
  }
  
  list(pop = pop,
       cost = cost,
       dcost = dcost,
       eff = eff,
       deff = deff)
}
