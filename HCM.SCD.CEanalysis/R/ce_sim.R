
#' Markov model cost-effectiveness simulation
#'
#' @param pop Number of individuals. List-matrix [[n_int]] n_sim x time
#' @param probs Transition probabilities. List-array [[n_int]] states x states x time x sim
#' @param c_unit Unit costs per state. List [[n_int]]
#' @param e_unit Health unit values per state. List [[n_int]]
#' @param pdecr Linear utility decrease per state
#' @param delta Discount rate; default 3.5\% as 0.035
#' @return List
#' @seealso [init_pop()]
#'
#' @importFrom purrr map
#'
#' @export
#'
ce_sim <- function(pop,
                   probs,
                   c_unit,
                   e_unit,
                   pdecr = NA,
                   delta = 0.035) {

  stopifnot(delta >= 0, delta <= 1)

  S <- dim(pop[[1]])[1]        # number of states
  tmax <- dim(pop[[1]])[2]     # time horizon
  n_sim <- dim(pop[[1]])[3]
  n_interv <- length(pop)

  if (any(is.na(pdecr))) {
    pdecr <- map(1:n_interv, ~rep(0, S))}

  # initialise empty output matrices
  out_mat <- map(1:n_interv,
                 ~matrix(NA,
                         nrow = n_sim,
                         ncol = tmax))
  cost <- out_mat
  eff <- out_mat
  dcost <- out_mat
  deff <- out_mat

  for (k in seq_len(n_interv)) {
    for (i in seq_len(n_sim)) {
      for (j in seq_len(tmax)) {

        # time-homogeneous dim 1
        t <- min(dim(probs[[1]])[3], j - 1)

        if (j > 1) {
          for (s in seq_len(S)) {

            pop[[k]][s, j, i] <-
              t(pop[[k]][, j - 1, i]) %*% probs[[k]][, s, t, i]
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

