
#' Create initial state population
#'
#' @param pop Number of individuals (state, time, n_sim)
#' @param n_sim Number of simulations (from jags)
#' @param tmax Time horizon
#' @param n_interv Number of interventions, including status-quo
#' @return List [[interv]] state x time x sim
#' @seealso [ce_sim()]
#'
#' @export
#'
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
