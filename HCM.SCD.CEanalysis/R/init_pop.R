
#' Create initial state population
#'
#' @param n_init Initial state populations; scalar or list
#' @param n_sim Number of simulations (from jags)
#' @param tmax Time horizon
#' @param n_interv Number of interventions, including status-quo
#'
#' @return List [[interv]] state x time x sim
#' @seealso [ce_sim()]
#'
#' @export
#'
init_pop <- function(n_init,
                     n_sim,
                     tmax,
                     n_interv = 2) {

  if (is.list(n_init) & length(n_init) != n_interv)
    stop("n_interv and length of n_init don't match", call. = FALSE)

  if (length(n_init) == 1)
    n_init <- purrr::map(1:n_interv, ~n_init)

  S <- length(n_init[[1]])
  pop <- list()

  for (n in seq_len(n_interv)) {
    pop[[n]] <- array(NA, c(S, tmax, n_sim))

    for (i in seq_len(n_sim)) {
      pop[[n]][, 1, i] <- n_init[[n]]
    }
  }

  pop
}


