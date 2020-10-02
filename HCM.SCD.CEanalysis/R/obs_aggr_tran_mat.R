
#' obs_aggr_tran_mat
#'
#' create aggregate observed transition counts
#' see Table 5.2 chapter 5 in BCEA book
#'
#' @param data new state arrivals at time step
#'             e.g. data = data_rule$`FALSE`
#'
obs_aggr_tran_mat <- function(data,
                              max_year = 5) {

  trans_mat <- list()
  state_names <- c("healthy", "scd_count", "non_scd_count")

  trans_mat[[1]] <-
    rbind(c(0,0,0),
          c(0,0,0),
          c(0,0,0)) %>%
    `colnames<-`(state_names)

  for (i in 1:max_year) {

    trans_mat[[i + 1]] <-
      trans_mat[[i]] +
      rbind(data[i, state_names],
            c(0,0,0),
            c(0,0,0)) +
      rbind(c(0,0,0),
            c(0, trans_mat[[i]][1, "scd_count"], 0),
            c(0,0, trans_mat[[i]][1, "non_scd_count"]))
  }

  names(trans_mat) <- c(0:max_year)

  cbind(
    trans_mat[[as.character(max_year)]],
    n =
      rowSums(
        trans_mat[[as.character(max_year)]]))
}

