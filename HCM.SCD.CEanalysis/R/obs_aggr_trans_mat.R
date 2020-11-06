
#' obs_aggr_trans_mat
#'
#' Create aggregate observed transition counts
#' see p.185 Table 5.2 chapter 5 in BCEA book
#'
#' @param data new state arrivals at time step
#'             e.g. data = data_rule$`FALSE`
#' @param max_year
#' @export
#'
obs_aggr_trans_mat <- function(data,
                               max_year = 5) {

  trans_mat <- list()
  state_names <- c("healthy", "scd_count", "non_scd_count")

  # initialise counts
  trans_mat[[1]] <-
    rbind(c(0,0,0),
          c(0,0,0),
          c(0,0,0)) %>%
    `colnames<-`(state_names) %>%
    `rownames<-`(state_names)

  for (i in seq_len(max_year)) {

    trans_mat[[i + 1]] <- trans_mat[[i]]

    # new transitions
      trans_mat[[i + 1]]["healthy", state_names] <-
        trans_mat[[i + 1]]["healthy", state_names] +
        unlist(data[i, state_names])

      trans_mat[[i + 1]]["scd_count", "scd_count"] <-
      trans_mat[[i + 1]]["scd_count", "scd_count"] +
        trans_mat[[i]]["healthy", "scd_count"]

      trans_mat[[i + 1]]["non_scd_count", "non_scd_count"] <-
      trans_mat[[i + 1]]["non_scd_count", "non_scd_count"] +
        trans_mat[[i]]["healthy", "non_scd_count"]
  }

  names(trans_mat) <- c(0:max_year)
  trans_mat_total <- trans_mat[[as.character(max_year)]]

  cbind(trans_mat_total,
        n = rowSums(trans_mat_total))
}

