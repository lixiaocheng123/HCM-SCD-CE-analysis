
#' annual_trans_counts
#'
#' Aggregate individual data to annual totals.
#'
annual_trans_counts <- function(data,
                                cycle_length = 1) {

  n_pop <- nrow(data)

  data %>%
    # create combined non-scd field
    mutate(
      non_scd = as.numeric((cvs_all & !scd) | non_cvs)) %>%
    select(time, non_scd, scd) %>%
    mutate(yr_grp = cut(time,
                        breaks = seq(0, 35, by = cycle_length),
                        right = FALSE)) %>%
    group_by(yr_grp) %>%
    summarise(scd_count = sum(scd),
              non_scd_count = sum(non_scd)) %>%
    mutate(cum_scd = cumsum(scd_count),
           cum_non_scd = cumsum(non_scd_count),
           healthy = n_pop - cum_scd - cum_non_scd,
           at_risk = lag(healthy, default = n_pop))
}

