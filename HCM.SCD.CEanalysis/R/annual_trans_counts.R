
#' annual_trans_counts
#'
#' Aggregate individual data to annual totals.
#'
#' @param data
#' @param cycle_length
#' @export
#'
annual_trans_counts <- function(data,
                                cycle_length = 1) {

  n_pop <- nrow(data)

  data %>%
    # create combined non-scd field
    mutate(
      non_scd = as.numeric((cvs_all & !scd) | non_cvs),
      cens = as.numeric(!non_scd & !scd)) %>%
    select(time, non_scd, scd, cens) %>%
    mutate(yr_grp = cut(time,
                        breaks = seq(0, 40, by = cycle_length),
                        right = FALSE)) %>%
    group_by(yr_grp) %>%
    summarise(scd_count = sum(scd),
              non_scd_count = sum(non_scd),
              cens_count = sum(cens)) %>%
    mutate(cum_scd = cumsum(scd_count),
           cum_non_scd = cumsum(non_scd_count),
           cum_cens = cumsum(cens_count),
           healthy = n_pop - cum_cens - cum_scd - cum_non_scd,
           at_risk0 = lag(healthy, default = n_pop),
           # actuarial assumption
           at_risk = round(at_risk0 - 0.5*cens_count, 0))
}

