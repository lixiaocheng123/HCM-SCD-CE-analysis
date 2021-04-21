
library(dplyr)
library(purrr)


pop_init <- c(1,0)

pop <- array(pop_init,
             dim = c(S, t_lim))

S <- length(pop_init)

probs <- array(c(0,1,
                 1,0),
               dim = c(S, S, 1))
t_lim <- 10

# time=dependent
c_unit <- function(t_lim) {
  function()
    purrr::map(1:t_lim, ~ c(1,1))
}

c_init <- function() c(1,1)

# interv > state
e_unit <- c(1,0)


list(pop = pop,
     probs = probs,
     c_unit = c_unit(t_lim),
     c_init = c_init,
     e_unit = e_unit) %>%
ce_markov()

test_that("basic inputs", {
  # expect_equal()
})





