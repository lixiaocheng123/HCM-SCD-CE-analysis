

S <- 2
n_sim <- 2
p <- array(c(1,0,
             0,1,
             0,1,
             1,0),
           dim = c(S, S, 1, n_sim))
t_lim <- 10

# interv
probs = list(p, p)

# time > interv > state
c_unit <- function(t_lim) {
  function()
    purrr::map(1:t_lim,
               ~list(c(1,1), c(1,1)))
}

# interv > state
c_init <- function()
  list(c(1,1), c(1,1))

# interv > state
e_unit <- list(c(1,1), c(0,0))

# interv > time > state
pop <- list(array(c(1,0,NA,NA),
                  dim = c(S, t_lim, 2)),
            array(c(1,0,NA,NA),
                  dim = c(S, t_lim, 2)))

test_that("basic inputs", {

  ce_sim(pop = pop,
         probs = probs,
         c_unit = c_unit(t_lim),
         e_unit = e_unit,
         c_init = c_init)
})






