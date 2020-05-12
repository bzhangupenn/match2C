#'Solve a network flow problem.
#'
#'This function solves network flow optimization problems
#'by calling the RELAX-IV algorithm implemented in FORTRAN
#'by Dimitri Bertsekas and Paul Tseng, and made available by
#'Sam Pimentel in the package rcbalance. This function is of
#'limited interest to users.
#'
#'
#'
#'@param net A list of five vectors: startn, endn, ucap, cost, b.
#'
#'
#'@return If the problem is feasible, function returns a list with the following elements:
#'crash: an integer, equal to zero if the algorithm ran correctly and equal to 1 if it crashed.
#'feasible: an integer, equal to zero if the problem is not feasible.
#'x: a vector equal in length to the number of arcs in argument problem net,
#'giving in each coordinate the number of units of flow passing across the
#'corresponding edge in the optimal network flow.
#'If the problem is not feasible, it returns "Not feasible."
#'
#'
#'@importFrom rcbalance callrelax
#'@export


solve_network_flow <- function(net){
  'Solve a network flow problem using callrelax.
  '
  res = rcbalance::callrelax(net)
  if (res$feasible == 0) return('Not feasible')
  return(res)
}
