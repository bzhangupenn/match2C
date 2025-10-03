#'Create a treate-to-control network to be solved via a network flow algorithm.
#'
#'This function takes in a list representation of distance matrix
#'and create a network structure to be solved.
#'
#'dist_list is a list consisting of the following three elements:
#'start_n:  the starting nodes for all edges,
#'end_n:    the ending nodes for all edges,
#'d: distance of all treated-control edges.
#' Function create_dist_list in this package constructs
#' such a list representation given a user-specified distance function.
#'
#'@param n_t Number of treated subjects.
#'@param n_c Number of controls.
#'@param dist_list A list representation of the distance matrix.
#'@param controls Number of controls matched to each treated.
#'
#'
#'@return  This function returns a list of five vectors:
#'startn, endn, ucap, cost, b.
#'
#'
#'@export


treated_control_net <- function(n_t, n_c, dist_list, controls = 1){

  start_n = dist_list$start_n
  end_n = dist_list$end_n
  d = dist_list$d

  # Sanity check I: lengths of start_n, end_n, and d need
  # to be the same.
  stopifnot(length(start_n) == length(end_n) && length(end_n) == length(d))

  # Create edges from each control to the sink, and from the source node
  # to each treated.

  start_n = c(rep(1, n_t), start_n + 1, seq(n_t + 2, n_t + n_c + 1, 1))
  end_n = c(seq(2, n_t + 1, 1), end_n + 1, rep(n_t + n_c + 2, n_c))
  d = c(rep(0, n_t), d, rep(0, n_c))
  num_edge = length(d)
  cap = c(rep(controls, n_t), rep(1, num_edge - n_t - n_c), rep(1, n_c))

  net = list(startn = start_n,
             endn = end_n,
             ucap = cap,
             cost = d,
             b =  c(controls*n_t, rep(0, n_t), rep(0, n_c), - controls*n_t))
}

