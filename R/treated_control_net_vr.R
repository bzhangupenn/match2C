#'Create a treate-to-control network to be solved via a network flow algorithm.
#'
#'
#'dist_list is a list consisting of the following three elements:
#'start_n:  the starting nodes for all edges,
#'end_n:    the ending nodes for all edges,
#'d:        distance of all treated-control edges.
#'
#'
#'@param n_t Number of treated subjects.
#'@param n_c Number of controls.
#'@param dist_list A list representation of the distance matrix.
#'@param maxcontrols Maximum number of controls matched to each treated.
#'@param mincontrols Minimum number of controls matched to each treated.
#'
#'
#'@return  This function returns a list of five vectors:
#'startn, endn, ucap, cost, b.
#'
#'
#'@export
#'
treated_control_net_vr <- function(n_t, n_c, dist_list, maxcontrols = NULL, mincontrols = 1){
  start_n = dist_list$start_n
  end_n = dist_list$end_n
  d = dist_list$d
  
  # lengths of start_n, end_n, and d need to be the same.
  stopifnot(length(start_n) == length(end_n) && length(end_n) == length(d))
  
  # The maximum number of control units to be matched to one treated unit
  if (is.null(maxcontrols)) maxcontrols = n_c
  
  # Create edges from each control to the sink, from the source node
  # to each treated and each added, from each treated to the Overflow
  ta <- unique(start_n)
  source_uf <- length(unique(start_n))
  sink_index <- max(c(start_n, end_n)) + 2
  start_n = c(rep(1, source_uf), start_n + 1, seq(n_t + 2, n_t + n_c + 1, 1), seq(2,n_t + 1, 1))
  end_n = c( ta + 1, end_n + 1, rep(sink_index, n_c), rep(sink_index+1, n_t))
  
  # The distance for each edge
  d = c(rep(0, source_uf), d, rep(0, n_c), rep(0, n_t))
  num_edge = length(d)
  
  # The upper capacity for each edge (the lower capacity are all zero)
  cap = c(rep(maxcontrols, n_t), rep(1, source_uf - n_t), 
          rep(1, num_edge - source_uf - n_t - n_c), 
          rep(1, n_c), rep(maxcontrols - mincontrols, n_t))
  
  # The supply or demand of each node 
  b = c(maxcontrols*n_t + source_uf - n_t, rep(0, n_t + n_c), rep(0, source_uf - n_t), - n_c, 
        - maxcontrols*n_t - source_uf + n_t + n_c)
  
  net = list(startn = start_n,
             endn = end_n,
             ucap = cap,
             cost = d,
             b =  b)
}
