#'Force including certain controls in the final matched samples.
#'
#'This function processes the given distance list by adding certain zero-cost
#'edges so that the user-specified controls are forced into the final matched
#'samples. This function is of little interest to most users.
#'
#'
#'@param dist_list A distance_list object.
#'@param Z A length-n vector of treatment indicator.
#'@param include A binary vector indicating which controls must be included (length(include) = sum(1-Z).
#'
#'
#'@return  This function returns a distance list object with added edges.
#'@export

force_control <- function(dist_list, Z, include){
  n_t = sum(Z)
  # node ID
  node_id = n_t + which(include == 1)

  # Edges whose end points are not in node_id
  bad_edge_id = which(!dist_list$end_n %in% node_id)

  dist_list$d[bad_edge_id] = dist_list$d[bad_edge_id]*100

  return(dist_list)
}
