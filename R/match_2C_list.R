#'Perform a pair matching using two user-specified list representations of
#'distance matrices.
#'
#'This function performs a pair-matching using (at most) two user-specified distance
#'matrices in their (possibly sparse) list representations. For more details
#'on ``list representations'' of a treatment-by-control distance matrix, see
#'the documentation of the function ``create_list_from_mat''.
#'
#'This function is designed for more experienced and sophisticated R users.
#'Instead of providing possibly dense treatment-by-control distance matrices
#'that take up a lot of memories, users may simply provide two lists that specifies
#'informations of edges: their starting points, ending points, capacity, and cost.
#'For more information on list representations of a distance matrix, see the
#'documentation of the function ``create_list_from_mat'' and ``create_list_from_scratch''.
#'Note that by setting dist_list_2 = NULL, the usual matching framework is restored.
#'
#'@param Z A length-n vector of treatment indicator.
#'@param dataset dataset to be matched.
#'@param dist_list_1 A (possibly sparse) list representation of
#'                   treatment-by-control distance matrix.
#'@param dist_list_2 A second (possibly sparse) list representation of
#'                   treatment-by-control distance matrix.
#'@param lambda A penalty that does a trade-off between two parts of the network.
#'@param controls Number of controls matched to each treated. Default is set to 1.
#'@param overflow A logical value indicating if overflow protection is turned on.

#'
#'@examples
#'
#'# We first prepare the input X, Z, propensity score
#'
#'attach(dt_Rouse)
#'X = cbind(female,black,bytest,dadeduc,momeduc,fincome)
#'Z = IV
#'propensity = glm(IV~female+black+bytest+dadeduc+momeduc+fincome,
#'                 family=binomial)$fitted.values
#'detach(dt_Rouse)
#'
#'# We next create two list representations of distance matrices
#'# using function create_list_from_scratch. See its focumentation
#'# for more details on using different methods and calipers.
#'
#'# Caveate: please switch the role of treatment and control when
#'# construcitg the second list. Simply let Z = 1 - Z.
#'
#' dist_list_pscore = create_list_from_scratch(Z, X, exact = NULL,
#'                    p = propensity, caliper_low = 0.03, k = 100,
#'                    method = 'maha')
#' matching_output = match_2C_list(Z, dt_Rouse,
#'                  dist_list_pscore)
#'
#' # Please refer to the RMarkdown tutorial for more examples.
#'
#'@return  This function returns the same object as function match_2C_mat.
#'@export

match_2C_list <- function(Z, dataset, dist_list_1, dist_list_2 = NULL,
                          lambda = 1000, controls = 1, overflow = FALSE){

  n_t = sum(Z)
  n_c = length(Z) - n_t
  net1 = treated_control_net(n_t, n_c, dist_list_1, controls)

  if (is.null(dist_list_2)) {
    #cat('Solving the network flow problem', '\n')
    res = solve_network_flow(net1)
    #cat('Finish solving the network flow problem', '\n')
  } else {
    net2 = treated_control_net(n_c, n_t, dist_list_2, controls)
    net = stitch_two_nets(net1, net2, lambda, controls, overflow)
    #cat('Solving the network flow problem', '\n')
    res = solve_network_flow(net)
    #cat('Finish solving the network flow problem', '\n')
  }
  return(construct_outcome(res, dist_list_1, Z, dataset, controls))
}
