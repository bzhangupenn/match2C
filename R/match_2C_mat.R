#'Perform a pair matching using two user-specified distance matrices.
#'
#'This function performs a pair-matching using two user-specified distance
#'matrices and two calipers. Typically one distance matrix is used to
#'minimize matched-pair differences, and a second distance matrix is used
#'to enforce constraints on marginal distributions of certain variables.
#'
#'This function performs a pair matching via a two-part network.
#'The first part is a network whose treatment-to-control distance
#'matrix is supplied by dist_mat_1. The second part of the
#'network is constructed using distance matrix specified by dist_mat_2. Often,
#'the first part of the network is used to minimize total treated-to-control matched pair
#'distances, and the second part is used to enforce certain marginal constraints.
#'
#'The function constructs two list representations of distance matrices, possibly
#'using the caliper. caliper_1 is applied to p_1 (caliper_2 applied to p_2) in order to
#'construct sparse list representations. For instance, a caliper equal to 0.2 (caliper_1 = 0.2)
#'applied to the propensity score (p_1).
#'
#'lambda is a penalty, or a tuning parameter, that balances these two objectives. When lambda is
#'very large, the network will first minimize the second part of network and then the first part.
#'
#'
#'@param Z A length-n vector of treatment indicator.
#'@param dataset The original dataset.
#'@param dist_mat_1 A user-specified treatment-by-control (n_t-by-n_c) distance matrix.
#'@param dist_mat_2 A second user-specified treatment-by-control (n_t-by-n_c) distance matrix.
#'@param lambda A penalty that controls the trade-off between two parts of the network.
#'@param controls Number of controls matched to each treated.
#'@param p_1 A length-n vector on which caliper_1 applies, e.g. a vector of propensity score.
#'@param caliper_1 Size of caliper_1.
#'@param k_1 Maximum number of controls each treated is connected to in the first network.
#'@param p_2 A length-n vector on which caliper_2 applies, e.g. a vector of propensity score.
#'@param caliper_2 Size of caliper_2.
#'@param k_2 Maximum number of controls each treated is connected to in the second network.
#'@param penalty Penalty for violating the caliper. Set to Inf by default.
#'@param overflow A logical value indicating if overflow protection is turned on.
#'
#'@examples
#'\dontrun{
#'To run the following code, one needs to first install
#'and load the package optmatch.
#'
#'# We first prepare the input X, Z, propensity score
#'
#' #attach(dt_Rouse)
#' #X = cbind(female,black,bytest,dadeduc,momeduc,fincome)
#' #Z = IV
#' #propensity = glm(IV~female+black+bytest+dadeduc+momeduc+fincome,
#' #family=binomial)$fitted.values
#' #n_t = sum(Z)
#' #n_c = length(Z) - n_t
#' #dt_Rouse$propensity = propensity
#' #detach(dt_Rouse)
#'
#'# Next, we use the match_on function in optmatch
#'to create two treated-by-control distance matrices.
#'
#'#library(optmatch)
#'# dist_mat_1 = match_on(IV~female+black+bytest+dadeduc+momeduc+fincome,
#'# method = 'mahalanobis', data = dt_Rouse)
#'
#'# dist_mat_2 = match_on(IV ~ female, method = 'euclidean', data = dt_Rouse)
#'
#'
#' # Feed two distance matrices to the function match_2C_mat without caliper
#' # and a large penalty lambda to enforce (near-)fine balance.
#'
#' #matching_output = match_2C_mat(Z, dt_Rouse, dist_mat_1, dist_mat_2,
#'#                               lambda = 10000, p_1 = NULL, p_2 = NULL)
#'
#' # For more examples, please consult the RMarkdown tutorial.
#'}
#'

#'
#'@return  This function returns a list of three objects including the feasibility
#'of the matching problem and the matched controls organized in different formats.
#'See the documentation of the function construct_outcome or the tutorial for more
#'details.
#'
#'
#'@export

match_2C_mat <- function(Z, dataset, dist_mat_1, dist_mat_2,
                         lambda, controls = 1,
                         p_1 = NULL, caliper_1 = NULL, k_1 = NULL,
                         p_2 = NULL, caliper_2 = NULL, k_2 = NULL,
                         penalty = Inf, overflow = FALSE){

  n_t = sum(Z)
  n_c = length(Z) - n_t

  # construct the list representation of the first distance matrix
  # If--else is not necessary. It is written this way just to be clear.
  if (is.null(p_1))
     dist_list_1 = create_list_from_mat(Z, dist_mat_1, p = NULL)
  else
    dist_list_1 = create_list_from_mat(Z, dist_mat_1, p_1, caliper_1, k_1, penalty)
  if (!is.list(dist_list_1)) return('Hard caliper fails. Please specify a soft caliper.')

  # construct the list representation of the second distance matrix
  # Note: we switch the role of treatment and control

  if (is.null(p_2))
    dist_list_before_reverting = create_list_from_mat(Z, dist_mat_2, p = NULL)
  else
    dist_list_before_reverting = create_list_from_mat(Z, dist_mat_2, p_2, caliper_2, k_2, penalty)
  if (!is.list(dist_list_before_reverting)) return('Hard caliper fails. Please specify a soft caliper.')
  else{
    dist_list_2 = revert_dist_list_cpp(n_t, n_c, dist_list_before_reverting$start_n,
                                       dist_list_before_reverting$end_n,
                                       dist_list_before_reverting$d)
    names(dist_list_2) = c('start_n', 'end_n', 'd')
  }

  cat('Finish converting distance matrices to lists', '\n')

  net1 = treated_control_net(n_t, n_c, dist_list_1, controls)
  net2 = treated_control_net(n_c, n_t, dist_list_2, controls)
  net = stitch_two_nets(net1, net2, lambda, controls, overflow)
  cat('Solving the network flow problem', '\n')
  res = solve_network_flow(net)
  cat('Finish solving the network flow problem', '\n')

  return(construct_outcome(res, dist_list_1, Z, dataset, controls))
}
