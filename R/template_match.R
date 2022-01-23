#'Optimal Matching with Two Criteria.
#'
#'This function takes as arguments a dataset to be matched and a template, and
#'outputs matched pairs that are closely matched, well balanced, and mimicking
#'the user-supplied template in covariates' distributions of the given template.
#'
#'Please refer to the vignette for reproducible examples.
#'
#'@param template A dataframe of template units.
#'@param X A n-by-p matrix of covariates with column names.
#'@param Z A length-n vector of treatment indicator.
#'@param dataset Dataset to be matched.
#'@param multiple Number of treated units matched to each template unit. Default is 1.
#'@param lambda A tuning parameter controlling the trade-off between internal and external validity. A large lambda favors resemblance to the template.
#'@param caliper_gscore Size of generalizability caliper.
#'@param k_gscore Connect each template unit to k_gscore treated units closest in the generalizability score.
#'@param penalty_gscore Penalty for violating the generalizability caliper. Set to Inf by default.
#'@param caliper_pscore Size of propensity score caliper.
#'@param k_pscore Connect each treated to k_pscore control units closest in the propensity score.
#'@param penalty_pscore Penalty for violating the propensity score caliper. Set to Inf by default.
#'
#'
#'
#'@return  This function returns a list of three objects: 1) feasible: 0/1 depending on the
#'feasibility of the matching problem; 2) match_treated: a data frame of the matched treated
#'units; 3) match_control: a data frame of the matched control units.
#'
#'@importFrom stats glm
#'@export

template_match <- function(template, X, Z,
                           dataset, multiple = 1, lambda = 1,
                           caliper_gscore = 1,
                           k_gscore = NULL,
                           penalty_gscore = Inf,
                           caliper_pscore = 1,
                           k_pscore = NULL,
                           penalty_pscore = Inf){



  n_template = dim(template)[1] # Number of units in the template
  d = dim(template)[2] # Number of covariates in the template
  n_t = sum(Z) # Number of treated units in the OBS dataset
  n_c = sum(1-Z) # Number of control units in the OBS dataset


  # Cast X into matrix if it is a vector
  if (is.vector(X)) X = matrix(X, ncol = 1)

  X_treated = X[Z == 1,]
  X_control = X[Z == 0,]


  # Below we construct the left-part of the network,
  # where units in the template are connected to the
  # treated units in the OBS dataset

  # Estimate the generalizability score
  Z_left = c(rep(1, n_template), rep(0, n_t)) # template membership
  X_left = rbind(template, X_treated[, 1:d])
  X_left = as.matrix(X_left)
  gscore = glm(Z_left ~ X_left, family = 'binomial')$fitted.values


  # Construct a distance list for the left-part
  dist_list_gscore_maha = create_list_from_scratch(Z_left, X_left, exact = NULL,
                                                   p = gscore,
                                                   caliper_low = caliper_gscore,
                                                   k = k_gscore,
                                                   method = 'robust maha',
                                                   penalty = penalty_gscore)

  net_left = treated_control_net(n_t = n_template, n_c = n_t,
                                 dist_list = dist_list_gscore_maha,
                                 controls = multiple)


  # Next, we construct the right-part of the network,
  # where controls units in the OBS dataset are connected
  # to the treated units in the OBS dataset

  # Estimate the propensity score
  Z_right = c(rep(1, n_t), rep(0, n_c)) # treatment membership
  X_right = rbind(X_treated, X_control)
  X_right = as.matrix(X_right)
  pscore = glm(Z_right ~ X_right, family = 'binomial')$fitted.values


  # Construct a distance list for the right part
  dist_list_pscore_maha = create_list_from_scratch(Z_right, X_right, exact = NULL,
                                                   p = pscore,
                                                   caliper_low = caliper_pscore,
                                                   k = k_pscore,
                                                   method = 'robust maha',
                                                   penalty = penalty_pscore)


  net_right = treated_control_net(n_t = n_t, n_c = n_c,
                                  dist_list = dist_list_pscore_maha,
                                  controls = 1)

  two_net = stitch_two_nets_template(net_left, net_right, n_c,
                                     lambda = lambda, multiple = multiple)
  res = solve_network_flow(two_net)



  num_edges_left = length(net_left$startn)
  return(construct_outcome_template(res, num_edges_left, Z, dataset))

}
