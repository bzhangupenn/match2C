#'Optimal Matching with Two Criteria.
#'
#'This function performs an optimal statistical matching that sequentially balances the nominal levels
#'(near-fine balance), the marginal distribution of the propensity score, and the total
#' within-matched-pair Mahalanobis distance.
#'
#'
#'@param Z A length-n vector of treatment indicator.
#'@param X A n-by-p matrix of covariates with column names.
#'@param propensity A vector of estimated propensity score (length(propensity) = length(Z)).
#'@param dataset Dataset to be matched.
#'@param method Method used to compute treated-control distance on the left. The default is the Mahalanobis distance.
#'@param exact A vector of strings indicating which variables need to be exactly matched.
#'@param caliper_left Size of caliper on the left network.
#'@param caliper_right Size of caliper on the right network.
#'@param k_left Connect each treated to k_left controls closest in the propensity score in the left network.
#'@param k_right Connect each treated to k_right controls closest in the propensity score in the right network.
#'@param fb_var A vector giving names of variables in matrix X to be finely balanced.
#'@param controls Number of controls matched to each treated. Default is 1.
#'@param include A binary vector indicating which controls must be included (length(include) = sum(1-Z)).
#'
#'@examples
#'\dontrun{
#'To run the example, one must first install the optmatch package
#'and agree to its terms of use.
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
#'matching_output_double_calipers = match_2C(Z = Z, X = X,
#'propensity = propensity,
#'caliper_left = 0.05, caliper_right = 0.05,
#'k_left = 100, k_right = 100,
#'dataset = dt_Rouse)
#'
#' # Please refer to the vignette for many more examples.
#'}
#'
#'@return  This function returns a list of three objects including the feasibility
#'of the matching problem and the matched controls organized in different formats.
#'See the documentation of the function construct_outcome or the vignette for more
#'details.
#'@export

match_2C <- function(Z, X, propensity,
                     dataset,
                     method = 'maha',
                     exact = NULL,
                     caliper_left = 1,
                     caliper_right = 1,
                     k_left = NULL,
                     k_right = NULL,
                     fb_var = NULL,
                     controls = 1,
                     include = NULL){


  # Construct distance list on the left: Maha subject to pscore caliper
  dist_list_left = create_list_from_scratch(Z = Z, X = X,
                                            exact = exact,
                                            p = propensity,
                                            caliper_low = caliper_left,
                                            k = k_left,
                                            method = method)

  # Construct distance list on the right: L1 distance on the pscore plus fine balance
  dist_list_right = create_list_from_scratch(Z = Z, X = propensity,
                                             p = propensity,
                                             caliper_low = caliper_right,
                                             k = k_right,
                                             method = 'L1')

  if (!is.null(fb_var)) {
    # Construct a fine balance distance list

    dist_list_right_fb = create_list_from_scratch(Z = Z, X = X[, fb_var],
                                                  p = propensity,
                                                  caliper_low = caliper_right,
                                                  k = k_right,
                                                  method = '0/1')
    dist_list_right$d = dist_list_right$d + 1000*dist_list_right_fb$d
  }

  # If we have to include certain controls, add them here
  if (!is.null(include)) {
    dist_list_left = force_control(dist_list_left, Z = Z, include = include)
    dist_list_right = force_control(dist_list_right, Z = Z, include = include)
  }


  if (is.na(dist_list_left)[1] | (is.na(dist_list_right))[1]) {
    cat('Matching is unfeasible. Please increase the caliper size or remove
        the exact matching constraints.')
    return(NA)
  }

  matching_output = match_2C_list(Z = Z, dataset = dataset,
                                  dist_list_1 = dist_list_left,
                                  dist_list_2 = dist_list_right,
                                  lambda = 1000,
                                  controls = controls)

  return(matching_output)
}
