#'Create a sparse list representation of treated-to-control distance
#'matrix with a fixed number caliper with L1-distance.
#'
#'This function takes in a n-by-p matrix of observed covariates,
#'a length-n vector of treatment indicator, a caliper, and construct
#'a possibly sparse list representation of the distance matrix with
#'Mahalanobis distance. Note that this function is of limited interest
#'to most users.
#'
#'
#'
#'@param Z A length-n vector of treatment indicator.
#'@param X A n-by-p matrix of covariates.
#'@param exact A vector of strings indicating which variables are to be exactly matched.
#'@param soft_exact If set to TRUE, the exact constraint is enforced up to a large penalty.
#'@param p A length-n vector on which a caliper applies, e.g. a vector of propensity score.
#'@param caliper_low Size of caliper_inf.
#'@param caliper_high Size of caliper_sup.
#'@param k Connect each treated to the nearest k controls
#'@param penalty Penalty for violating the caliper. Set to Inf by default.
#'@param dist_func A function used to calculate distance
#'
#'
#'@return  This function returns a list of three objects: start_n, end_n, and d.
#'         See documentation of function ``create_list_from_mat'' for more details.
#'
#'@examples
#'# We first prepare the input X, Z, propensity score
#'
#'attach(dt_Rouse)
#'X = cbind(female,black,bytest,dadeduc,momeduc,fincome)
#'Z = IV
#'propensity = glm(IV~female+black+bytest+dadeduc+momeduc+fincome,
#'                 family=binomial)$fitted.values
#'detach(dt_Rouse)
#'
#'# Define a Mahalanobis-distance function
#'
#' cov_matrix = chol(cov(X))
#' compute_maha_dist <- function(X_control, X_treated_i){
#'  return(mvnfast::maha(X_control, t(as.matrix(X_treated_i)), cov_matrix, isChol=TRUE))
#'}
#'
#'# create a distance list using distance function compute_maha_dist
#'output = create_list_from_scratch_overall(Z, X, p = propensity,
#'                                caliper_low = 0.05, k = 100,
#'                                dist_func = compute_maha_dist)
#'
#'
#'# More examples, including how to use a user-supplied
#'# distance function, can be found in the accompanying RMarkdown tutorial.
#'
#'@importFrom mvnfast maha
#'@importFrom stats cov var
#'@export

create_list_from_scratch_overall <- function(Z, X, exact = NULL, soft_exact = FALSE,
                                          p = NULL, caliper_low = NULL, caliper_high = NULL, k = NULL,
                                          penalty = Inf, dist_func = NULL){

  if (is.null(caliper_high)) caliper_high = caliper_low
  n_t = sum(Z) # n_t is number of treated
  n_c = length(Z) - n_t
  # Cast X into matrix if it is a vector
  if (is.vector(X)) X = matrix(X, ncol=1)

  X_treated = X[Z == 1,]
  X_control = X[Z == 0,]

  if (is.vector(X_treated)) X_treated = matrix(X_treated, ncol=1)
  if (is.vector(X_control)) X_control = matrix(X_control, ncol=1)

  # No exact matching requirement or pscore
  # Create a fully-connected graph.
  if (is.null(exact) & is.null(p)){
    start_n = rep(seq(1,n_t,1), each = n_c)
    end_n = rep(seq(n_t+1, n_t+n_c, 1), n_t)
    d = numeric(n_t*n_c)
    point_start = 1
    point_end = 1
    for (i in 1:n_t){
      temp_d = dist_func(X_control, X_treated[i,])
      point_end = point_start + length(temp_d) - 1
      d[point_start:point_end] = temp_d
      point_start = point_end + 1
    }
  }

  # Exact matching constraint and
  # ``hard'' caliper constraint apply
  else if ((is.infinite(penalty)) & (!soft_exact)){

    start_n = numeric(n_t*k*1.5)
    end_n = numeric(n_t*k*1.5)
    d = numeric(n_t*k*1.5)
    point_start = 1
    point_end = 1
    for (i in 1:n_t){

      ind_control_exact = seq(1, n_c, 1)
      ind_control_within_caliper = seq(1, n_c, 1)

      if (!is.null(exact)){
        # All those controls with the same ``exact'' variable
        #dt_logic = X_control[,exact] == X_treated[i, exact]
        #if (is.vector(dt_logic)) dt_logic = matrix(dt_logic, ncol=1)
        #ind_control_exact = which(apply(dt_logic, 1, all))
        X_control_exact_cols = X_control[,exact]
        if (is.vector(X_control_exact_cols)) X_control_exact_cols = data.frame(X_control_exact_cols)
        ind_control_exact = which(apply(X_control_exact_cols, 1,
                                        function(x) identical(unname(x), unname(X_treated[i, exact]))))

      }

      if (!is.null(p)){

        p_treated = p[which(Z == 1)]
        p_control = p[which(Z == 0)]

        # All those controls within the caliper
        ind_control_within_caliper = which((p_control >= p_treated[i] - caliper_low) & (p_control <= p_treated[i] + caliper_high))

      }

      # Find those who are exactly matched on ``exact''
      # and within caliper
      ind_control = intersect(ind_control_within_caliper, ind_control_exact)
      #cat(i, length(ind_control_within_caliper), '\n')

      #If ind_control_within_caliper = NULL, add three smallest in p_diff
      if (length(ind_control) < 1){
        return('Hard caliper fails. Please specify a soft caliper.')
      }



      # Obtain k closest controls if there are still too
      # many after applying the caliper

      if (length(ind_control) > k) {
        p_diff = abs(p_treated[i] - p_control)
        p_diff_smallest_k = sort(p_diff)[1:k]
        ind_control = which(p_diff %in% p_diff_smallest_k, arr.ind = TRUE)
      }

      point_end = point_start + length(ind_control) - 1
      #cat(i, '\n')
      start_n[point_start:point_end] = rep(i, length(ind_control))
      end_n[point_start:point_end] = n_t + ind_control

      controls_within_caliper = as.matrix(X_control[ind_control,], ncol = dim(X_treated)[2])
      # Compute distance
      temp_d = dist_func(controls_within_caliper, X_treated[i,])
      d[point_start:point_end] = temp_d
      point_start = point_end + 1
    }
    start_n = head(start_n, point_end)
    end_n = head(end_n, point_end)
    d = head(d, point_end)
  }

  # Exact matching constraint applies,
  # and soft caliper constraint applies.
  else if ((!is.infinite(penalty)) & (!soft_exact)){
    start_n = numeric(n_t*k*1.5)
    end_n = numeric(n_t*k*1.5)
    d = numeric(n_t*k*1.5)
    point_start = 1
    point_end = 1
    for (i in 1:n_t){

      ind_control_exact = seq(1, n_c, 1)
      if (!is.null(exact)){
        # All those controls with the same ``exact'' variable
        #dt_logic = X_control[,exact] == X_treated[i, exact]
        #if (is.vector(dt_logic)) dt_logic = matrix(dt_logic, ncol=1)
        #ind_control_exact = which(apply(dt_logic, 1, all))
        X_control_exact_cols = X_control[,exact]
        if (is.vector(X_control_exact_cols)) X_control_exact_cols = data.frame(X_control_exact_cols)
        ind_control_exact = which(apply(X_control_exact_cols, 1,
                                        function(x) identical(unname(x), unname(X_treated[i, exact]))))

      }

      point_end = point_start + length(ind_control_exact) - 1
      #cat(i, '\n')
      start_n[point_start:point_end] = rep(i, length(ind_control_exact))
      end_n[point_start:point_end] = n_t + ind_control_exact

      controls_exact_match = as.matrix(X_control[ind_control_exact,], ncol = dim(X_treated)[2])
      # Compute Mahalanobis distance
      p_treated = p[which(Z == 1)]
      p_control = p[which(Z == 0)]
      p_control_exact = p_control[ind_control_exact]

      # Whether or not each control is within the caliper
      control_outside_caliper = ((p_control_exact < p_treated[i] - caliper_low) | (p_control_exact > p_treated[i] + caliper_high)) + 0

      temp_d = dist_func(controls_exact_match, X_treated[i,]) + control_outside_caliper*penalty
      d[point_start:point_end] = temp_d
      point_start = point_end + 1
    }
    start_n = head(start_n, point_end)
    end_n = head(end_n, point_end)
    d = head(d, point_end)
  }
  # soft exact constraints plus soft caliper
  else if (soft_exact){
    start_n = rep(seq(1,n_t,1), each = n_c)
    end_n = rep(seq(n_t+1, n_t+n_c, 1), n_t)
    d = numeric(n_t*n_c)
    point_start = 1
    point_end = 1
    for (i in 1:n_t){
      temp_d = dist_func(X_control, X_treated[i,])

      # find ind_exact
      X_control_exact_cols = X_control[,exact]
      if (is.vector(X_control_exact_cols)) X_control_exact_cols = data.frame(X_control_exact_cols)
      ind_control_exact = which(apply(X_control_exact_cols, 1,
                                      function(x) identical(unname(x), unname(X_treated[i, exact]))))
      ind_control_not_exact = setdiff(seq(1, n_c, 1), ind_control_exact)
      temp_d[ind_control_not_exact] = temp_d[ind_control_not_exact] + 1000

      if (!is.null(p)){ # only do sof caliper for now
        # find within/outside the caliper
        p_treated = p[which(Z == 1)]
        p_control = p[which(Z == 0)]
        control_outside_caliper = ((p_control < p_treated[i] - caliper_low) | (p_control > p_treated[i] + caliper_high)) + 0
        temp_d = temp_d + control_outside_caliper*penalty
      }


      point_end = point_start + length(temp_d) - 1
      d[point_start:point_end] = temp_d
      point_start = point_end + 1
    }
  }

  if (any(d < 0)){
    offset_d = -min(d)
    d = d + offset_d
  }

  return(list(start_n = unname(start_n),
              end_n = unname(end_n),
              d = unname(d)))
}
