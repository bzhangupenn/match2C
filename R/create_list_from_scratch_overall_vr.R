#'Create a sparse list representation of treatment-to-control distance matrix
#'
#'
#'@param Z A length-n vector of treatment indicator.
#'@param X A n-by-p numerical data.frame of covariates, including the covariate(s) used for fine balance. 
#'         The nominal covariate used for fine balance can be kept as a single factor; there is no need to create dummy variables.
#'@param fine The column names of the covariate(s) used for fine balance.
#'@param p A length-n vector on which a caliper applies, e.g. a vector of propensity score.
#'@param caliper_low Size of caliper_inf.
#'@param caliper_high Size of caliper_sup.
#'@param k Connect each treated to the nearest k controls.
#'@param penalty Penalty for violating the caliper. Set to Inf by default.
#'@param ad Adjustment coefficient applied to \eqn{\kappa_{\max}}{kappa_max}.
#'@param dist_func A function used to calculate distance.
#'
#'
#'@return  This function returns a list of three objects: start_n, end_n, and d.
#'
#'@importFrom mvnfast maha
#'@importFrom stats cov var
#'@export
#'
#'

create_list_from_scratch_overall_vr <- function(Z, X, fine = NULL,
                                                p = NULL, caliper_low = NULL, caliper_high = NULL, k = NULL,
                                                penalty = Inf, ad = NULL, dist_func = NULL){
  
  if (is.null(caliper_high)) caliper_high = caliper_low
  n_t = sum(Z)  # number of treated
  n_c = length(Z) - n_t  # number of controls
  
  if (is.vector(X)) X = matrix(X, ncol=1)  # ensure X is a matrix
  X_treated = X[Z == 1,]
  X_control = X[Z == 0,]
  
  if (is.vector(X_treated)) X_treated = matrix(X_treated, ncol=1)
  if (is.vector(X_control)) X_control = matrix(X_control, ncol=1)
  
  p_treated = if (!is.null(p)) p[Z == 1] else NULL
  p_control = if (!is.null(p)) p[Z == 0] else NULL
  
  X_treated_nofine = X_treated[, setdiff(colnames(X_treated), fine), drop = FALSE]
  X_control_nofine = X_control[, setdiff(colnames(X_control), fine), drop = FALSE]
  X_treated_nofine <- as.matrix(X_treated_nofine)
  X_control_nofine <- as.matrix(X_control_nofine)
  # If fine is NULL, X_treated_nofine is exactly X_treated which should be numerical as an input
  # So is X_control_nofine
  
  if (is.null(p)){
    # Fully connected graph without constraints
    start_n = rep(seq(1,n_t), each = n_c)
    end_n = rep(seq(n_t+1, n_t+n_c), n_t)
    d = numeric(n_t*n_c)
    point_start = 1
    point_end = 1
    for (i in 1:n_t){
      temp_d = dist_func(X_control_nofine, X_treated_nofine[i,])
      point_end = point_start + length(temp_d) - 1
      d[point_start:point_end] = temp_d
      point_start = point_end + 1
    }
  }else if (is.infinite(penalty)){
    
    start_n = numeric(n_t*k*1.5)
    end_n = numeric(n_t*k*1.5)
    d = numeric(n_t*k*1.5)
    point_start = 1
    point_end = 1
    for (i in 1:n_t){
      # ind_control_within_caliper = seq(1, n_c, 1)
      ind_control_within_caliper = which((p_control >= p_treated[i] - caliper_low) & (p_control <= p_treated[i] + caliper_high))
      
      #If ind_control_within_caliper = NULL, add three smallest in p_diff
      if (length(ind_control_within_caliper) < 1){
        return('Hard caliper fails. Please specify a soft caliper.')
      }
      
      # Obtain k closest controls if there are still too
      # many after applying the caliper
      if (length(ind_control_within_caliper) > k) {
        p_diff = abs(p_treated[i] - p_control[ind_control_within_caliper])
        p_diff_smallest_k = sort(p_diff)[1:k]
        ind_control_within_caliper = which(p_diff %in% p_diff_smallest_k, arr.ind = TRUE)
        # sorted_indices = order(p_diff)
        # ind_control_within_caliper = ind_control_within_caliper[sorted_indices[1:k]]
      }
      
      point_end = point_start + length(ind_control_within_caliper) - 1
      #cat(i, '\n')
      start_n[point_start:point_end] = rep(i, length(ind_control_within_caliper))
      end_n[point_start:point_end] = n_t + ind_control_within_caliper
      
      controls_within_caliper = as.matrix(X_control_nofine[ind_control_within_caliper,], ncol = dim(X_treated_nofine)[2])
      # Compute distance
      temp_d = dist_func(controls_within_caliper, X_treated_nofine[i,])
      d[point_start:point_end] = temp_d
      point_start = point_end + 1
    }
    start_n = head(start_n, point_end)
    end_n = head(end_n, point_end)
    d = head(d, point_end)
  } else {
    start_n = rep(seq(1,n_t), each = n_c)
    end_n = rep(seq(n_t+1, n_t+n_c), n_t)
    d = numeric(n_t*n_c)
    point_start = 1
    point_end = 1
    for (i in 1:n_t){
      control_outside_caliper = ((p_control < p_treated[i] - caliper_low) | (p_control > p_treated[i] + caliper_high)) + 0
      
      temp_d = dist_func(X_control_nofine, X_treated_nofine[i,]) + control_outside_caliper*penalty
      point_end = point_start + length(temp_d) - 1
      d[point_start:point_end] = temp_d
      point_start = point_end + 1
    }
  }
  
  if (!is.null(fine)) {
    # Extract level indicators for the nominal covariate
    # level_treated = X_treated[, fine]
    # level_control = X_control[, fine]
    # To accommodate two or more covariates
    level_treated <- interaction(
      as.data.frame(X_treated[, fine, drop = FALSE]),
      drop = TRUE
    )
    
    level_control <- interaction(
      as.data.frame(X_control[, fine, drop = FALSE]), drop = TRUE
    )
    
    # levels <- union(levels(level_treated), levels(level_control))   
    levels = sort(unique(c(level_treated, level_control)))
    
    n_b = sapply(levels, function(l) sum(level_treated == l))  # treated counts per level
    N_b = sapply(levels, function(l) sum(level_control == l))  # control counts per level
    
    if (any(n_b == 0)) {
      stop("Some strata have zero treated units. Fine balance is infeasible. ")
    }
    
    if (!all(N_b >= n_b)){
      if (all(n_b >= N_b)){
        stop("Please exchange the position of treated and control group.")
      } else {
        stop("Fine balance is infeasible.")
      }
    }
    
    # We also need to consider near-fine balance here
    
    kappa = min(N_b / n_b)
    
    if (!is.null(ad)) {
      # ad \in (0,1)
      ad_kappa <- ad * kappa
      if (is.finite(ad_kappa) && ad_kappa >= 1) {
        kappa <- ad_kappa
      }
    }
    
    M_b = N_b - floor(kappa * n_b)  # number of nodes to be added per level
    total_added_nodes = sum(M_b)
    
    added_index = n_t + n_c + 1
    
    for (b in seq_along(levels)) {
      n_an <- M_b[b]
      if (n_an <= 0) next
      
      control_level_b <- which(level_control == levels[b])
      for (m in 1:n_an) {
        start_n = c(start_n, rep(added_index, length(control_level_b)))
        end_n   = c(end_n, n_t + control_level_b)
        d       = c(d, rep(0, length(control_level_b)))
        added_index <- added_index + 1
      }
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
