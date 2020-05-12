#'Create a list representation of a distance matrix.
#'
#'This function creates a ``list representations''
#'of a treatment-by-control distance matrix.
#'
#'This function creates a list representation of a treatment-by-control
#'network. The list representation can be made sparse using a user-specified
#'caliper. A list representation of a treatment-by-control distance matrix
#'consists of the following arguments:
#'\itemize{
#'\item start_n: a vector containing the node numbers of
#'              the start nodes of each arc in the network.
#'\item end_n: a vector containing the node numbers of
#'            the end nodes of each arc in the network.
#'\item d: a vector containing the integer cost of
#'            each arc in the network.
#'}
#'Node 1,2,...,n_t are n_t treatment nodes; n_t + 1, n_t + 2, ..., n_t + n_c
#'are n_c control nodes. start_n, end_n, and d should have the same lengths,
#'all of which equal to the number of edges.
#'
#'There are two options for users to make a network sparse. First, caliper
#'is a value applied to the vector p to avoid connecting treated to controls
#'whose covariate/propensity score defined by p is outside p +/- caliper.
#'Second, within a specified caliper, sometimes there are too many controls
#'connected to the treated, and we can further trim down this number up to k
#'with resctricting attention to the k nearest (in p) to each treated.
#'
#'@param Z A length (n = n_t + n_c) vector of treatment indicators.
#'@param dist_mat A treatment-by-control (n_t-by-n_c) distance matrix.
#'@param p A vector of length (n_t + n_c) on which caliper applies (e.g. propensity scores)
#'@param caliper Size of the caliper.
#'@param k Connect each treated to the nearest k controls
#'@param penalty Penalty for violating the caliper. Set to Inf by default.
#'
#'
#'@return  This function returns a list that consists of three arguments: start_n, end_n, and d,
#'         as described above.
#'
#'@examples
#'\dontrun{
#'#To run the following code, one needs to first install
#'#and load the package optmatch.
#'
#'# We first prepare the input X, Z, propensity score
#'
#' attach(dt_Rouse)
#' X = cbind(female,black,bytest,dadeduc,momeduc,fincome)
#' Z = IV
#' propensity = glm(IV~female+black+bytest+dadeduc+momeduc+fincome,
#' family=binomial)$fitted.values
#' n_t = sum(Z)
#' n_c = length(Z) - n_t
#' dt_Rouse$propensity = propensity
#' detach(dt_Rouse)
#'
#'# Next, we use the match_on function in optmatch
#'# to create two treated-by-control distance matrices.
#'
#'library(optmatch)
#' dist_mat_1 = match_on(IV~female+black+bytest+dadeduc+momeduc+fincome,
#' method = 'mahalanobis', data = dt_Rouse)
#'
#' # Convert the distance matrix to a distance list
#' dist_list_1 = create_list_from_mat(Z, dist_mat_1, p = NULL)
#'
#' # For more examples, please consult the RMarkdown tutorial.
#'}
#'
#'@export

create_list_from_mat <- function(Z, dist_mat, p = NULL, caliper = NULL, k = NULL, penalty = Inf){
  'This function creates a list representation of a distance matrix.
   Z: A vector of length n (n = n_t + n_c) of treatment indicators.
   dist_mat: A treatment-by-control (n_t-by-n_c) distance matrix.
   p: A vector of length (n_t + n_c) on which caliper applies (e.g. propensity scores).
   caliper: Size of the caliper.
  '
  n_t = nrow(dist_mat)
  n_c = ncol(dist_mat)

  # Do not use caliper
  if (is.null(p)){
    start_n = rep(seq(1,n_t,1), each = n_c)
    end_n = rep(seq(n_t+1, n_t+n_c, 1), n_t)
    d = c(t(dist_mat))
  } else if (is.infinite(penalty)){
    if (is.null(k)) k = n_c

    start_n = numeric(n_t*k*1.5)
    end_n = numeric(n_t*k*1.5)
    d = numeric(n_t*k*1.5)

    p_treated = p[which(Z == 1)]
    p_control = p[which(Z == 0)]

    point_start = 1
    point_end = 1
    for (i in 1:n_t){
      p_diff = abs(p_treated[i] - p_control)

      # All those controls within the caliper
      ind_control_within_caliper = which(p_diff <= caliper)

      #If ind_control_within_caliper = NULL, add three smallest in p_diff
      if (length(ind_control_within_caliper) < 1){
        message('Hard caliper fails. Please specify a soft caliper.', '\n')
        return(NA)
      }

      # Obtain k closest controls if there are still too
      # many after applying the caliper

      if (length(ind_control_within_caliper) > k) {
        p_diff_smallest_k = sort(p_diff)[1:k]
        ind_control_within_caliper = which(p_diff %in% p_diff_smallest_k, arr.ind = TRUE)
      }

      point_end = point_start + length(ind_control_within_caliper) - 1
      start_n[point_start:point_end] = rep(i, length(ind_control_within_caliper))
      end_n[point_start:point_end] = n_t + ind_control_within_caliper

      # Compute Mahalanobis distance
      temp_d = dist_mat[i, ind_control_within_caliper]
      d[point_start:point_end] = temp_d
      point_start = point_end + 1
    }
    start_n = head(start_n, point_end)
    end_n = head(end_n, point_end)
    d = head(d, point_end)
  } else{
    # Soft caliper
    if (is.null(k)) k = n_c
    d = numeric(n_t*k*1.5)

    p_treated = p[which(Z == 1)]
    p_control = p[which(Z == 0)]

    point_start = 1
    point_end = 1
    for (i in 1:n_t){
      p_diff = abs(p_treated[i] - p_control)

      # All those controls outside caliper
      control_outside_caliper = (p_diff > caliper) + 0

      point_end = point_start + n_c - 1

      temp_d = dist_mat[i, ] + control_outside_caliper*penalty
      d[point_start:point_end] = temp_d
      point_start = point_end + 1
    }
    start_n = rep(seq(1,n_t,1), each = n_c)
    end_n = rep(seq(n_t+1, n_t+n_c, 1), n_t)
    d = head(d, point_end)
  }
  return(list(start_n = unname(start_n), end_n = unname(end_n), d = unname(d)))
}
