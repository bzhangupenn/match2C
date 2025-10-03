#'Create a sparse list representation of treated-to-control distance matrix
#'
#'
#'@param Z A length-n vector of treatment indicator.
#'@param X A n-by-p numerical data.frame of covariates, including the covariate(s) used for fine balance.
#'@param fine The column names of the covariate(s) used for fine balance.
#'@param p A length-n vector on which a caliper applies, e.g. a vector of propensity score.
#'@param caliper_low Size of caliper_inf.
#'@param caliper_high Size of caliper_sup.
#'@param k Connect each treated to the nearest k controls.
#'@param alpha Tuning parameter.
#'@param penalty Penalty for violating the caliper. Set to Inf by default.
#'@param ad Adjustment coefficient applied to \eqn{\kappa_{\max}}{kappa_max}.
#'@param method Method used to compute treated-control distance.
#'@param dist_func A function used to calculate distance.
#'                 See documentation of "create_list_from_scratch" for more details.
#'
#'
#'@return  This function returns a list of three objects: start_n, end_n, and d.
#'
#'@importFrom mvnfast maha
#'@importFrom stats cov var
#'@export

create_list_from_scratch_vr <- function(Z, X, fine = NULL,
                                        p = NULL, caliper_low = NULL, caliper_high = NULL,
                                        k = NULL, alpha = 1,
                                        penalty = Inf,ad = NULL, method = 'maha', dist_func = NULL){
  
  if (is.null(k)) k = length(Z) - sum(Z)
  
  # Ensure X is matrix/numerical
  # Cast X into matrix if it is a vector
  if (is.vector(X)) X = matrix(X, ncol=1)
  Xco = X[, setdiff(colnames(X), fine), drop = FALSE]
  Xco = as.matrix(Xco)
  
  if (method == 'maha'){
    cov_matrix = chol(stats::cov(Xco))
    # Costomized function computing Maha distance
    compute_maha_dist <- function(X_control, X_treated_i){
      return(mvnfast::maha(X_control, t(as.matrix(X_treated_i)), cov_matrix, isChol=TRUE))
    }
    output = create_list_from_scratch_overall_vr(Z, X, fine, p, caliper_low, caliper_high, k,
                                                 penalty, ad, dist_func = compute_maha_dist)
  }
  
  if (method == 'Hamming'){
    compute_hamming_dist <- function(X_control, X_treated_i){
      return(ncol(X_control) - rowSums(sweep(X_control, 2, as.matrix(X_treated_i)) == 0))
    }
    output = create_list_from_scratch_overall_vr(Z, X, fine, p, caliper_low, caliper_high, k,
                                                 penalty, ad, dist_func = compute_hamming_dist)
  }
  
  if (method == 'L1') {
    compute_L1_dist <- function(X_control, X_treated_i){
      return(rowSums(abs(sweep(X_control, 2, as.matrix(X_treated_i)))))
    }
    output = create_list_from_scratch_overall_vr(Z, X, fine, p, caliper_low, caliper_high, k,
                                                 penalty, ad, dist_func = compute_L1_dist)
  }
  
  if (method == 'L1_convex') {
    compute_L1_convex_dist <- function(X_control, X_treated_i){
      return(alpha*(-rowSums(sweep(X_control, 2, as.matrix(X_treated_i)))) - 0)
    }
    output = create_list_from_scratch_overall_vr(Z, X, fine, p, caliper_low, caliper_high, k,
                                                 penalty, ad, dist_func = compute_L1_convex_dist)
  }
  
  
  if (method == 'vanilla_directional') {
    compute_vanilla_dir_dist <- function(X_control, X_treated_i){
      return(alpha*(-rowSums(sweep(X_control, 2, as.matrix(X_treated_i))) - 0))
    }
    output = create_list_from_scratch_overall_vr(Z, X, fine, p, caliper_low, caliper_high, k,
                                                 penalty, ad, dist_func = compute_vanilla_dir_dist)
  }
  
  if (method == 'hockey_stick') {
    compute_hockey_stick_dist <- function(X_control, X_treated_i){
      d_1 = pmax(-rowSums(sweep(X_control, 2, as.matrix(X_treated_i))), 0)
      #d_2 = pmax(rowSums(sweep(X_control, 2, as.matrix(X_treated_i))), 0)
      return(alpha*(d_1 - 0.01))
    }
    output = create_list_from_scratch_overall_vr(Z, X, fine, p, caliper_low, caliper_high, k,
                                                 penalty, ad, dist_func = compute_hockey_stick_dist)
  }
  
  if (method == '0/1/directional'){
    compute_0_1_dir_dist <- function(X_control, X_treated_i){
      d_1 = (((-rowSums(sweep(X_control, 2, as.matrix(X_treated_i)))) > 0) + 0)
      #d_2 = (((rowSums(sweep(X_control, 2, as.matrix(X_treated_i)))) >= 0) + 0)
      return(alpha*(d_1  - 0.01))
    }
    output = create_list_from_scratch_overall_vr(Z, X, fine, p, caliper_low, caliper_high, k,
                                                 penalty, ad, dist_func = compute_0_1_dir_dist)
  }
  
  if (method == '0/1') {
    compute_01_dist <- function(X_control, X_treated_i){
      return(1 - (rowSums(sweep(X_control, 2, X_treated_i) == 0) == dim(X_control)[2]))
    }
    output = create_list_from_scratch_overall_vr(Z, X, fine, p, caliper_low, caliper_high, k,
                                                 penalty, ad, dist_func = compute_01_dist)
  }
  
  if (method == 'robust maha') {
    if (is.vector(X)) X = matrix(X, ncol=1)
    Xco = X[, setdiff(colnames(X), fine), drop = FALSE]
    Xco = as.matrix(Xco)
    
    # X <- as.matrix(X)
    n<-dim(Xco)[1]
    rownames(Xco) <- 1:n
    for (j in 1:dim(Xco)[2]) Xco[,j]<-rank(Xco[,j])
    cv<-stats::cov(Xco)
    vuntied<-stats::var(1:n)
    rat<-sqrt(vuntied/diag(cv))
    cv<-diag(rat)%*%cv%*%diag(rat)
    cov_matrix = chol(cv) # Cholesky decomp of cov matrix
    
    compute_maha_dist <- function(X_control, X_treated_i){
      return(mvnfast::maha(X_control, t(as.matrix(X_treated_i)), cov_matrix, isChol=TRUE))
    }
    output = create_list_from_scratch_overall_vr(Z, X, fine, p, caliper_low, caliper_high, k,
                                                 penalty,ad,  dist_func = compute_maha_dist)
  }
  
  
  if (method == 'other')
    output = create_list_from_scratch_overall_vr(Z, X, fine, p, caliper_low, caliper_high,
                                                 k, penalty, ad, dist_func = dist_func)
  
  
  if (is.character(output)) {
    cat("Hard caliper fails. Please specify a soft caliper.", '\n')
    return(NA)
  }
  else {
    start_n = output[[1]]
    end_n = output[[2]]
    d = output[[3]]
    
    return(list(start_n = unname(start_n),
                end_n = unname(end_n),
                d = unname(d)))
  }
}
