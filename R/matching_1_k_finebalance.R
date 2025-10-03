#' create_list_from_scratch_fm is the variant of create_list_from_scratch.
#' It is designed for 1-to-k matching with fine balance.
#' 
#'@param Z A length-n vector of treatment indicator.
#'@param X A n-by-p numerical data.frame of covariates, including the covariate(s) used for fine balance.
#'@param fine The column names of the covariate(s) used for fine balance.
#'@param k The number of control units matched to each treated unit.
#'@param alpha Tuning parameter.
#'@param method Method used to compute treated-control distance.
#'@param dist_func A function used to calculate distance.
#'                 See documentation of "create_list_from_scratch" for more details.
#'
#'@return  This function returns a list of three objects: start_n, end_n, and d.
#'
create_list_from_scratch_fm <- function(Z, X, fine, k, alpha = 1, method = 'maha', dist_func = NULL){
  if (is.vector(X)) X = matrix(X, ncol=1)
  
  Xco = X[, setdiff(colnames(X), fine), drop = FALSE]
  Xco = as.matrix(Xco)
  
  if (method == 'maha'){
    cov_matrix = chol(stats::cov(Xco))
    # Costomized function computing Maha distance
    compute_maha_dist <- function(X_control, X_treated_i){
      return(mvnfast::maha(X_control, t(as.matrix(X_treated_i)), cov_matrix, isChol=TRUE))
    }
    output = create_list_from_scratch_overall_fm(Z, X, fine, k, dist_func = compute_maha_dist)
  }
  
  if (method == 'Hamming'){
    compute_hamming_dist <- function(X_control, X_treated_i){
      return(ncol(X_control) - rowSums(sweep(X_control, 2, as.matrix(X_treated_i)) == 0))
    }
    output = create_list_from_scratch_overall_fm(Z, X, fine, k, dist_func = compute_hamming_dist)
  }
  
  if (method == 'L1') {
    compute_L1_dist <- function(X_control, X_treated_i){
      return(rowSums(abs(sweep(X_control, 2, as.matrix(X_treated_i)))))
    }
    output = create_list_from_scratch_overall_fm(Z, X, fine, k, dist_func = compute_L1_dist)
  }
  
  if (method == 'L1_convex') {
    compute_L1_convex_dist <- function(X_control, X_treated_i){
      return(alpha*(-rowSums(sweep(X_control, 2, as.matrix(X_treated_i)))) - 0)
    }
    output = create_list_from_scratch_overall_fm(Z, X, fine, k,  dist_func = compute_L1_convex_dist)
  }
  
  
  if (method == 'vanilla_directional') {
    compute_vanilla_dir_dist <- function(X_control, X_treated_i){
      return(alpha*(-rowSums(sweep(X_control, 2, as.matrix(X_treated_i))) - 0))
    }
    output = create_list_from_scratch_overall_fm(Z, X, fine, k,  dist_func = compute_vanilla_dir_dist)
  }
  
  if (method == 'hockey_stick') {
    compute_hockey_stick_dist <- function(X_control, X_treated_i){
      d_1 = pmax(-rowSums(sweep(X_control, 2, as.matrix(X_treated_i))), 0)
      #d_2 = pmax(rowSums(sweep(X_control, 2, as.matrix(X_treated_i))), 0)
      return(alpha*(d_1 - 0.01))
    }
    output = create_list_from_scratch_overall_fm(Z, X, fine, k,  dist_func = compute_hockey_stick_dist)
  }
  
  if (method == '0/1/directional'){
    compute_0_1_dir_dist <- function(X_control, X_treated_i){
      d_1 = (((-rowSums(sweep(X_control, 2, as.matrix(X_treated_i)))) > 0) + 0)
      #d_2 = (((rowSums(sweep(X_control, 2, as.matrix(X_treated_i)))) >= 0) + 0)
      return(alpha*(d_1  - 0.01))
    }
    output = create_list_from_scratch_overall_fm(Z, X, fine, k, dist_func = compute_0_1_dir_dist)
  }
  
  if (method == '0/1') {
    compute_01_dist <- function(X_control, X_treated_i){
      return(1 - (rowSums(sweep(X_control, 2, X_treated_i) == 0) == dim(X_control)[2]))
    }
    output = create_list_from_scratch_overall_fm(Z, X, fine, k, dist_func = compute_01_dist)
  }
  
  if (method == 'robust maha') {
    # Need further refinement
    if (is.vector(X)) X = matrix(X, ncol=1)
    
    X <- as.matrix(X)
    n<-dim(X)[1]
    rownames(X) <- 1:n
    for (j in 1:dim(X)[2]) X[,j]<-rank(X[,j])
    cv<-stats::cov(X)
    vuntied<-stats::var(1:n)
    rat<-sqrt(vuntied/diag(cv))
    cv<-diag(rat)%*%cv%*%diag(rat)
    cov_matrix = chol(cv) # Cholesky decomp of cov matrix
    
    compute_maha_dist <- function(X_control, X_treated_i){
      return(mvnfast::maha(X_control, t(as.matrix(X_treated_i)), cov_matrix, isChol=TRUE))
    }
    output = create_list_from_scratch_overall_fm(Z, X, fine, k, dist_func = compute_maha_dist)
  }
  
  
  if (method == 'other')
    output = create_list_from_scratch_overall_fm(Z, X, fine, k,  dist_func = dist_func)
  
  
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

#' create_list_from_scratch_overall_fm is the variant of create_list_from_overall_scratch.
#' It is designed for 1-to-k matching with fine balance.
#' 
#'@param Z A length-n vector of treatment indicator.
#'@param X A n-by-p numerical data.frame of covariates, including the covariate(s) used for fine balance.
#'@param fine The column names of the covariate(s) used for fine balance.
#'@param k The number of control units matched to each treated unit.
#'@param dist_func A function used to calculate distance.
#'                 See documentation of "create_list_from_scratch" for more details.
#'
#'@return  This function returns a list of three objects: start_n, end_n, and d.
#'
create_list_from_scratch_overall_fm <- function(Z, X, fine, k, dist_func){
  if (!is.numeric(k) || length(k) != 1 || k != as.integer(k)) {
    stop("Parameter 'k' must be a single integer value.")
  }
  
  if (is.null(fine)){
    stop("Please select covariates for fine balance")
  }
  
  if (!all(fine %in% colnames(X))) {
    stop("Some 'fine' variables are not found in the data.")
  }
  
  if (is.vector(X)) X = matrix(X, ncol=1)
  X_treated = X[Z == 1,]
  X_control = X[Z == 0,]
  
  if (is.vector(X_treated)) X_treated = matrix(X_treated, ncol=1)
  if (is.vector(X_control)) X_control = matrix(X_control, ncol=1)
  
  level_treated <- interaction(
    as.data.frame(X_treated[, fine, drop = FALSE]),
    drop = TRUE
  )
  
  level_control <- interaction(
    as.data.frame(X_control[, fine, drop = FALSE]), drop = TRUE
  )
  
  levels = sort(unique(c(level_treated, level_control)))
  
  treated_counts = sapply(levels, function(l) sum(as.character(level_treated) == l))
  control_counts = sapply(levels, function(l) sum(as.character(level_control) == l))
  
  X_treated_nofine = X_treated[, setdiff(colnames(X_treated), fine), drop = FALSE]
  X_control_nofine = X_control[, setdiff(colnames(X_control), fine), drop = FALSE]
  X_treated_nofine <- as.matrix(X_treated_nofine)
  X_control_nofine <- as.matrix(X_control_nofine)
  
  comparison <- control_counts >= k * treated_counts
  if (any(!comparison)) {
    stop("Please try a smaller value of 'k', or ensure that there are more units 
            in the control group than in the treated group for each level of the covariate.")
  }
  
  n_t <- sum(Z)
  n_c <- length(Z) - sum(Z)
  
  start_n <- rep(seq(1,n_t), each = n_c)
  end_n <- rep(seq(n_t+1, n_t+n_c), n_t)
  d <- numeric(n_t*n_c)
  point_start <- 1
  point_end <- 1
  for (i in 1:n_t){
    temp_d <- dist_func(X_control_nofine, X_treated_nofine[i,])
    point_end <- point_start + length(temp_d) - 1
    d[point_start:point_end] <- temp_d
    point_start <- point_end + 1
  }
  
  added_index <- n_t + n_c + 1
  difference <- control_counts - k*treated_counts 
  for (lvl in seq_along(levels)) {
    n_an <- difference[lvl]
    if (n_an <= 0) next
    
    control_lvl <- which(level_control == levels[lvl])
    for (m in 1:n_an) {
      start_n <- c(start_n, rep(added_index, length(control_lvl)))
      end_n <- c(end_n, n_t + control_lvl)
      d <- c(d, rep(0, length(control_lvl)))
      added_index <- added_index + 1
    }
  }
  
  if (any(d < 0)){
    offset_d <- -min(d)
    d <- d + offset_d
  }
  
  return(list(start_n = unname(start_n),
              end_n = unname(end_n),
              d = unname(d)))
  
}

#' treated_control_net_fm is the variant of treated_control_net.
#' It is designed for 1-to-k matching with fine balance.
#' 
#'@param n_t Number of treated units.
#'@param n_c Number of control units.
#'@param dist_list A list representation of the distance matrix.
#'@param k Number of control units matched to each treated unit.
treated_control_net_fm <- function(n_t, n_c, dist_list, k){
  start_n = dist_list$start_n
  end_n = dist_list$end_n
  d = dist_list$d
  # lengths of start_n, end_n, and d need to be the same.
  stopifnot(length(start_n) == length(end_n) && length(end_n) == length(d))
  
  # Create edges from each control to the sink, from the source node
  # to each treated and each added
  ta <- unique(start_n)
  source_uf <- length(ta)
  sink_index <- max(c(start_n, end_n)) + 2
  start_n = c(rep(1, source_uf), start_n + 1, seq(n_t + 2, n_t + n_c + 1, 1))
  end_n = c(ta + 1, end_n + 1, rep(sink_index, n_c))
  d = c(rep(0, source_uf), d, rep(0, n_c))
  num_edge = length(d)
  cap = c(rep(k, n_t), rep(1, source_uf - n_t), 
          rep(1, num_edge - source_uf - n_c), 
          rep(1, n_c))
  b = c(k*n_t + source_uf - n_t, rep(0, n_t + n_c), rep(0, source_uf - n_t), - n_c)
  
  net = list(startn = start_n,
             endn = end_n,
             ucap = cap,
             cost = d,
             b =  b)
}
