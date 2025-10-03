#'Perform 1-to-k matching with fine balance.
#'
#'
#'@param Z A length-n vector of treatment indicator.
#'@param X A n-by-p numerical data.frame of covariates, including the covariate(s) used for fine balance.
#'@param fine The column names of the covariate(s) used for fine balance.
#'@param k The number of control units matched to each treated unit.
#'       (Or the number of treated units matched to each control unit if there are more treated units.)
#'@param data The original dataset used to construct final results.
#'@param method Method used to compute treated-control distance.
#'
#'@return  This function returns a data frame that is the same as the input "data", 
#'except that a column called "distance" and a column called "matched_to_index" are added to it. 
#'Or NULL if such matching is infeasible.
#'
#'@examples
#'\dontrun{
#'data("dt_rhc")
#'data("dt_rhc_ac")
#'p <- glm(z ~ . - ninsclasMedicare 
#'               - `ninsclasMedicare & Medicaid` 
#'               - `ninsclasNo insurance` 
#'               - `ninsclasPrivate` 
#'               - `ninsclasPrivate & Medicare`,
#'               data = dt_rhc_ac, family = binomial())
#'dt_rhc_ac$p <-  predict(p, type = "response")
#'Group_index <- assign_int(dt_rhc_ac$p,beta = 4)
#'dt_rhc_ac$subgroup <- Group_index
#'
#'dt_rhc$p <- dt_rhc_ac$p
#'dt_rhc$subgroup <- dt_rhc_ac$subgroup
#'dt_rhc_1 <- dt_rhc[dt_rhc$subgroup==1,]
#'dt_rhc_ac_1 <- dt_rhc_ac[dt_rhc_ac$subgroup == 1,] 
#'dt_rhc_ac_1_cov <- dt_rhc_ac_1[, setdiff(colnames(dt_rhc_ac_1), c("z","subgroup")), drop = FALSE]
#'
#'A_matched_1 <- matching_fm(Z = dt_rhc_ac_1$z, X = dt_rhc_ac_1_cov, 
#'                           fine = c("ninsclasMedicare",
#'                                    "ninsclasMedicare & Medicaid",
#'                                    "ninsclasNo insurance",
#'                                    "ninsclasPrivate",
#'                                    "ninsclasPrivate & Medicare"),
#'                           data = dt_rhc_1, k = 1, method = "maha")
#'}
#'
#'
#'@importFrom mvnfast maha
#'@importFrom stats cov var
#'@export

matching_fm <- function(Z,X,fine,data,k,method = "maha"){
  X <- mahapd(X, fine)
  out <- create_list_from_scratch_fm(Z = Z, X = X, fine = fine, k = k, method = method)
  out$d <- 10000*out$d
  net <- treated_control_net_fm(n_t = sum(Z), n_c = length(Z) - sum(Z), dist_list = out, k = k)
  sol <- solve_network_flow(net)
  sim <- construct_outcome_vr(res = sol, dist_list = out, Z = Z, dataset = data)
  A <- sim$data_without_unmatched_controls
  A <- A[, setdiff(colnames(A), "matched_to"), drop = FALSE]
  return(A)
}



#' Perform a pair matching with near fine balance.
#' Arguments are identical to those of matching_fm except that "lambda", a penalty 
#' that does a trade-off between two parts of the network, is added and "k" is removed.
#'@param Z A length-n vector of treatment indicator.
#'@param X A n-by-p numerical data.frame of covariates, including the covariate(s) used for fine balance.
#'@param fine The column names of the covariate(s) used for fine balance.
#'@param data The original dataset used to construct final results.
#'@param method Method used to compute treated-control distance.
#'@param lambda A penalty that does a trade-off between two parts of the network.
#'
#'@return This function returns a data frame that is the same as the input "data", 
#'except that a column called "distance" and a column called "matched_to_index" are added to it. 
#'
matching_nf <-  function(Z,X,fine,data,method = "maha", lambda = 1000){
  X <- mahapd(X)
  dist_list_1 <- create_list_from_scratch(Z = Z, 
                                          X = as.matrix(X[, !(colnames(X) %in% fine)]), method = "maha")
  dist_list_2 <- create_list_from_scratch(Z = Z, 
                                          X = as.matrix(X[, fine]), method = "0/1")
  matching_out <- match_2C_list(Z = Z, dataset = data, dist_list_1 = dist_list_1,
                                dist_list_2 = dist_list_2, lambda = lambda)
  A <- matching_out$matched_data_in_order
  A <- A[!is.na(A$matched_set), , drop = FALSE]
  A$matched_to_index <- NA_integer_
  for (ms in unique(A$matched_set)) {
    rows <- which(A$matched_set == ms)
    control_row <- rows[A$Z[rows] == 0]
    treated_row <- rows[A$Z[rows] == 1]
    
    if (length(control_row) == 1 && length(treated_row) == 1) {
      A$matched_to_index[treated_row] <- as.integer(rownames(A)[control_row])
    }
  }
  
  A <- A[, setdiff(colnames(A), "matched_set"), drop = FALSE]
}



#'Performs 1-to-k matching with fine balance, or pair matching with near-fine balance 
#'when fine balance is infeasible, depending on the data.
#'
#'
#'@param Z A length-n vector of treatment indicator.
#'@param X A n-by-p numerical data.frame of covariates, including the covariate(s) used for fine balance.
#'@param fine The column names of the covariate(s) used for fine balance.
#'@param data The original dataset used to construct final results.
#'@param kmax The maximum number of control units matched to each treated unit. 
#'            Its value is identical with the number of strata or subgroups
#'@param method Method used to compute treated-control distance.
#'@param lambda A penalty that does a trade-off between two parts of the network.
#'
#'@return  This function returns a data frame that is the same as the input "data", 
#'except that a column called "distance" and a column called "matched_to_index" are added to it. 
#'
#'@examples
#'\dontrun{
#'data("dt_rhc")
#'data("dt_rhc_ac")
#'p <- glm(z ~ . - ninsclasMedicare 
#'               - `ninsclasMedicare & Medicaid` 
#'               - `ninsclasNo insurance` 
#'               - `ninsclasPrivate` 
#'               - `ninsclasPrivate & Medicare`,
#'               data = dt_rhc_ac, family = binomial())
#'dt_rhc_ac$p <-  predict(p, type = "response")
#'Group_index <- assign_int(dt_rhc_ac$p,beta = 4)
#'dt_rhc_ac$subgroup <- Group_index
#'
#'dt_rhc$p <- dt_rhc_ac$p
#'dt_rhc$subgroup <- dt_rhc_ac$subgroup
#'dt_rhc_2 <- dt_rhc[dt_rhc$subgroup==2,]
#'dt_rhc_ac_2 <- dt_rhc_ac[dt_rhc_ac$subgroup == 2,] 
#'dt_rhc_ac_2_cov <- dt_rhc_ac_2[, setdiff(colnames(dt_rhc_ac_2), c("z","subgroup")), drop = FALSE]
#'
#'A_matched_2 <- match_subgroup_auto(Z = dt_rhc_ac_2$z, X = dt_rhc_ac_2_cov, 
#'                           fine = c("ninsclasMedicare",
#'                                    "ninsclasMedicare & Medicaid",
#'                                    "ninsclasNo insurance",
#'                                    "ninsclasPrivate",
#'                                    "ninsclasPrivate & Medicare"),
#'                           data = dt_rhc_2, kmax = 4, method = "maha")
#'}
#'
#'@importFrom mvnfast maha
#'@importFrom stats cov var
#'@export
match_subgroup_auto <- function(Z, X, fine, data, kmax, method = "maha", lambda = 1000) {
  f <- interaction(X[, fine, drop = FALSE], drop = TRUE, lex.order = TRUE)
  kc <- kcalculator(f, Z)
  
  if (is.null(kc$k)) {
    if (sum(Z) <= length(Z) - sum(Z)) {
      A <- matching_nf(Z = Z, X = X, fine = fine, data = data, method = method, lambda = lambda)
    } else {
      Z_pseudo <- 1 - Z
      A <- matching_nf(Z = Z_pseudo, X = X, fine = fine, data = data, method = method, lambda = lambda)
      }
  } else {
    k <- min(kc$k, kmax)   
    if (kc$flip) {
      Z_pseudo <- 1 - Z
      A <- matching_fm(Z = Z_pseudo, X = X, fine = fine, data = data, k = k, method = method)
      } else {
      A <- matching_fm(Z = Z, X = X, fine = fine, data = data, k = k, method = method)
    }
  }
  return(A)
}
