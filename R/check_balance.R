#'Check balance after matching.
#'
#'This function checks the overall balance after statistical matching and plots
#'the distribution of the propensity score in the treated group, the control group,
#'and the matched control group.
#'
#'
#'
#'@param Z A vector of treatment indicator.
#'@param match_object An object returned by match_2C or match_2C_mat or match_2C_list.
#'@param cov_list A vector of names of covariates as appeared in the original dataset.
#'@param plot_propens Post-matching distribution of the estimated propensity scores in
#'                    two groups is plotted if TRUE; FALSE by default.
#'@param propens NULL by default. If plot_propens = TRUE, then a vector of
#'               estimated propensity scores satisfying length(propens) = length(Z)
#'               needs to be supplied.
#'
#'
#'@return  This function returns a data frame of the overall balance after statistical
#'matching. We tabulate the mean of each covariate in the cov_list in the treated group
#'and control groups after matching, and calculate their standardized differences.
#'Standardized difference is defined as the mean difference divided by the pooled
#'standard error before matching.
#'
#'@importFrom ggplot2 ggplot aes geom_density theme_bw theme
#'@importFrom stats sd
#'@export
#'

check_balance <- function(Z, match_object, cov_list, plot_propens, propens){

  # Compute before-matching balance statistics
  dt_treated_before <- match_object$data_with_matched_set_ind[Z==1, cov_list]
  dt_control_before <- match_object$data_with_matched_set_ind[Z==0, cov_list]

  # Add pscore
  #dt_treated_before$propensity_score = propens[Z==1]
  #dt_control_before$propensity_score = propens[Z==0]

  mean_treated_before = apply(dt_treated_before, 2, mean)
  mean_control_before = apply(dt_control_before, 2, mean)

  mean_diff_before = mean_treated_before - mean_control_before

  sd_treated_before = apply(dt_treated_before, 2, stats::sd)
  sd_control_before = apply(dt_control_before, 2, stats::sd)

  pooled_sd = sqrt(sd_treated_before^2 + sd_control_before^2)

  std_before = mean_diff_before/pooled_sd

  # Compute after-matching balance statistics
  dt_matched = match_object$data_with_matched_set_ind
  #dt_matched$propensity_score = propens
  Z_matched = Z[!is.na(dt_matched$matched_set)]
  dt_matched = dt_matched[!is.na(dt_matched$matched_set), cov_list]
  dt_treated_after <- dt_matched[Z_matched==1, ]
  dt_control_after <- dt_matched[Z_matched==0, ]

  mean_treated_after = apply(dt_treated_after, 2, mean)
  mean_control_after = apply(dt_control_after, 2, mean)

  mean_diff_after = mean_treated_after - mean_control_after

  std_after = mean_diff_after/pooled_sd

  # Tabulate all results
  balance_table = data.frame(mean_treated_before, mean_control_before, std_before,
                             mean_control_after, std_after)
  rownames(balance_table) <= cov_list
  colnames(balance_table) <- c('Z = 1', 'Z = 0 (Bef)',
                               'Std. Diff (Bef)',
                               'Z = 0 (Aft)',
                               'Std. Diff (Aft)')

  # Plot the propensity score distribution if plot_propens == TRUE
  if (plot_propens) {
    propens_treated = propens[Z == 1]
    propens_all_control = propens[Z == 0]
    propens_matched_control = propens[Z == 0 & (!is.na(match_object$data_with_matched_set_ind$matched_set))]


    group = c(rep('Treated', length(propens_treated)),
              rep('All controls', length(propens_all_control)),
              rep('Matched controls', length(propens_matched_control)))


    propensity = c(propens_treated, propens_all_control, propens_matched_control)

    propens_df = data.frame(group, propensity)
    pt = ggplot2::ggplot(data = propens_df, ggplot2::aes(x = propensity, color = group)) +
      ggplot2::geom_density(size = 1.5) + ggplot2::theme_bw(base_size = 20) +
      ggplot2::theme(legend.position = 'top')
    print(pt)
  }

  return(balance_table)
}

