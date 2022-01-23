#'Check balance after template matching.
#'
#'This function checks the overall balance after template matching and returns a
#'dataframe with 7 columns: (1) mean of all covariates in the treated group,
#'(2) mean of all covariates in the control group,
#'(3) standardized mean differences of (1) and (2),
#'(4) mean of all covariates in the matched treated group,
#'(5) mean of all covariates in the matched control group,
#'(6) standardized mean differences of (4) and (5),
#'(7) mean of covariates in the template
#'
#'
#'@param dataset The original dataset.
#'@param template A data frame of the template.
#'@param template_match_object An object returned by template_match.
#'@param cov_list A vector of names of covariates as appeared in the original dataset and the template.
#'
#'
#'@return  This function returns a data frame of the overall balance after template
#'matching. We tabulate the mean and SMD of each covariate in the cov_list in the template,
#'the matched treated group, and the matched control group.
#'
#'@importFrom stats sd
#'@export
#'

check_balance_template <- function(dataset, template, template_match_object, cov_list){

  # Compute before-matching balance statistics
  dt_treated_before <- dataset[dataset$Z == 1, cov_list]
  dt_control_before <- dataset[dataset$Z == 0, cov_list]
  n_t = dim(dt_treated_before)[1]
  n_c = dim(dt_control_before)[1]

  mean_treated_before = apply(dt_treated_before, 2, mean)
  mean_control_before = apply(dt_control_before, 2, mean)

  mean_diff_before = mean_treated_before - mean_control_before

  sd_treated_before = apply(dt_treated_before, 2, stats::sd)
  sd_control_before = apply(dt_control_before, 2, stats::sd)

  pooled_sd = sqrt(sd_treated_before^2 + sd_control_before^2)

  std_before = mean_diff_before/pooled_sd

  # Compute after-matching balance statistics

  dt_treated_after <- template_match_object$match_treated[, cov_list]
  dt_control_after <- template_match_object$match_control[, cov_list]
  n_match = dim(dt_treated_after)[1]

  mean_treated_after = apply(dt_treated_after, 2, mean)
  mean_control_after = apply(dt_control_after, 2, mean)

  mean_diff_after = mean_treated_after - mean_control_after

  std_after = mean_diff_after/pooled_sd


  # Compute template statistics
  template_relevant = template[, intersect(cov_list, colnames(template))]
  n_template = dim(template_relevant)[1]
  mean_template = apply(template_relevant, 2, mean)
  mean_template_df = data.frame(mean_template)
  rownames(mean_template_df) = intersect(cov_list, colnames(template))

  # Tabulate all results
  balance_table = data.frame(mean_treated_before, mean_control_before, std_before,
                             mean_treated_after, mean_control_after, std_after)
  rownames(balance_table) <= cov_list

  final_table = merge(balance_table, mean_template_df, by = 0, all = TRUE)
  row.names(final_table) = final_table$Row.names
  final_table = final_table[2:length(final_table)]
  #final_table = rbind(n = c(n_t, n_c, NA, n_match, n_match, NA, n_template),
  #                    final_table)
  colnames(final_table) <- c('Z = 1 (Bef)', 'Z = 0 (Bef)',
                               'Std. Diff (Bef)',
                               'Z = 1 (Aft)', 'Z = 0 (Aft)',
                               'Std. Diff (Aft)',
                               'Template')
  #final_table[1,] = as.integer(final_table[1,])

  return(final_table)
}

