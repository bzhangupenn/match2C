#'Construct an output for matching.
#'
#'This function constructs the output given the relaxsolution
#'to the associated network flow problem and the original dataset.
#'
#'
#'
#'@param res A callrelax output.
#'@param dist_list_1 A possibly sparse representation of the first distance matrix.
#'@param Z A vector of treatment status.
#'@param dataset The original dataset.
#'@param controls Number of controls matched to each treated.
#'
#'
#'@return  This function returns a list of three objects: 1) feasible: 0/1 depending on the
#'feasibility of the matching problem; 2) data_with_matched_set_ind: a dataframe that is the
#'same as the original dataframe, except that a column called ``matched_set'' and a column
#'called ``distance'' are added to it. ``matched_set'' column assigns 1,2,...,n_t to each matched
#'set, and NA to those not matched to any treated. Variable ``distance'' records the distance
#'(as specified in the left network) between each matched control and the treated, and assigns NA
#'to all treated and cotnrols that are left unmatched. If matching is not feasible, NULL will be returned;
#'3) matched_data_in_order:a dataframe organized in the order of matched sets and otherwise the
#'same as data_with_matched_set_ind. Note that the matched_set column assigns 1,2,...,n_t for
#' as indices for matched sets, and NA for those controls that are not paired. Null will be returned
#' if the matching is infeasible.
#'
#'
#'@export
#'

construct_outcome <- function(res, dist_list_1, Z, dataset, controls = 1){
  if (identical(res, 'Not feasible')){
    cat('The matching is not feasible. Please specify a large caliper and/or k.')
    return(list(feasible = 0,
                data_with_matched_set_ind = NULL,
                matched_data_in_order = NULL,
                balance_table = NULL))
  } else{

    # Extract the matched controls and their distances
    n_t = sum(Z)
    n_c = length(Z) - n_t
    dataset_control = dataset[Z == 0,]

    route = res$x[seq(n_t + 1, length(dist_list_1$d) + n_t, 1)]
    id_c = dist_list_1$end_n[route == 1]
    d_c = dist_list_1$d[route == 1]

    dataset_matched_controls = dataset_control[id_c - n_t,]

    # Add matched-set and distance information to matched controls
    dataset_matched_controls$matched_set = rep(seq(1, n_t, 1), each = controls)
    dataset_matched_controls$distance = d_c

    # Add matched-set and distance information to the treated
    dataset_treated = dataset[Z == 1,]
    dataset_treated$matched_set = seq(1, n_t, 1)
    dataset_treated$distance = rep(NA, n_t)

    # Add NA to those not matched
    dataset_not_matched = dataset_control[setdiff(seq(1,n_c,1), id_c - n_t),]
    dataset_not_matched$matched_set = NA
    dataset_not_matched$distance = NA

    # combine dataset_treated, dataset_matched_control, and dataset_not_matched
    # and return dataset_with_matched_set
    # Note dataset_with_matched_set is in the orignal order
    dataset_with_matched_set = rbind(dataset_treated, dataset_matched_controls, dataset_not_matched)
    dataset_with_matched_set = dataset_with_matched_set[order(as.numeric(row.names(dataset_with_matched_set))),]

    # The same dataframe except that ordered by the matched_set indeces
    dataset_with_matched_set_in_order = dataset_with_matched_set[order(as.numeric(dataset_with_matched_set$matched_set)),]


    return(list(feasible = 1,
                data_with_matched_set_ind = dataset_with_matched_set,
                matched_data_in_order = dataset_with_matched_set_in_order))
  }
}
