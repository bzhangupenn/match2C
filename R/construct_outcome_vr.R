#'Construct an output for matching
#'
#'
#'@param res A callrelax output.
#'@param dist_list A possibly sparse representation of the first distance matrix.
#'@param Z A vector of treatment status.
#'@param dataset The original dataset.
#'
#'
#'@return a list comprising the following three elements: 
#'feasibility: 0/1 depend on the feasibility of the matching problem.
#'
#'data_with_matched_to_ind: a data frame that is the same as the original data frame, 
#'except that a column called "matched_to", a column called "distance" and a column called 
#'"matched_to_index" are added to it. The variable "matched_to" records the order of the 
#'treated unit within the treatment group for each matched control, and assigns NA 
#'to controls not matched to any treated unit. The variable "matched_to_ind" then converts 
#'this order into the original observation index. Variable "distance" records the 
#'control-to-treated distance between one matched control and the corresponding 
#'treated unit, and assigns NA to controls that are left unmatched and all treated 
#'units. If matching is not feasible, NULL will be returned.
#'
#'data_without_unmatched_controls: a data frame that is the same as data_with_matched_to_ind, 
#'except that the unmatched controls are removed. If matching is not feasible, NULL will be returned.
#'
#'
#'@export
#'
#'
construct_outcome_vr <- function(res, dist_list, Z, dataset){
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
    
    source_uf = length(unique(dist_list$start_n))
    route = res$x[seq(source_uf + 1, length(dist_list$d) + source_uf, 1)]
    id_t = dist_list$start_n[route == 1]
    id_c = dist_list$end_n[route == 1]
    d_c = dist_list$d[route == 1]
    
    dataset_treated = dataset[Z == 1,]
    dataset_treated$matched_to = NA
    dataset_treated$distance = NA
    dataset_treated$matched_to_index = NA
    
    dataset_control$matched_to = NA
    dataset_control$distance = NA
    dataset_control$matched_to[id_c - n_t] = id_t
    dataset_control$distance[id_c - n_t] = d_c
    dataset_control$matched_to_index = NA
    is.real.treated <- dataset_control$matched_to < n_t + 1
    dataset_control$matched_to_index[is.real.treated] = 
      rownames(dataset_treated)[dataset_control$matched_to[is.real.treated]]
    
    dataset_matched_treated = dataset_treated
    dataset_matched_treated$matched_to_index = NA
    
    dataset_matched_controls <- dataset_control[dataset_control$matched_to < n_t + 1,]
    dataset_matched_controls$matched_to_index = rownames(dataset_treated)[dataset_matched_controls$matched_to]
    
    
    # combine dataset_treated, dataset_matched_control, and dataset_not_matched
    # and return dataset_with_matched_set
    # Note dataset_with_matched_set is in the orignal order
    dataset_with_matched_set = rbind(dataset_treated, dataset_control)
    dataset_with_matched_set = dataset_with_matched_set[order(as.numeric(row.names(dataset_with_matched_set))),]
    
    # The same dataframe except that ordered by the matched_set indeces
    dataset_exclude_unmatched_controls = rbind(dataset_matched_treated, dataset_matched_controls)
    dataset_exclude_unmatched_controls = dataset_exclude_unmatched_controls[order(as.numeric(row.names(dataset_exclude_unmatched_controls))),]
    
    return(list(feasible = 1,
                data_with_matched_to_ind = dataset_with_matched_set,
                data_without_unmatched_controls = dataset_exclude_unmatched_controls))
  }
}
