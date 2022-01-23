#'Construct an output for template matching.
#'
#'This function constructs the output for template matching given
#'the relaxsolution to the network flow problem,
#'number of edges in the template-to-treated network,
#'a vector of treatment status, and the original dataset. This function
#'is of little interest to users.
#'
#'
#'
#'@param res A callrelax output.
#'@param num_edges_left Number of edges in the template-to-treatment network.
#'@param Z A vector of treatment status.
#'@param dataset The original dataset.
#'
#'
#'@return  This function returns a list of three objects: 1) feasible: 0/1 depending on the
#'feasibility of the matching problem; 2) match_treated: a data frame of the matched treated
#'units; 3) match_control: a data frame of the matched control units.
#'
#'
#'@export
#'

construct_outcome_template <- function (res, num_edges_left, Z, dataset) {
  if (identical(res, "Not feasible")) {
    cat("The matching is not feasible. Please specify a large caliper and/or k.")
    return(list(feasible = 0,
                match_treated = NULL,
                match_control = NULL))
  }
  else {
    n_t = sum(Z)
    n_c = length(Z) - n_t

    dataset_treated = dataset[Z == 1, ]
    dataset_control = dataset[Z == 0, ]

    # Obtain matched treated
    middle_edges = res$x[seq(num_edges_left - n_t + 1, num_edges_left, 1)]
    treated_ind = which(middle_edges == 1)
    dataset_treated_matched = dataset_treated[treated_ind, ]
    dataset_treated_not_matched = dataset_treated[-treated_ind, ]

    num_pairs_matched = dim(dataset_treated_matched)[1]
    num_treated_not_matched = dim(dataset_treated_not_matched)[1]
    dataset_treated_matched$matched_set = seq(1, num_pairs_matched, 1)
    dataset_treated_not_matched$matched_set = NA

    # Obtain matched control
    rightmost_edges = tail(res$x, n_c)
    control_ind = which(rightmost_edges == 1)
    dataset_control_matched = dataset_control[control_ind, ]
    dataset_control_not_matched = dataset_control[-control_ind, ]

    num_control_not_matched = dim(dataset_control_not_matched)[1]
    dataset_control_matched$matched_set = seq(1, num_pairs_matched, 1)
    dataset_control_not_matched$matched_set = NA

    return(list(feasible = 1,
                match_treated = dataset_treated_matched,
                match_control = dataset_control_matched))
  }
}

