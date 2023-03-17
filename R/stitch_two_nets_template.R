#'Stitch a template-to-treated network and a treated-to-control network
#'into one two-part network.
#'
#'This function takes as inputs a template-to-treated network,
#'one treated-to-control network, a tuning parameter lambda,
#'and number of controls, and constructs one two-part network out of them.
#'
#'This function is of limited interest to users. Parameter lambda is a
#'weight given to the first part of the network, and a large lambda value
#'emphasizes resemblance to the template. Parameter multiple could be taken
#'as any integer number between 1 and floor(treated size / template size).
#'
#'@param net1 A list of five vectors: startn, endn, ucap, cost, b.
#'@param net2 A list of five vectors: startn, endn, ucap, cost, b.
#'@param n_c Number of control units.
#'@param lambda A penalty.
#'@param multiple Number of treated units matched to each unit in the template
#'
#'
#'@return  This function returns a list of five vectors:
#'startn, endn, ucap, cost, b.
#'
#'
#'@importFrom utils head tail
#'
#'@export


stitch_two_nets_template <- function(net1, net2, n_c, lambda, multiple = 1) {
  construct_startn <- function(n_temp, n_t, n_c, net1, net2) {
    first_part = net1$startn
    second_part = tail(net2$startn, -n_t) + n_temp + n_t
    return(c(first_part, second_part))
  }
  construct_endn <- function(n_temp, n_t, n_c, net1, net2) {
    first_part = head(net1$endn, -n_t)
    second_part = net2$endn + n_temp + n_t
    return(c(first_part, second_part))
  }
  construct_capn <- function(n_temp, n_t, n_c, net1, net2) {
    return(c(net1$ucap, tail(net2$ucap, -n_t)))
  }
  construct_costn <- function(n_temp, n_t, n_c, net1, net2) {
    first_part = 100*head(net1$cost, -n_t)
    second_part = rep(0, n_t)
    third_part = 100*lambda*tail(net2$cost, -n_t)
    return(c(first_part, second_part, third_part))
  }
  n_temp = net1$b[1]/multiple
  n_t = net2$b[1]
  n_c = n_c
  num_node = 1 + n_temp + 2 * n_t + n_c + 1
  net = list(startn = construct_startn(n_temp, n_t, n_c, net1, net2),
             endn = construct_endn(n_temp, n_t, n_c, net1, net2),
             ucap = construct_capn(n_temp, n_t, n_c, net1, net2),
             cost = construct_costn(n_temp, n_t, n_c, net1, net2),
             b = c(multiple * n_temp, rep(0, num_node - 2), -multiple * n_temp))
  return(net)
}

