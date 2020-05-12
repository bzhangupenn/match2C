#'Stitch two treated-to-control networks into one two-part networks.
#'
#'This function takes as inputs two networks and one penalty lambda,
#'and constructs one two-part network out of them.
#'
#'This function is of limited interest to users. Once overflow is set to
#'TRUE, each cotnrol in the first network will be directly connected to
#'the sink at a large cost, so that the network flow problem is feasible
#'as long as the first part is feasible.
#'
#'@param net1 A list of five vectors: startn, endn, ucap, cost, b.
#'@param net2 A list of five vectors: startn, endn, ucap, cost, b.
#'@param lambda A penalty.
#'@param controls Number of controls matched to each treated.
#'@param overflow A logical value indicating if overflow protection is turned on.
#'
#'
#'@return  This function returns a list of five vectors:
#'startn, endn, ucap, cost, b.
#'
#'
#'@importFrom utils head tail
#'
#'@export


stitch_two_nets <- function(net1, net2, lambda, controls = 1, overflow = FALSE){

  # Construct startn
  construct_startn <- function(n_t, n_c, net1, net2, overflow = FALSE){
    first_part = net1$startn # source -- treated (net1) -- control (net1) -- treated (net2)
    second_part = tail(net2$startn, -n_c) + n_t + n_c # treated (net2) -- control (net2) -- sink
    # If overflow is turned on, we further connect each treated in net2
    # directly to the sink
    if (overflow == TRUE)
      overflow_part = seq(2+n_t+n_c, n_c + 1 + n_t + n_c)
    else overflow_part = NULL
    return(c(first_part, second_part, overflow_part))
  }

  construct_endn <- function(n_t, n_c, net1, net2, overflow = FALSE){
    first_part = head(net1$endn, -n_c)# source -- treated (net1) -- control (net1)
    second_part = net2$endn + n_t + n_c # control (net1) -- treated (net2) -- control (net2) -- sink
    if (overflow == TRUE)
      overflow_part = rep(2+2*n_t+2*n_c, n_c)
    else overflow_part = NULL
    return(c(first_part, second_part, overflow_part))
  }

  construct_capn <- function(n_t, n_c, net1, net2, overflow = FALSE){
    if (overflow == TRUE) overflow_part = rep(1, n_c)
    else overflow_part = NULL
    return(c(net1$ucap, tail(net2$ucap, -n_c), overflow_part))
  }

  construct_costn <- function(n_t, n_c, net1, net2, overflow = FALSE){
    first_part = head(net1$cost, -n_c) # source -- treated (net1) -- control (net1)
    second_part = rep(0, n_c) # control (net1) -- treated (net2)
    third_part = lambda*tail(net2$cost, -n_c) # treated (net2) -- control (net2) -- sink
    if (overflow == TRUE) overflow_part = rep(lambda*100, n_c)
    else overflow_part = NULL
    return(c(first_part, second_part, third_part, overflow_part))
  }

  n_t = net1$b[1]/controls # number of treated in net1 == number of control in net2
  n_c = net2$b[1]/controls # number of treated in net2 == number of control in net1
  num_node = 1 + 2*n_t + 2*n_c + 1

  net = list(startn = construct_startn(n_t, n_c, net1, net2, overflow),
             endn = construct_endn(n_t, n_c, net1, net2, overflow),
             ucap = construct_capn(n_t, n_c, net1, net2, overflow),
             cost = construct_costn(n_t, n_c, net1, net2, overflow),
             b = c(controls*n_t, rep(0, num_node-2), - controls*n_t))
  return(net)
}
