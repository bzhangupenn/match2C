## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  fig.width=7.5, 
  fig.height=5, 
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
options(scipen = 99)
options(digits = 3)
library(match2C)
library(ggplot2)
library(mvtnorm)

## ---- echo=TRUE, warning=FALSE, message=FALSE---------------------------------
attach(dt_Rouse)
X = cbind(female,black,bytest,dadeduc,momeduc,fincome) # covariates to be matched
Z = IV # IV-defined exposure in this dataset

# Fit a propensity score model
propensity = glm(IV~female+black+bytest+dadeduc+momeduc+fincome,
                 family=binomial)$fitted.values

# Number of treated and control
n_t = sum(Z) # 1,122 treated
n_c = length(Z) - n_t # 1,915 control

dt_Rouse$propensity = propensity
detach(dt_Rouse)

## ----intro example, echo=FALSE------------------------------------------------
dist_list_pscore = create_list_from_scratch(Z, X, exact = NULL, p = propensity, 
                                       caliper_low = 0.008, 
                                       k = NULL, 
                                       method = 'robust maha')

matching_output_example = match_2C_list(Z, dt_Rouse, 
                                        dist_list_pscore, 
                                        dist_list_2 = NULL, 
                                        controls = 1)

## ----intro example feasible, echo=TRUE----------------------------------------
# Check feasibility
matching_output_example$feasible

## ----intro example data1, echo=TRUE-------------------------------------------
# Check the original dataset with two new columns
head(matching_output_example$data_with_matched_set_ind, 6)

## ----intro example data2, echo=TRUE-------------------------------------------
# Check dataframe organized in matched set indices
head(matching_output_example$matched_data_in_order, 6)

## ----intro example check balance table, echo=TRUE-----------------------------
tb_example = check_balance(Z, matching_output_example, 
              cov_list = c('female', 'black', 'bytest', 'fincome', 'dadeduc', 'momeduc', 'propensity'),
              plot_propens = FALSE)
print(tb_example)

## ----intro example check balance plot, echo=TRUE------------------------------
tb_example = check_balance(Z, matching_output_example, 
              cov_list = c('female', 'black', 'bytest', 'fincome', 
                           'dadeduc', 'momeduc', 'propensity'),
              plot_propens = TRUE, propens = propensity)

## ----match2C no caliper, echo=TRUE--------------------------------------------
# Perform a matching with minimal input
matching_output = match_2C(Z = Z, X = X, method = 'robust maha',
                           propensity = propensity, 
                           dataset = dt_Rouse)
tb = check_balance(Z, matching_output, 
                   cov_list = c('female', 'black', 'bytest', 'fincome', 'dadeduc', 'momeduc', 'propensity'),
                   plot_propens = TRUE, propens = propensity)
print(tb)

## ----match2c exact------------------------------------------------------------
# Perform a matching with minimal input
matching_output_with_exact = match_2C(Z = Z, X = X, exact = c('dadeduc', 'momeduc'),
                           propensity = propensity, 
                           dataset = dt_Rouse)

# Check exact matching
head(matching_output_with_exact$matched_data_in_order[, c('female', 'black', 'bytest', 
                                      'fincome', 'dadeduc', 'momeduc', 
                                      'propensity', 'IV', 'matched_set')])

# Check overall balance
tb = check_balance(Z, matching_output_with_exact, 
                   cov_list = c('female', 'black', 'bytest', 'fincome', 'dadeduc', 'momeduc', 'propensity'),
                   plot_propens = TRUE, propens = propensity)


## ----match2C fine balance, echo=TRUE------------------------------------------
# Perform a matching with fine balance
matching_output2 = match_2C(Z = Z, X = X, 
                            propensity = propensity, 
                            dataset = dt_Rouse,
                            fb_var = c('dadeduc'))

## ----match2C fine balance check, echo=TRUE------------------------------------
# Perform a matching with fine balance
tb2 = check_balance(Z, matching_output2, 
                   cov_list = c('female', 'black', 'bytest', 'fincome', 'dadeduc', 'momeduc', 'propensity'),
                   plot_propens = TRUE, propens = propensity)
print(tb2)

## ----match2C fine balance 2, echo=TRUE----------------------------------------
# Perform a matching with fine balance on dadeduc and moneduc
matching_output3 = match_2C(Z = Z, X = X, 
                            propensity = propensity, 
                            dataset = dt_Rouse,
                            fb_var = c('dadeduc', 'momeduc'))
tb3 = check_balance(Z, matching_output2, 
                   cov_list = c('female', 'black', 'bytest', 'fincome', 'dadeduc', 'momeduc', 'propensity'),
                   plot_propens = FALSE)
print(tb3)

## ----match2C with or without caliper speed, echo=TRUE-------------------------
# Timing the vanilla match2C function
ptm <- proc.time()
matching_output2 = match_2C(Z = Z, X = X, 
                            propensity = propensity, 
                            dataset = dt_Rouse)
time_vanilla = proc.time() - ptm

# Timing the match2C function with caliper on the left
ptm <- proc.time()
matching_output_one_caliper = match_2C(Z = Z, X = X, propensity = propensity, 
                            caliper_left = 0.05, caliper_right = 0.05, 
                            k_left = 100,
                            dataset = dt_Rouse)
time_one_caliper = proc.time() - ptm

# Timing the match2C function with caliper on the left and right
ptm <- proc.time()
matching_output_double_calipers = match_2C(Z = Z, X = X, 
                            propensity = propensity, 
                            caliper_left = 0.05, caliper_right = 0.05, 
                            k_left = 100, k_right = 100,
                            dataset = dt_Rouse)
time_double_caliper = proc.time() - ptm

rbind(time_vanilla, time_one_caliper, time_double_caliper)[,1:3]

## ----match2C small caliper fail, echo=TRUE------------------------------------
# Perform a matching with fine balance on dadeduc and moneduc
matching_output_unfeas = match_2C(Z = Z, X = X, propensity = propensity, 
                                  dataset = dt_Rouse,
                                  caliper_left = 0.001)

## ----match2C force control, echo=TRUE-----------------------------------------

# Create a binary vector with 1's in the first 100 entries and 0 otherwise
# length(include_vec) = n_c

include_vec = c(rep(1, 100), rep(0, n_c - 100))
# Perform a matching with minimal input
matching_output_force_include = match_2C(Z = Z, X = X, 
                           propensity = propensity, 
                           dataset = dt_Rouse, 
                           include = include_vec)

## ----match2C force control 2, echo=TRUE---------------------------------------

matched_data = matching_output_force_include$data_with_matched_set_ind
matched_data_control = matched_data[matched_data$IV == 0,]
head(matched_data_control) # Check the matched_set column

## ----mat no caliper, echo=TRUE------------------------------------------------
# Construct a distance matrix based on Mahalanobis distance
dist_mat_1 = optmatch::match_on(IV~female+black+bytest+dadeduc+momeduc+fincome, 
                      method = 'mahalanobis', data = dt_Rouse)

# Construct a second distance matrix based on variable dadeduc
dist_mat_2 = optmatch::match_on(IV ~ dadeduc, method = 'euclidean', 
                                data = dt_Rouse)
matching_output_mat = match_2C_mat(Z, dt_Rouse, dist_mat_1, dist_mat_2, 
                               lambda = 10000, controls = 1,
                               p_1 = NULL, p_2 = NULL)

# Examine the balance after matching
tb_mat = check_balance(Z, matching_output_mat, 
                   cov_list = c('female', 'black', 'bytest', 'fincome', 'dadeduc', 'momeduc', 'propensity'),
                   plot_propens = FALSE)

print(tb_mat)

## ----mat with caliper, echo=TRUE----------------------------------------------
matching_output_mat_caliper = match_2C_mat(Z, dt_Rouse, 
                                 dist_mat_1, dist_mat_2, 
                                 lambda = 100000, controls = 1,
                                 p_1 = propensity, 
                                 caliper_1 = 0.05, k_1 = 100)

## ----list create list from mat no caliper, echo = TRUE------------------------
dist_mat_1 = optmatch::match_on(IV ~ female + black + bytest + 
                                dadeduc + momeduc + fincome, 
                                method = 'mahalanobis', data = dt_Rouse)

list_0 = create_list_from_mat(Z, dist_mat_1, p = NULL) 
length(list_0$start_n) # number of edges in the network
identical(length(list_0$start_n), n_t*n_c) # Check # of edges is n_t * n_c

## ----list create list from mat with caliper, echo = TRUE----------------------
list_1 = create_list_from_mat(Z, dist_mat_1, 
                              p = propensity, 
                              caliper = 0.05)
length(list_1$start_n) # Number of edges is almost halved

## ----list from scratch ex1, echo=TRUE-----------------------------------------
# Mahalanobis distance on all variables
dist_list_vanilla_maha = create_list_from_scratch(Z, X, exact = NULL, 
                                                  method = 'maha') 

# Hamming distance on all variables
dist_list_vanilla_Hamming = create_list_from_scratch(Z, X, exact = NULL, 
                                                      method = 'Hamming') 

# Robust Mahalanobis distance on all variables
dist_list_vanilla_robust_maha = create_list_from_scratch(Z, X, exact = NULL, 
                                                      method = 'robust maha') 


## ----list from scratch ex2, echo=TRUE-----------------------------------------
# Mahalanobis distance on all variables with pscore caliper
dist_list_pscore_maha = create_list_from_scratch(Z, X, exact = NULL, 
                                       p = propensity, 
                                       caliper_low = 0.05, 
                                       k = 100, 
                                       method = 'maha') 

# Hamming distance on all variables with pscore caliper
dist_list_pscore_Hamming = create_list_from_scratch(Z, X, exact = NULL, 
                                       p = propensity, 
                                       caliper_low = 0.05, 
                                       k = 100, 
                                       method = 'Hamming') 

# Robust Mahalanobis distance on all variables with pscore caliper
dist_list_pscore_robust_maha = create_list_from_scratch(Z, X, exact = NULL, 
                                       p = propensity, 
                                       caliper_low = 0.05, 
                                       k = 100, 
                                       method = 'robust maha') 

## ----list from scratch ex small caliper fail, echo=TRUE-----------------------
dist_list_pscore_maha_hard = create_list_from_scratch(Z, X, exact = NULL, 
                                       p = propensity, 
                                       caliper_low = 0.001, 
                                       method = 'maha') 


## ----list from scratch ex soft caliper, echo=TRUE-----------------------------
dist_list_pscore_maha_soft = create_list_from_scratch(Z, X, exact = NULL, 
                                       p = propensity, 
                                       caliper_low = 0.001, 
                                       method = 'maha', 
                                       penalty = 1000) 


## ----list from scratch ex3, echo=TRUE-----------------------------------------
dist_list_exact_dadeduc_maha = create_list_from_scratch(Z, X, 
                                                        exact = c('dadeduc'), 
                                                        method = 'maha') 

## ----list from scratch ex4, echo=TRUE-----------------------------------------
dist_list_exact_dad_mom_with_caliper = create_list_from_scratch(Z, X, 
                                                exact = c('dadeduc', 'momeduc'), 
                                                p = propensity, 
                                                caliper_low = 0.05, 
                                                caliper_high = 0.1,
                                                k = 100, 
                                                method = 'maha') 

## ----dist list ex1, echo=TRUE-------------------------------------------------
# Construct a distance list representing the network structure on the left.
dist_list_pscore = create_list_from_scratch(Z, X, exact = NULL, 
                                       p = propensity, 
                                       caliper_low = 0.008, 
                                       k = NULL, 
                                       method = 'maha')

# Perform matching. Set dist_list_2 = NULL as we are 
# performing a bipartite matching.
matching_output_ex1 = match_2C_list(Z, dt_Rouse, dist_list_pscore, 
                                    dist_list_2 = NULL, 
                                    controls = 1)

## ----dist list ex2, echo=TRUE-------------------------------------------------
# Mahalanobis distance on all variables; no caliper
dist_list_no_caliper = create_list_from_scratch(Z, X, exact = NULL, 
                                                p = NULL, 
                                                method = 'maha')

# Connect treated to controls within a stringent propensity score caliper.
# We use a soft caliper here to ensure feasibility.
dist_list_2 = create_list_from_scratch(Z = Z, X = rep(1, length(Z)), 
                                       exact = NULL,
                                       p = propensity, 
                                       caliper_low = 0.002, 
                                       method = 'L1', 
                                       k = NULL,
                                       penalty = 100)

matching_output_ex2 = match_2C_list(Z, dt_Rouse, 
                                    dist_list_no_caliper, 
                                    dist_list_2, 
                                    lambda = 1000, controls = 1)

## ----dist list ex3_1, echo=TRUE-----------------------------------------------
# Mahalanobis distance with exact matching on dadeduc and momeduc
dist_list_1 = create_list_from_scratch(Z, X, exact = c('dadeduc', 'momeduc'), 
                                       p = propensity, caliper_low = 0.05, 
                                       method = 'maha')

matching_output_ex3_1 = match_2C_list(Z, dt_Rouse, dist_list_1, 
                                  dist_list_2 = NULL, lambda = NULL)

## ----dist list ex3_2, echo=TRUE-----------------------------------------------
# Maha distance with exact matching on dadeduc and momeduc
dist_list_1 = create_list_from_scratch(Z, X, 
                                       exact = c('dadeduc', 'momeduc'), 
                                       method = 'maha')

# Maha distance on all other variables
dist_list_2 = create_list_from_scratch(Z, X[, c('female', 'black', 'bytest', 'fincome')], 
                                       p = propensity, 
                                       caliper_low = 0.05, 
                                       method = 'maha')

matching_output_ex3_2 = match_2C_list(Z, dt_Rouse, dist_list_1, dist_list_2, lambda = 100)

## ----dist list ex3 check balance, echo=TRUE-----------------------------------
tb_ex3_1 = check_balance(Z, matching_output_ex3_1, 
                        cov_list = c('female', 'black', 'bytest', 'fincome', 'dadeduc', 'momeduc', 'propensity'),
                        plot_propens = TRUE, propens = propensity)

print(tb_ex3_1)

tb_ex3_2 = check_balance(Z, matching_output_ex3_2, 
                        cov_list = c('female', 'black', 'bytest', 'fincome', 'dadeduc', 'momeduc', 'propensity'),
                        plot_propens = TRUE, propens = propensity)
print(tb_ex3_2)


## ----dist list ex4, echo=TRUE-------------------------------------------------
dist_list_1 = create_list_from_scratch(Z = Z, X = X, 
                                       exact = c('female', 'black'),
                                       p = propensity,
                                       caliper_low = 0.15,
                                       method = 'maha')

dist_list_2 = create_list_from_scratch(Z = Z, X = X[, c('dadeduc', 'momeduc')],
                                                  method = '0/1')
matching_output_ex4 = match_2C_list(Z, dt_Rouse, dist_list_1, dist_list_2, lambda = 1000)

tb_ex4 = check_balance(Z, matching_output_ex4, 
                        cov_list = c('female', 'black', 'bytest', 'fincome', 'dadeduc', 'momeduc', 'propensity'),
                        plot_propens = TRUE, propens = propensity)
print(tb_ex4)


## ----dist list ex5, echo=TRUE-------------------------------------------------

# The first 50 controls must be included
include_vec = c(rep(1, 50), rep(0, n_c - 50))

# dist_list_1 and dist_list_2 are from example 4 above
dist_list_1_update = force_control(dist_list_1, Z = Z, include = include_vec)
dist_list_2_update = force_control(dist_list_2, Z = Z, include = include_vec)

matching_output_ex5 = match_2C_list(Z, dt_Rouse, 
                                    dist_list_1_update, 
                                    dist_list_2_update, lambda = 1000)

tb_ex5 = check_balance(Z, matching_output_ex5, 
                        cov_list = c('female', 'black', 'bytest', 'fincome', 'dadeduc', 'momeduc', 'propensity'),
                        plot_propens = TRUE, propens = propensity)
print(tb_ex5)

## ----dist list ex6, echo=TRUE-------------------------------------------------
# Construct distance list on the left: Maha subject to pscore caliper
dist_list_left = create_list_from_scratch(Z = Z, X = X,
                                            exact = c('black'),
                                            p = propensity,
                                            caliper_low = 0.2,
                                            method = 'maha')

# Construct distance list on the right: L1 distance on the pscore plus fine balance
dist_list_right_pscore = create_list_from_scratch(Z = Z, X = propensity,
                                             p = propensity,
                                             caliper_low = 1,
                                             method = 'L1')

dist_list_right_fb = create_list_from_scratch(Z = Z, X = X[, c('dadeduc')],
                                                  p = propensity,
                                                  caliper_low = 1,
                                                  method = '0/1')

dist_list_right = dist_list_right_pscore
dist_list_right$d = dist_list_right$d + 100*dist_list_right_fb$d 


# The first 50 controls must be included
include_vec = c(rep(1, 50), rep(0, n_c - 50))

# dist_list_1 and dist_list_2 are from example 4 above
dist_list_left_update = force_control(dist_list_left, Z = Z, include = include_vec)
dist_list_right_update = force_control(dist_list_right, Z = Z, include = include_vec)

matching_output_ex6 = match_2C_list(Z = Z, dataset = dt_Rouse,
                                  dist_list_1 = dist_list_left_update,
                                  dist_list_2 = dist_list_right_update,
                                  lambda = 100)

tb_ex6 = check_balance(Z, matching_output_ex6, 
                        cov_list = c('female', 'black', 'bytest', 'fincome', 'dadeduc', 'momeduc', 'propensity'),
                        plot_propens = TRUE, propens = propensity)
print(tb_ex6)

## ----template match generate data, echo = TRUE--------------------------------
set.seed(123)
ratio = 3 # Control-to-treated ratio
n_t = 500 # 500 treated units
n_c = n_t * ratio # 1500 control units
p = 10 # Number of covariates

# Generate covariates for the treated and control units
X_treated = rmvnorm(n_t, mean = c(1, rep(0, p - 1)), sigma = diag(p))
X_control = rmvnorm(n_c, mean = rep(0.5, p), sigma = diag(p))

X = rbind(X_treated, X_control) # 2000-by-10 matrix of covariates
colnames(X) = c('V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'V9', 'V10')
Z = c(rep(1, n_t), rep(0, n_c)) # Length-2000 vector of treatment status
dataset = data.frame(X, Z) # The original dataset

head(X)

## ----template match generate template, echo = TRUE----------------------------
n_template = 100
beta = 0.4
d = 5
template = as.data.frame(rmvnorm(n_template, mean = c(beta, rep(0,d-1))))
head(template)

## ----template basic 1, echo = TRUE--------------------------------------------
template_match_res = template_match(template = template, X = X, Z = Z, 
                                     dataset = dataset, multiple = 1, lambda = 0,
                                     caliper_pscore = 0.1, penalty_pscore = Inf)

check_balance_template(dataset = dataset,
                       template = template, 
                      template_match_object = template_match_res, 
                      cov_list = c('V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'V9', 'V10'))

## ----template basic 1 300, echo = TRUE----------------------------------------
template_match_res1 = template_match(template = template, X = X, Z = Z, 
                                     dataset = dataset, multiple = 3, lambda = 0,
                                     caliper_pscore = 0.1, penalty_pscore = Inf)

check_balance_template(dataset = dataset,
                       template = template, 
                      template_match_object = template_match_res1, 
                      cov_list = c('V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'V9', 'V10'))

## ----template basic 2, echo = TRUE--------------------------------------------
template_match_res2 = template_match(template = template, X = X, Z = Z, 
                                     dataset = dataset, multiple = 1, lambda = 1000,
                                     caliper_pscore = 0.2, penalty_pscore = Inf)

check_balance_template(dataset = dataset,
                       template = template, 
                      template_match_object = template_match_res2, 
                      cov_list = c('V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'V9', 'V10'))

## ----template basic 3, echo = TRUE--------------------------------------------
template_match_res3 = template_match(template = template, X = X, Z = Z, 
                                     dataset = dataset, multiple = 1, lambda = 100,
                                     caliper_gscore = 0.02, penalty_gscore = 100)


check_balance_template(dataset = dataset,
                       template = template, 
                      template_match_object = template_match_res3, 
                      cov_list = c('V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'V9', 'V10'))

