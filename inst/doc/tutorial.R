## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  fig.width=7.5, 
  fig.height=5, 
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
options(scipen = 99)
options(digits = 2)
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
matching_output = match_2C(Z = Z, X = X, 
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

