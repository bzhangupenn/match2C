#' Assign subgroup indicators
#'
#' Partitions x into beta strata using:
#' - S_1:     (1/3, 1]
#' - S_k:     (1/(k+2), 1/(k+1)] for k = 2, ..., beta-1
#' - S_beta:  [0, 1/(beta + 1)]
#'
#' @param x Numeric vector to stratify. Values outside [0, 1] are returned as NA.
#' @param beta Integer (>= 2): number of subgroups.
#'
#' @return An integer vector of subgroup indicators (1..beta), with NA where x is NA
#'   or outside [0, 1].
#' @export
assign_int <- function(x, beta) {
  if (!is.numeric(x)) stop("`x` must be numeric.")
  if (length(beta) != 1L || !is.finite(beta) || beta < 2 || beta != as.integer(beta)) {
    stop("`beta` must be a single integer >= 2.")
  }
  beta <- as.integer(beta)

  if (any(is.na(x)) || any(x < 0 | x > 1)) {
    stop("`x` must be between 0 and 1 and non-NA.")
  }
  
  result <- rep(NA_integer_, length(x))
  
  for (i in seq_along(x)) {
    xi <- x[i]
    # S_1: (1/3, 1]
    if (xi > 1/3 && xi <= 1) {
      result[i] <- 1
    }
    
    assigned <- FALSE
    # S_k: (1/(k+2), 1/(k+1)] for k = 2 to beta - 1
    if (beta > 2){ 
      for (k in 2:(beta - 1)) {
        lower <- 1 / (k + 2)
        upper <- 1 / (k + 1)
        if (xi > lower && xi <= upper) {
          result[i] <- k
          assigned <- TRUE
          break
        }
      }
    }
    
    # S_beta: [0, 1/(beta + 1)]
    if (!assigned && xi >= 0 && xi <= 1 / (beta + 1)) {
      result[i] <- beta
    }
  }
  return(result)
}



#' Calculate feasible k for 1-to-k matching
#'
#' Computes the largest feasible matching ratio \code{k} between treated and
#' control units within the levels of a categorical covariate. The function
#' reports whether controls or treated units dominate, and whether roles
#' should be flipped to achieve 1-to-k matching.
#'
#' @param factorco A categorical variable used for fine balance.
#' @param Z A binary (0/1) treatment indicator vector of the same length as \code{factorco}.
#' @param verbose Logical; if \code{TRUE}, prints diagnostic messages. Default is \code{FALSE}.
#'
#' @return A list with components:
#' \item{k}{Maximum feasible matching ratio \code{k}, or \code{NULL} if infeasible.}
#' \item{role}{Dominance type: \code{"control_dominant"}, \code{"treated_dominant"},
#'   \code{"mixed"}, or \code{"infeasible"}.}
#' \item{flip}{Logical; \code{TRUE} if roles should be flipped, otherwise \code{FALSE}.
#'            \code{FALSE} means more controls}
#' \item{levels}{Character vector of levels with both treated and control units.}
#' 
#' @export
kcalculator <- function(factorco, Z, verbose = FALSE) {

if (!is.factor(factorco)) factorco <- factor(factorco)
levels_all <- levels(factorco)

# counts per level
n_t_all <- sapply(levels_all, function(l) sum(Z == 1 & factorco == l))
n_c_all <- sapply(levels_all, function(l) sum(Z == 0 & factorco == l))

control_dominant <- any(n_c_all > n_t_all)
treated_dominant <- any(n_t_all > n_c_all)  

if (control_dominant && treated_dominant) {
  if (verbose) message("Infeasible: mixed dominance across original levels.")
  return(list(k = NULL, role = "mixed", flip = FALSE, levels = character(0)))
}

# keep only levels with both treated and control
valid <- which(n_t_all > 0 & n_c_all > 0)
if (length(valid) == 0) {
  if (verbose) message("Infeasible: no levels with both treated and control.")
  return(list(k = NULL, role = "infeasible", flip = FALSE, levels = character(0)))
}

n_t <- n_t_all[valid]
n_c <- n_c_all[valid]

all_control_dominant <- all(n_c >= n_t)
all_treated_dominant <- all(n_t >= n_c)

if (all_control_dominant) {
  k_max <- min(floor(n_c / n_t))
  role  <- "control_dominant"
  flip  <- FALSE
} else if (all_treated_dominant) {
  # To perform k:1 with *more controls*, we suggest flipping roles.
  k_max <- min(floor(n_t / n_c))
  role  <- "treated_dominant"
  flip  <- TRUE
} else {
  if (verbose) message("Mixed dominance across levels; fine balance infeasible.")
  return(list(k = NULL, role = "mixed", flip = FALSE, levels = names(n_t)))
}

# guard: if k_max < 1 then effectively infeasible
if (is.na(k_max) || k_max < 1) {
  return(list(k = NULL, role = role, flip = flip, levels = names(n_t)))
}

return(list(k = as.integer(k_max), role = role, flip = flip, levels = names(n_t)))
}



#' simple covariate data preprocessing
#' @param X A n-by-p numerical data frame.
#' @param fine Column name of covariates used for fine balance
#' @param var_tol A non-negative numeric tolerance for variance filtering.
#' @param qr_tol A numeric tolerance for QR decomposition rank determination.
#' 
#' @return A data frame containing only the retained covariates 
#'   (non-constant and linearly independent columns).
mahapd <- function(X, fine, var_tol = 0, qr_tol = 1e-10) {
  proc <- setdiff(names(X), fine)
  keep <- character(0)
  
  if (length(proc) > 0L) {
    Xco <- as.matrix(X[proc])
    v <- apply(Xco, 2, var, na.rm = TRUE)
    keep_var <- which(is.finite(v) & v > var_tol)
    if (length(keep_var) > 0L) {
      Xco <- Xco[, keep_var, drop = FALSE]
      q <- qr(Xco, tol = qr_tol)  # keep linearly independent columns
      keep <- colnames(Xco)[q$pivot[seq_len(q$rank)]]
    }
  }
  X[, names(X) %in% c(fine, keep), drop = FALSE] 
}
