## =========================================================================

## Application of the Rank Ordered Logit (ROL) model

## Data generation for the Monte Carlo simulation in:
## Fok, Paap & van Dijk (2010), "A Rank-Ordered Logit Model with
## Unobserved Heterogeneity in Ranking Capabilities", Section 3.

## =========================================================================

rm(list = ls())
set.seed(123)  

## -------------------------------------------------------------------------
## Fixed settings of the simulation
## -------------------------------------------------------------------------

J <- 4      # number of alternatives
N <- 1000   # number of individuals (as in Table I of the paper)

## True beta parameters (Table I, "True" column).
## Category J = 4 is the base category, so beta_{.,4} = 0.
beta0 <- c(1.00,  0.25, -0.25, 0)   # intercepts,      j = 1,2,3,4
beta1 <- c(0.75, -0.50,  1.00, 0)   # coefficient x1,  j = 1,2,3,4
beta2 <- c(-0.30, 0.45,  0.80, 0)   # coefficient x2,  j = 1,2,3,4

## Segment (ranking-capability) probabilities.
## Segment k = 0,...,J-1 : individual can correctly rank the k
## most-preferred items; the remaining J-k items are ranked randomly.
## p0 = 0 : nobody is completely unable to rank (as chosen in the paper).
p_seg <- c(p0 = 0.0, p1 = 0.3, p2 = 0.3, p3 = 0.4)
stopifnot(abs(sum(p_seg) - 1) < 1e-10)

## -------------------------------------------------------------------------
## Function to generate one Monte Carlo data set
## -------------------------------------------------------------------------

generate_lcrol_data <- function(N, J, beta0, beta1, beta2, p_seg) {
  
  ## --- explanatory variables ---
  x1 <- rnorm(N, mean = 0, sd = 1)                 # continuous covariate
  x2 <- rbinom(N, size = 1, prob = 0.5)            # 0/1 dummy
  
  ## --- latent utilities U_ij = beta0_j + x1_i*beta1_j + x2_i*beta2_j + eps_ij ---
  ## eps_ij follows a type-I extreme value (standard Gumbel) distribution
  eps <- matrix(-log(-log(runif(N * J))), nrow = N, ncol = J)
  
  U <- matrix(NA_real_, nrow = N, ncol = J)
  for (j in 1:J) {
    U[, j] <- beta0[j] + x1 * beta1[j] + x2 * beta2[j] + eps[, j]
  }
  
  ## --- true ranking implied by the utilities (most preferred first) ---
  true_order <- t(apply(U, 1, function(u) order(u, decreasing = TRUE)))
  colnames(true_order) <- paste0("true_rank_pos", 1:J)
  
  ## --- assign each individual to a ranking-capability segment ---
  segment_labels <- 0:(J - 1)
  segment <- sample(segment_labels, size = N, replace = TRUE, prob = p_seg)
  
  ## --- construct the OBSERVED ranking for each individual ---
  ## The first k positions (k = segment) equal the true utility order;
  ## the remaining J-k positions are a random permutation of the
  ## remaining items (the respondent cannot rank these correctly).
  observed_order <- matrix(NA_integer_, nrow = N, ncol = J)
  for (i in 1:N) {
    k <- segment[i]
    ord_i <- true_order[i, ]
    if (k > 0) {
      observed_order[i, 1:k] <- ord_i[1:k]
    }
    if (k < J) {
      remaining_items <- ord_i[(k + 1):J]
      if (length(remaining_items) > 1) {
        remaining_items <- sample(remaining_items)  # random ordering
      }
      observed_order[i, (k + 1):J] <- remaining_items
    }
  }
  colnames(observed_order) <- paste0("obs_rank_pos", 1:J)
  
  ## --- also provide the "rank given to each item" representation ---
  ## y_ij = rank that individual i gives to item j (1 = most preferred)
  y <- matrix(NA_integer_, nrow = N, ncol = J)
  for (i in 1:N) {
    for (pos in 1:J) {
      item <- observed_order[i, pos]
      y[i, item] <- pos
    }
  }
  colnames(y) <- paste0("rank_item", 1:J)
  
  data.frame(
    id = 1:N,
    x1 = x1,
    x2 = x2,
    segment = segment,
    U,                     # latent utilities (usually unobserved in practice)
    true_order,            # true preference order (for checking purposes)
    observed_order,        # observed order (survey response, item numbers)
    y                      # observed rank given to each item
  )
}

## -------------------------------------------------------------------------
## Generate one Monte Carlo replication
## -------------------------------------------------------------------------

sim_data <- generate_lcrol_data(N, J, beta0, beta1, beta2, p_seg)

## =========================================================================
## Function: rol_probability
##
## Computes P(r_i; beta) for every individual i, where r_i = (r_i1,...,r_iJ)
## is the vector giving the item number that received rank j (rank 1 =
## most preferred), as in equation (2) of Fok, Paap & van Dijk (2010):
##
##   P(r_i; beta) = prod_{j=1}^{J-1} exp(V_i,r_ij) / sum_{l=j}^{J} exp(V_i,r_il)
##
## Linear index:  V_ij = beta0_j + x_i1*beta1_j + x_i2*beta2_j,  j = 1,...,J
## with beta_J = (0,0,0) fixed for identification (base category).
##
## Arguments:
##   data     : data frame containing columns x1, x2, and obs_rank_pos1,...,
##              obs_rank_posJ (the observed ranking, item numbers per rank
##              position), exactly as produced by generate_lcrol_data().
##   beta_vec : numeric vector of length 3*(J-1), ordered as
##              (beta0_1, beta1_1, beta2_1, beta0_2, beta1_2, beta2_2, ...,
##               beta0_{J-1}, beta1_{J-1}, beta2_{J-1})
##              i.e. parameters for alternatives 1,...,J-1 stacked in order;
##              alternative J is the base category (beta_J = 0) and must
##              NOT be included in beta_vec.
##   J        : number of alternatives. If NULL (default), it is inferred
##              from the number of obs_rank_pos* columns in data.
##
## Returns:
##   A numeric vector of length N with P(r_i; beta) for each individual.
## =========================================================================

rol_probability <- function(data, beta_vec, J = NULL) {
  
  ## --- infer J if not supplied ---
  if (is.null(J)) {
    J <- sum(grepl("^obs_rank_pos", names(data)))
  }
  
  N <- nrow(data)
  
  ## --- sanity check on parameter vector length ---
  if (length(beta_vec) != 3 * (J - 1)) {
    stop(sprintf("beta_vec must have length 3*(J-1) = %d, got %d",
                 3 * (J - 1), length(beta_vec)))
  }
  
  ## --- build (J x 3) matrix of parameters, adding the base category (0,0,0) ---
  beta_mat <- matrix(beta_vec, nrow = J - 1, ncol = 3, byrow = TRUE)
  beta_mat <- rbind(beta_mat, c(0, 0, 0))   # alternative J = base category
  colnames(beta_mat) <- c("beta0", "beta1", "beta2")
  
  ## --- design matrix (N x 3): intercept, x1, x2 ---
  X <- cbind(1, data$x1, data$x2)
  
  ## --- linear index V (N x J): V[,j] = beta0_j + x1*beta1_j + x2*beta2_j ---
  V <- X %*% t(beta_mat)
  
  ## --- extract the observed ranking (item number per rank position) ---
  obs_cols  <- paste0("obs_rank_pos", 1:J)
  order_mat <- as.matrix(data[, obs_cols])
  
  ## --- reorder V for each individual according to their reported ranking:
  ##     V_ordered[i, j] = V_i,r_ij  (utility of the item ranked j-th by i) ---
  V_ordered <- matrix(NA_real_, nrow = N, ncol = J)
  for (i in 1:N) {
    V_ordered[i, ] <- V[i, order_mat[i, ]]
  }
  
  expV <- exp(V_ordered)
  
  ## --- denom[,j] = sum_{l=j}^{J} exp(V_ordered[,l])  (reverse cumulative sum) ---
  denom <- t(apply(expV, 1, function(row) rev(cumsum(rev(row)))))
  
  ## --- P(r_i) = prod_{j=1}^{J-1} expV[,j] / denom[,j] ---
  ratio <- expV[, 1:(J - 1), drop = FALSE] / denom[, 1:(J - 1), drop = FALSE]
  P <- apply(ratio, 1, prod)
  
  return(P)
}

## -------------------------------------------------------------------------
## Quick check using the true parameter values from Table I
## (run after sourcing lcrol_simulation_data.R so that sim_data exists)
## -------------------------------------------------------------------------

true_beta_vec <- c(1.00, 0.75, -0.30,    # beta_1
                   0.25, -0.50, 0.45,   # beta_2
                   -0.25, 1.00, 0.80)   # beta_3

P_true <- rol_probability(sim_data, true_beta_vec)
summary(P_true) ; length(P_true)


## =========================================================================
## Function: rol_loglik
##
## Computes the (full-sample) ROL log-likelihood, eq. in the "Rank Ordered
## Logit from Scratch" note:
##
##   ell(beta) = sum_i sum_{j=1}^{J-1} x'_{i,r_ij} beta
##             - sum_i sum_{j=1}^{J-1} log( sum_{l=j}^{J} exp(x'_{i,r_il} beta) )
##
## which is equivalent to sum_i log P(r_i; beta), with P(r_i; beta) as in
## rol_probability(). Working directly in log-space (rather than taking
## log(prod(...))) avoids unnecessary underflow when multiplying many
## probabilities together.
##
## Arguments:
##   data     : data frame containing columns x1, x2, and obs_rank_pos1,...,
##              obs_rank_posJ (item number receiving rank j), as produced by
##              generate_lcrol_data().
##   beta_vec : numeric vector of length 3*(J-1), ordered as
##              (beta0_1, beta1_1, beta2_1, beta0_2, beta1_2, beta2_2, ...,
##               beta0_{J-1}, beta1_{J-1}, beta2_{J-1}).
##              Alternative J is the base category (beta_J = 0) and must
##              NOT be included in beta_vec.
##   J        : number of alternatives. If NULL (default), inferred from
##              the number of obs_rank_pos* columns in data.
##
## Returns:
##   A single numeric value: the log-likelihood evaluated at beta_vec,
##   summed over all N individuals in data.
## =========================================================================

rol_loglik <- function(data, beta_vec, J = NULL) {
  
  ## --- infer J if not supplied ---
  if (is.null(J)) {
    J <- sum(grepl("^obs_rank_pos", names(data)))
  }
  
  N <- nrow(data)
  
  ## --- sanity check on parameter vector length ---
  if (length(beta_vec) != 3 * (J - 1)) {
    stop(sprintf("beta_vec must have length 3*(J-1) = %d, got %d",
                 3 * (J - 1), length(beta_vec)))
  }
  
  ## --- build (J x 3) matrix of parameters, adding the base category (0,0,0) ---
  beta_mat <- matrix(beta_vec, nrow = J - 1, ncol = 3, byrow = TRUE)
  beta_mat <- rbind(beta_mat, c(0, 0, 0))   # alternative J = base category
  
  ## --- design matrix (N x 3): intercept, x1, x2 ---
  X <- cbind(1, data$x1, data$x2)
  
  ## --- linear index V (N x J): V[,j] = beta0_j + x1*beta1_j + x2*beta2_j ---
  V <- X %*% t(beta_mat)
  
  ## --- extract the observed ranking (item number per rank position) ---
  obs_cols  <- paste0("obs_rank_pos", 1:J)
  order_mat <- as.matrix(data[, obs_cols])
  
  ## --- reorder V for each individual according to their reported ranking:
  ##     V_ordered[i, j] = x'_{i,r_ij} beta ---
  V_ordered <- matrix(NA_real_, nrow = N, ncol = J)
  for (i in 1:N) {
    V_ordered[i, ] <- V[i, order_mat[i, ]]
  }
  
  expV <- exp(V_ordered)
  
  ## --- denom[,j] = sum_{l=j}^{J} exp(V_ordered[,l])  (reverse cumulative sum) ---
  denom <- t(apply(expV, 1, function(row) rev(cumsum(rev(row)))))
  
  ## --- ell = sum_{i,j<=J-1} V_ordered[i,j]  -  sum_{i,j<=J-1} log(denom[i,j]) ---
  term1 <- sum(V_ordered[, 1:(J - 1), drop = FALSE])
  term2 <- sum(log(denom[, 1:(J - 1), drop = FALSE]))
  
  ell <- term1 - term2
  return(ell)
}

## -------------------------------------------------------------------------
## Quick checks (run after sourcing lcrol_simulation_data.R and
## rol_probability.R, so that sim_data and rol_probability() exist)
## -------------------------------------------------------------------------
ll_true <- rol_loglik(sim_data, true_beta_vec) ; ll_true

# cross-check: should be numerically identical (up to float precision)
#  to summing the log of rol_probability()
all.equal(ll_true, sum(log(rol_probability(sim_data, true_beta_vec)))) # OK
# log-likelihood at beta = 0 should be worse (lower) than at the truth
rol_loglik(sim_data, rep(0, 9)) #OK

## =========================================================================
## Function: rol_gradient
##
## Computes the analytical gradient of the ROL log-likelihood:
##
##   grad(beta) = sum_i sum_{j=1}^{J-1} ( x_{i,r_ij} - sum_{l=j}^{J} P_{i,r_il|j} x_{i,r_il} )
##
##   where P_{i,r_il|j} = exp(x'_{i,r_il} beta) / sum_{k=j}^{J} exp(x'_{i,r_ik} beta)
##
## Here x_{ij} is the K-vector of covariates for alternative j, built with
## the usual alternative-specific dummy-interaction scheme (as in the
## drugChoice / ModeCanada example): for j = 1,...,J-1, x_{ij} places
## (1, x_i1, x_i2) in the block of coordinates belonging to alternative j
## and zero elsewhere; for the base alternative J, x_{iJ} = 0 (consistent
## with beta_J = 0 fixed, so V_iJ = 0 regardless of beta).
##
## Arguments:
##   data     : data frame containing columns x1, x2, and obs_rank_pos1,...,
##              obs_rank_posJ (item number receiving rank j), as produced by
##              generate_lcrol_data().
##   beta_vec : numeric vector of length K = 3*(J-1), ordered as
##              (beta0_1, beta1_1, beta2_1, beta0_2, beta1_2, beta2_2, ...,
##               beta0_{J-1}, beta1_{J-1}, beta2_{J-1}).
##   J        : number of alternatives. If NULL (default), inferred from
##              the number of obs_rank_pos* columns in data.
##
## Returns:
##   A numeric vector of length K = 3*(J-1): the gradient evaluated at
##   beta_vec, summed over all N individuals in data.
## =========================================================================

rol_gradient <- function(data, beta_vec, J = NULL) {
  
  ## --- infer J if not supplied ---
  if (is.null(J)) {
    J <- sum(grepl("^obs_rank_pos", names(data)))
  }
  
  N <- nrow(data)
  K <- 3 * (J - 1)
  
  if (length(beta_vec) != K) {
    stop(sprintf("beta_vec must have length 3*(J-1) = %d, got %d", K, length(beta_vec)))
  }
  
  ## --- build (J x 3) matrix of parameters, adding the base category (0,0,0) ---
  beta_mat <- matrix(beta_vec, nrow = J - 1, ncol = 3, byrow = TRUE)
  beta_mat <- rbind(beta_mat, c(0, 0, 0))   # alternative J = base category
  
  ## --- base covariates (N x 3): intercept, x1, x2 ---
  X <- cbind(1, data$x1, data$x2)
  
  ## --- linear index V (N x J) ---
  V <- X %*% t(beta_mat)
  
  ## --- observed ranking: item number receiving rank j ---
  obs_cols  <- paste0("obs_rank_pos", 1:J)
  order_mat <- as.matrix(data[, obs_cols])
  
  ## --- helper: K-vector x_{i,alt} for a given alternative index ---
  get_x_alt <- function(xi, alt) {
    v <- numeric(K)
    if (alt <= (J - 1)) {
      idx <- ((alt - 1) * 3 + 1):((alt - 1) * 3 + 3)
      v[idx] <- xi
    }
    v   # zero vector if alt == J (base category)
  }
  
  grad <- numeric(K)
  
  for (i in 1:N) {
    
    r_i  <- order_mat[i, ]     # r_i[j] = alternative ranked j-th by individual i
    xi   <- X[i, ]             # (1, x1_i, x2_i)
    V_i  <- V[i, ]             # utilities for alternatives 1,...,J
    
    V_ordered    <- V_i[r_i]   # V reordered according to the ranking
    expV_ordered <- exp(V_ordered)
    
    for (j in 1:(J - 1)) {
      
      idx_remaining <- j:J                                  # positions l = j,...,J
      denom_j       <- sum(expV_ordered[idx_remaining])
      P_l_given_j   <- expV_ordered[idx_remaining] / denom_j # conditional probs
      
      x_rij <- get_x_alt(xi, r_i[j])   # covariate vector of the item ranked j-th
      
      ## weighted sum over remaining alternatives: sum_l P_{l|j} * x_{i,r_il}
      x_weighted_sum <- numeric(K)
      for (m in seq_along(idx_remaining)) {
        l     <- idx_remaining[m]
        alt_l <- r_i[l]
        x_weighted_sum <- x_weighted_sum + P_l_given_j[m] * get_x_alt(xi, alt_l)
      }
      
      grad <- grad + (x_rij - x_weighted_sum)
    }
  }
  
  return(grad)
}

## -------------------------------------------------------------------------
## Quick check: the analytical gradient should match a numerical gradient
## of rol_loglik() (run after sourcing lcrol_simulation_data.R, rol_loglik.R)
## -------------------------------------------------------------------------
library(numDeriv)
g_analytical <- rol_gradient(sim_data, true_beta_vec)
g_analytical
g_numerical  <- grad(function(b) rol_loglik(sim_data, b), true_beta_vec)
cbind(analytical = g_analytical, numerical = g_numerical)
round(max(abs(g_analytical - g_numerical)),5)   # close to 0, OK. 


## =========================================================================
## Function: rol_hessian
##
## Computes the analytical Hessian H(beta) = -grad_{beta beta'} ell(beta)
## of the ROL log-likelihood:
##
##   H(beta) = sum_i sum_{j=1}^{J-1} sum_{l=j}^{J} P_{i,r_il|j} *
##                 ( x_{i,r_il} - sum_{k=j}^{J} x_{i,r_ik} P_{i,r_ik|j} ) x'_{i,r_il}
##
## where P_{i,r_il|j} = exp(x'_{i,r_il} beta) / sum_{k=j}^{J} exp(x'_{i,r_ik} beta),
## and x_{ij} is the K-vector of covariates for alternative j, built with the
## same alternative-specific dummy-interaction scheme used in rol_gradient():
## for j = 1,...,J-1, x_{ij} places (1, x_i1, x_i2) in the block of
## coordinates belonging to alternative j and zero elsewhere; for the base
## alternative J, x_{iJ} = 0.
##
## H(beta) is the negative Hessian of the log-likelihood; since ell(beta) is
## globally concave, H(beta) is symmetric positive semi-definite.
##
## Arguments:
##   data     : data frame containing columns x1, x2, and obs_rank_pos1,...,
##              obs_rank_posJ (item number receiving rank j), as produced by
##              generate_lcrol_data().
##   beta_vec : numeric vector of length K = 3*(J-1), ordered as
##              (beta0_1, beta1_1, beta2_1, beta0_2, beta1_2, beta2_2, ...,
##               beta0_{J-1}, beta1_{J-1}, beta2_{J-1}).
##   J        : number of alternatives. If NULL (default), inferred from
##              the number of obs_rank_pos* columns in data.
##
## Returns:
##   A K x K numeric matrix (K = 3*(J-1)): the negative Hessian evaluated
##   at beta_vec, summed over all N individuals in data.
## =========================================================================

rol_hessian <- function(data, beta_vec, J = NULL) {
  
  ## --- infer J if not supplied ---
  if (is.null(J)) {
    J <- sum(grepl("^obs_rank_pos", names(data)))
  }
  
  N <- nrow(data)
  K <- 3 * (J - 1)
  
  if (length(beta_vec) != K) {
    stop(sprintf("beta_vec must have length 3*(J-1) = %d, got %d", K, length(beta_vec)))
  }
  
  ## --- build (J x 3) matrix of parameters, adding the base category (0,0,0) ---
  beta_mat <- matrix(beta_vec, nrow = J - 1, ncol = 3, byrow = TRUE)
  beta_mat <- rbind(beta_mat, c(0, 0, 0))   # alternative J = base category
  
  ## --- base covariates (N x 3): intercept, x1, x2 ---
  X <- cbind(1, data$x1, data$x2)
  
  ## --- linear index V (N x J) ---
  V <- X %*% t(beta_mat)
  
  ## --- observed ranking: item number receiving rank j ---
  obs_cols  <- paste0("obs_rank_pos", 1:J)
  order_mat <- as.matrix(data[, obs_cols])
  
  ## --- helper: K-vector x_{i,alt} for a given alternative index ---
  get_x_alt <- function(xi, alt) {
    v <- numeric(K)
    if (alt <= (J - 1)) {
      idx <- ((alt - 1) * 3 + 1):((alt - 1) * 3 + 3)
      v[idx] <- xi
    }
    v   # zero vector if alt == J (base category)
  }
  
  H <- matrix(0, K, K)
  
  for (i in 1:N) {
    
    r_i <- order_mat[i, ]      # r_i[j] = alternative ranked j-th by individual i
    xi  <- X[i, ]              # (1, x1_i, x2_i)
    V_i <- V[i, ]              # utilities for alternatives 1,...,J
    
    V_ordered    <- V_i[r_i]   # V reordered according to the ranking
    expV_ordered <- exp(V_ordered)
    
    ## covariate vectors for each rank position, precomputed once per i
    x_ordered <- lapply(r_i, function(alt) get_x_alt(xi, alt))
    
    for (j in 1:(J - 1)) {
      
      idx_remaining <- j:J                                    # positions l = j,...,J
      denom_j       <- sum(expV_ordered[idx_remaining])
      P_l_given_j   <- expV_ordered[idx_remaining] / denom_j   # conditional probs
      
      ## weighted mean covariate vector at stage j: xbar_j = sum_k P_{k|j} x_k
      xbar_j <- numeric(K)
      for (m in seq_along(idx_remaining)) {
        xbar_j <- xbar_j + P_l_given_j[m] * x_ordered[[idx_remaining[m]]]
      }
      
      ## accumulate sum_l P_{l|j} (x_l - xbar_j) x_l'
      for (m in seq_along(idx_remaining)) {
        l   <- idx_remaining[m]
        x_l <- x_ordered[[l]]
        H <- H + P_l_given_j[m] * outer(x_l - xbar_j, x_l)
      }
    }
  }
  
  return(H)
}

## -------------------------------------------------------------------------
## Quick check: the analytical Hessian should match a numerical Hessian
## of rol_loglik() (run after sourcing lcrol_simulation_data.R, rol_loglik.R)
## -------------------------------------------------------------------------

H_analytical <- rol_hessian(sim_data, true_beta_vec)
H_numerical  <- -hessian(function(b) rol_loglik(sim_data, b), true_beta_vec)
#
max(abs(H_analytical - H_numerical)) # OK
isSymmetric(H_analytical)              # OK


## =========================================================================
## Function: rol_newton_raphson
##
## Newton-Raphson maximum likelihood estimation of the ROL model, using
## rol_loglik(), rol_gradient() and rol_hessian() defined earlier.
##
## IMPORTANT SIGN NOTE:
##   H(beta) = -grad_{beta beta'} ell(beta)  (as returned by rol_hessian)
## is the NEGATIVE Hessian of the log-likelihood (positive semi-definite,
## since ell is globally concave). The Newton step that maximizes ell is
## therefore
##       beta_new = beta - [grad_{beta beta'} ell(beta)]^{-1} g(beta)
##                = beta + [H(beta)]^{-1} g(beta)
## i.e. a PLUS, not a minus, when using H as defined here.
##
## Arguments:
##   data  : data frame with x1, x2, obs_rank_pos1,...,obs_rank_posJ
##   start : starting values for beta (length K = 3*(J-1))
##   J     : number of alternatives (inferred from data if NULL)
##   tol   : convergence tolerance on max absolute parameter change
##   maxit : maximum number of iterations
##
## Returns a list with:
##   beta       : estimated parameter vector (named)
##   loglik     : log-likelihood at the estimate
##   Hessian    : H(beta_hat), i.e. -Hessian of ell (used for std. errors)
##   iterations : number of iterations used
## =========================================================================

rol_newton_raphson <- function(data, start, J = NULL, tol = 1e-10, maxit = 100) {
  
  if (is.null(J)) {
    J <- sum(grepl("^obs_rank_pos", names(data)))
  }
  
  beta <- start
  
  for (it in 1:maxit) {
    
    g <- rol_gradient(data, beta, J)
    H <- rol_hessian(data, beta, J)
    
    step     <- solve(H, g)     # H^{-1} g
    beta_new <- beta + step     # PLUS, since H = -Hessian(ell)
    
    if (max(abs(beta_new - beta)) < tol) {
      beta <- beta_new
      cat("Converged in", it, "iterations\n")
      break
    }
    beta <- beta_new
  }
  
  ll <- rol_loglik(data, beta, J)
  
  list(beta = beta, loglik = ll, Hessian = H, iterations = it)
}

## -------------------------------------------------------------------------
## Example usage on sim_data
## (run after sourcing lcrol_simulation_data.R, rol_loglik.R,
##  rol_gradient.R, rol_hessian.R, so that sim_data and the three
##  functions above are available)
## -------------------------------------------------------------------------
J ; K <- 3 * (J - 1)

# ## descriptive names for the parameter vector: (intercept, x1, x2) per
# ## non-base alternative, matching the order used throughout
par_names <- paste0(rep(c("intercept_alt", "x1_alt", "x2_alt"), J - 1),
                    rep(1:(J - 1), each = 3))

start <- rep(0, K) # starting values for the search
names(start) <- par_names # associate the name of each parameter

fit_nr <- rol_newton_raphson(sim_data, start) # start the search
names(fit_nr$beta) <- par_names
fit_nr

se_nr <- sqrt(diag(solve(fit_nr$Hessian)))
names(se_nr) <- par_names
se_nr

# ## quick comparison against the true beta values used to simulate the data
cbind(estimate = fit_nr$beta, true = true_beta_vec, se = se_nr) 


## =========================================================================
## Monte Carlo study: 10,000 replications
##
## For each replication we (i) simulate a fresh data set with
## generate_lcrol_data(), (ii) estimate beta via rol_newton_raphson()
## (starting at 0 each time), and (iii) store the resulting estimates.
## At the end we report, for each parameter, the mean and standard
## deviation of the estimates across replications, together with the
## true value used to generate the data -- mirroring Table I of the paper.
##
## NOTE ON RUNTIME: rol_gradient() and rol_hessian() use explicit loops
## over the N = 1000 individuals (and, within each, over the J-1 ranking
## stages), so each Newton-Raphson iteration is not fully vectorized.
## With ~5-10 NR iterations per replication and 10,000 replications, this
## can take a while to run. 
## =========================================================================

## -------------------------------------------------------------------------
## Redefine rol_newton_raphson with a 'verbose' switch, so we don't print
## "Converged in X iterations" 10,000 times. Behavior is unchanged when
## verbose = TRUE (the default), so any earlier calls in the script are
## unaffected.
## -------------------------------------------------------------------------

rol_newton_raphson <- function(data, start, J = NULL, tol = 1e-10,
                               maxit = 100, verbose = TRUE) {
  
  if (is.null(J)) {
    J <- sum(grepl("^obs_rank_pos", names(data)))
  }
  
  beta <- start
  
  for (it in 1:maxit) {
    
    g <- rol_gradient(data, beta, J)
    H <- rol_hessian(data, beta, J)
    
    step     <- solve(H, g)     # H^{-1} g
    beta_new <- beta + step     # PLUS, since H = -Hessian(ell)
    
    if (max(abs(beta_new - beta)) < tol) {
      beta <- beta_new
      if (verbose) cat("Converged in", it, "iterations\n")
      break
    }
    beta <- beta_new
  }
  
  ll <- rol_loglik(data, beta, J)
  
  list(beta = beta, loglik = ll, Hessian = H, iterations = it)
}

## -------------------------------------------------------------------------
## Monte Carlo loop
## -------------------------------------------------------------------------

n_reps <- 10000

mc_results <- matrix(NA_real_, nrow = n_reps, ncol = K)
colnames(mc_results) <- par_names

set.seed(1234)  # separate seed dedicated to the Monte Carlo study

t0 <- proc.time()

for (r in 1:n_reps) {
  
  dat_r <- generate_lcrol_data(N, J, beta0, beta1, beta2, p_seg)
  
  fit_r <- tryCatch(
    rol_newton_raphson(dat_r, start, verbose = FALSE),
    error = function(e) NULL
  )
  
  if (!is.null(fit_r)) {
    mc_results[r, ] <- fit_r$beta
  }
  
  if (r %% 500 == 0) {
    cat("Replication", r, "of", n_reps,
        "-- elapsed:", round((proc.time() - t0)["elapsed"], 1), "sec\n")
  }
}

cat("Total time:", round((proc.time() - t0)["elapsed"] / 60, 1), "minutes\n")

n_failed <- sum(is.na(mc_results[, 1]))
cat("Number of failed replications (dropped):", n_failed, "\n")

## -------------------------------------------------------------------------
## Summary: mean and sd of the estimates across replications,
## compared to the true parameter values
## -------------------------------------------------------------------------

mc_mean <- colMeans(mc_results, na.rm = TRUE)
mc_sd   <- apply(mc_results, 2, sd, na.rm = TRUE)

mc_summary <- data.frame(
  parameter = par_names,
  true      = true_beta_vec,
  mc_mean   = mc_mean,
  mc_sd     = mc_sd,
  bias      = mc_mean - true_beta_vec
)

print(mc_summary, digits = 4)

save(mc_summary, file = "mc_summary_ROL")




