#' Calculate a covariance matrix
#' 
#' Calculate a covariance matrix of the design matrix \code{x} using one of the methods 
#' \code{complete_observations}, \code{pairwise_covariance}, \code{loh_wainwright_bias_correction} or \code{average_imputation}.
#' 
#' @param x A n x p design matrix
#' @param NA_method The method that should be applied to estimate the covariance structure.
#' @param min_points Minimal number of available observations of a variable necessary such that an estimation of the variance of that 
#' variable is reliable. An integer not smaller than 2.
#' @return A p x p matrix with an estimate for the covariance structure of \code{x}. If less than \code{min_points} observations are
#' available for some variable, the corresponding rows and columns are NA.
#' @importFrom Matrix nearPD
#' @importFrom stats cov
#' @importFrom stats complete.cases
#' @export
get_cov_mat <- function(x, NA_method = c("complete_observations", "pairwise_covariance", 
                                         "loh_wainwright_bias_correction", "average_imputation"), min_points = 2) {
  stopifnot(min_points >= 2)
  NA_mth <- match.arg(NA_method)
  p <- ncol(x)
  cov_mat <- matrix(NA, ncol = p, nrow = p)
  
  if (NA_mth == "complete_observations") {
    
    n_obs <- sum(stats::complete.cases(x))
    
    if (n_obs < 2) {
      return(list(mat = cov_mat, inds = rep(F, p), n_eff_obs = 0))
    } else {
      return(list(mat = (n_obs - 1)/n_obs * cov(x, use = "na.or.complete"), 
                  inds = rep(T, ncol(x)), n_eff_obs = n_obs))
    }
    
    
  } else {
    
    # calculate available entries for each column
    available_obs <- apply(!is.na(x), 2, sum)
    
    # we can only estimate variance for columns with at least 2 available
    # entries
    n <- nrow(x)
    inds <- (available_obs >= max(2, min_points))
    
    # if no variable has enough observations, return NA
    if (!any(inds)) {
      return(list(mat = cov_mat, inds = rep(F, p), n_eff_obs = 0))
    }
    
    n_eff_obs <- sum(available_obs)/sum(inds)
    
    # pairwise covariance method. Calculates pairwise covariances and
    # impute leftover covariances with 0.
    if (NA_mth == "pairwise_covariance") {
      
      temp <- stats::cov(x[, inds, drop = F], use = "pairwise")
      temp[is.na(temp)] <- 0
      
      # As the result might not be positive semi-definite, which is required
      # for e.g. glasso, project onto closest (wrt Frobenius) positive
      # definite matrix
      cov_mat[inds, inds] <- as.matrix(Matrix::nearPD(temp)$mat)
      
      
      # Loh-Wainwright bias correction method.  Tries to remove bias of
      # (co-)variance estimation if missing values are present.
    } else if (NA_mth == "loh_wainwright_bias_correction") {
      
      # first center the observations and impute the missing values with the
      # mean of the corresponding column
      z <- scale(x[, inds, drop = F], T, F)
      z[is.na(z)] <- 0
      
      # estimate local missingness probability per predictor
      miss_frac <- 1 - available_obs[inds, drop = F]/nrow(x)
      stopifnot(length(miss_frac) == ncol(z))
      
      # Initiate helper matrix for LW bias correction
      M <- matrix(rep((1 - miss_frac), ncol(z)), ncol = ncol(z))
      M <- t(M) * M
      diag(M) <- 1 - miss_frac
      
      cov_mat[inds, inds] <- as.matrix(Matrix::nearPD(stats::cov(z) * (n - 
                                                                  1)/n/M)$mat)
      
      # to avoid double debiasing, use n_eff_obs = n_obs = nrow(x) n_eff_obs
      # <- nrow(z)
      
    } else if (NA_mth == "average_imputation") {
      # center the observations and impute the missing values with the mean
      # of the corresponding column
      z <- scale(x[, inds, drop = F], T, F)
      z[is.na(z)] <- 0
      
      cov_mat[inds, inds] <- stats::cov(z) * (n_eff_obs - 1)/n_eff_obs
    }
    
    list(mat = cov_mat, inds = inds, n_eff_obs = n_eff_obs)
  }
}


#' Negative loglikelihood of a multivariate normal
#' 
#' @param x A design matrix
#' @param mu Mean of the multivariate Gaussian distribution
#' @param cov_mat_inv precision matrix of the multivariate Gaussian distribution 
#' @return An array with the negative log-likelihoods of the observations in \code{x} given the Gaussian model
#' specified.
#' @export
loglikelihood <- function(x, mu, cov_mat_inv) {
  if (!is.matrix(x)) {
    x <- as.matrix(x)
    warning("non matrix input into loglikelihood")
  }
  
  n <- nrow(x)
  p <- ncol(x)
  loss <- rep(NA, n)
  
  # to avoid multiple inverseions of cov_mat[inds_cur, inds_cur], order
  # observations by missingness structure and handle in groups where
  # missingness structure is same
  ranks <- data.table::frank(data.table::data.table(is.na(x)))
  
  for (j in unique(ranks)) {
    
    # x_cur is a group of observations with the same missingness structure
    x_cur <- x[ranks == j, , drop = F]
    inds_cur <- !is.na(x_cur[1, ])
    
    # are there any missing values and at least one observed value?
    if (any(!inds_cur) & any(inds_cur)) {
      log_det <- -determinant(cov_mat_inv[inds_cur, inds_cur, drop = F], 
                              logarithm = TRUE)$modulus
      v <- t(x_cur[, inds_cur, drop = F]) - mu[inds_cur, drop = F]
      # distance <- diag(t(v) %*% solve(cov_mat[inds_cur, inds_cur, drop =
      # F], v))
      distance <- diag(t(v) %*% cov_mat_inv[inds_cur, inds_cur, drop = F] %*% 
                         v)
      loss[ranks == j] <- distance + log_det + sum(inds_cur) * log(2 * pi)
      # are there no missing values?
    } else if (!any(!inds_cur)) {
      log_det <- -determinant(cov_mat_inv, logarithm = TRUE)$modulus
      v <- t(x_cur) - mu
      distance <- diag(t(v) %*% cov_mat_inv %*% v)
      loss[ranks == j] <- distance + log_det + sum(inds_cur) * log(2 * 
                                                                     pi)
      # are there no observed values?
    } else if (all(!inds_cur)) {
      loss[ranks == j] <- 0
    }
  }
  # return loss vector
  loss/2
}

#' Get a glasso fit
#' 
#' @param x A design matrix
#' @param lambda Regularisation parameter for the glasso function.
#' @param control An object of class \code{hdcd_control} generated by \link{hdcd_control}.
#' @return A list with values \code{cov_mat}, the original estimate of the covariance structure,
#' \code{inds} an array indicating the variables for which a variance could be estimated, 
#' the precision matrix estimated with the glasso \code{wi} and its inverse \code{w}. 
#' @importFrom glasso glasso
#' @export
get_glasso_fit <- function(x, lambda, control) {
  
  penalize_diagonal <- control$glasso_penalize_diagonal
  standardize <- control$glasso_standardize
  threshold <- 1e-4
  min_points <- 2
  NA_method <- control$glasso_NA_method
  n_obs <- nrow(x)
  
  stopifnot(!any(is.null(penalize_diagonal), is.null(standardize), is.null(threshold), 
                 is.null(min_points), is.null(NA_method), is.null(n_obs)))
  
  # Get (estimate of) covariance matrix
  cov_mat_output <- get_cov_mat(x, NA_method, min_points = min_points)
  
  # Prepare list to save results in
  out <- list()
  out$cov_mat <- cov_mat_output$mat
  out$inds <- cov_mat_output$inds
  out$w <- matrix(NA, nrow = ncol(x), ncol = ncol(x))
  out$wi <- matrix(NA, nrow = ncol(x), ncol = ncol(x))
  
  p_cur <- sum(out$inds)
  if (p_cur == 0) {
    return(out)
  }
  
  # lambda should be rescaled by sqrt of length of segment corresponding
  # to asymptotic theory
  obs_share <- cov_mat_output$n_eff_obs/n_obs
  
  if (standardize) {
    rho <- lambda/sqrt(obs_share) * diag(out$cov_mat[out$inds, out$inds, 
                                                     drop = F])
  } else {
    rho <- lambda/sqrt(obs_share)
  }
  
  glasso_output <- glasso::glasso(out$cov_mat[out$inds, out$inds, drop = F], 
                                  rho = rho, penalize.diagonal = penalize_diagonal, thr = threshold)
  
  out$w[out$inds, out$inds] <- glasso_output$w
  out$wi[out$inds, out$inds] <- glasso_output$wi
  out$mu <- colMeans(x, na.rm = T)
  
  # rescaling of loglikelihood returned from glasso
  if (!standardize) {
    out$loglik <- ((-2/p_cur) * glasso_output$l - sum(abs(lambda/sqrt(obs_share) * 
                                                            glasso_output$wi))) * obs_share
  } else {
    out$loglik <- ((-2/p_cur) * glasso_output$l - sum(abs(lambda/sqrt(obs_share) * 
                                                            sqrt((diag(out$cov_mat[out$inds, out$inds, drop = F]) %*% t(diag(out$cov_mat[out$inds, 
                                                                                                                                         out$inds, drop = F])))) * glasso_output$wi))) * obs_share
  }
  
  out
}
