#' Closure generating function to calculate gains for splits using a classifier
#'
#' @inheritParams classifier_best_split_function
#' @inheritParams hdcd
#' @param start The start of the investigated interval
#' @param end The end of the investigated interval
#'
#' @return A closure with parameters \code{x}, \code{start} and \code{end}, that when evaluated
#' will itself return a closure with arguments \code{split_point} and \code{split_candidates}.
#' Evaluated this function returns the gain when splitting the segment \code{(start, end]} of \code{x}
#'  at \code{split_point}.
classifier_gain_function <-
  function(x, start, end, classifier, control) {
    function(split_point, split_candidates = split_point) {
      n1 <- split_point - start
      n2 <- end - split_point
      n <- n1 + n2
      stopifnot(n1 > 0 & n2 > 0)
      
      # prepare labels
      y <- rep.int(1:2, times = c(n1, n2))
      
      predictions <-
        classifier(x_train = x[(start + 1):end,],
                   y_train = y,
                   control = control)$train
      
      if (is.null(dim(predictions)) || dim(predictions) != c(n, 2)) {
        stop('classifier$train needs to return a n x k matrix with predictions for each class')
      }
      
      stopifnot(dim(predictions) == c(n, 2))
      
      likelihood_matrix <-
        classifier_loglikelihood(
          y_train = y,
          predictions = predictions,
          oob = control$classifier_oob
        )$likelihood_matrix
      
      likelihood_sum <- sum(likelihood_matrix[, 2])
      likelihood_difference <-
        likelihood_matrix[, 1] - likelihood_matrix[, 2]
      
      gain <- array(NA, nrow(x))
      gain[(start + 1):end] <-
        (likelihood_sum + cumsum(likelihood_difference)) / nrow(x)
      best_split <-
        split_candidates[which.max(gain[split_candidates])]
      max_gain <- gain[best_split]
      
      n_perm <- control$permutation_test_n
      
      set.seed(0)
      permutation_test <- array(NA, dim = n_perm)
      
      for (i in 1:n_perm) {
        sigma <- sample(1:nrow(x))
        sigma <- sigma[sigma >= (start + 1) & sigma <= end] - start
        permutation_test[i] <-
          (likelihood_sum + max(cumsum(likelihood_difference[sigma])[split_candidates - start])) / nrow(x)
      }
      
      list(
        gain = gain,
        best_split = best_split,
        max_gain = max_gain,
        permutation_test = permutation_test,
        pval = (1 + sum(permutation_test >= max_gain)) / (length(permutation_test) + 1)
      )
    }
  }
# }
# classifier_gain_function <- function(x, start, end, classifier, control){
#
#   function(split_point, split_candidates = split_point){
#
#     n1 <- split_point - start
#     n2 <- end - split_point
#     n <- n1 + n2
#     stopifnot(n1 > 0 & n2 > 0)
#
#     # prepare labels
#     y <- rep.int(1 : 2, times = c(n1, n2))
#
#     predictions <- classifier(x_train = x[(start + 1) : end, ], y_train = y, control = control)$train
#
#     if(is.null(dim(predictions)) || dim(predictions) != c(n, 2)){
#       stop('classifier$train needs to return a n x k matrix with predictions for each class')
#     }
#
#     stopifnot(dim(predictions) == c(n, 2))
#
#     likelihood_matrix <- classifier_loglikelihood(y_train = y, predictions = predictions, oob = control$classifier_oob)$likelihood_matrix
#
#     gain <- sum(likelihood_matrix[, 2]) + cumsum(likelihood_matrix[, 1] - likelihood_matrix[, 2])
#
#     gain[split_candidates - start]
#   }
#   #   gain[split_candidates - start]
#   #   if(evaluate_gain == 'single'){
#   #     sum(log_eps(n / (n1 - 1) * predictions[1 : n1])) + sum(log_eps(n / (n2 - 1) * predictions[(n1 + 1) : n])) / nrow(x)
#   #   } else if (evaluate_gain == 'all'){
#   #     gain <- rep(NA, nrow(x))
#   #     logits <- log_eps( n / (n1 - (y == 1)) * predictions) - log_eps( n / (n2 - (y == 2)) * (1 - predictions))
#   #     gain[(start + 1) : end] <- cumsum(logits) + sum(pmin(6, pmax(-6, log(1 - predictions)))) + n1 * log((n - 1) / n2) + n2 * log((n - 1) / (n2 - 1))
#   #     gain / nrow(x)
#   #   } else if (evaluate_gain == 'adjusted'){
#   #     gain <- rep(NA, nrow(x))
#   #     logits <- log_eps( n / (n1 - (y == 1)) * predictions) - log_eps( n / (n2 - (y == 2)) * (1 - predictions))
#   #     gain[(start + 1) : end] <- cumsum(logits) + sum(pmin(6, pmax(-6, log(1 - predictions)))) + n1 * log((n - 1) / n2) + n2 * log((n - 1) / (n2 - 1))
#   #     gain[(start + 1) : end] <- gain[(start + 1) : end] - c(seq(from = 1, to = 0, length.out = n1)^2 * gain[start + 1],  seq(from = 1 / n2, to = 1, length.out = n2)^2 * gain[end])
#   #     gain / nrow(x)
#   #   } else if (evaluate_gain == 'shift_in_mean'){
#   #     gain <- rep(NA, nrow(x))
#   #     logits <- log_eps( n / (n1 - (y == 1)) * predictions) - log_eps( n / (n2 - (y == 2)) * (1 - predictions))
#   #     gain[(start + 1) : end] <- shift_in_mean_and_variance(logits) * (end - start) / nrow(x)
#   #     gain
#   #   }
#   # }
# }


#' Closure generating function to calculate gains when splitting and learning a glasso model
#'
#' @param control an object of type \code{hdcd_control} returned by \link{hdcd_control}
#'
#' @return A closure with parameters \code{x}, \code{start} and \code{end}, that when evaluated
#' will itself return a closure with parameter \code{split_point}. This calculates the gain when
#' splitting the segment \code{(start, end]} of \code{x} at \code{split_point}. If the closure is
#' additionally supplied with \code{evaluate_all = TRUE}, an array of length \code{nrow(x)} is
#' returned with differences of loglikelihoods for each observation in \code{(start, end]} when split at
#' \code{split_point}.
# glasso_gain_function <- function(x, start, end, lambda, control = hdcd_control()) {
#
#   if(is.null(lambda)){
#     stop('Please supply a value for lambda to the glasso gain function')
#   }
#
#   control$n_obs <- nrow(x)
#
#   fit_global <- get_glasso_fit(x[(start + 1) : end, , drop = F], lambda = lambda, control = control)
#
#   function(split_point, split_candidates = split_point){
#
#     fit_left <- get_glasso_fit(x[(start + 1) : split_point, , drop = F], lambda = lambda, control = control)
#     fit_right <- get_glasso_fit(x[(split_point + 1) : end, , drop = F], lambda = lambda, control = control)
#
#     if(is.null(split_candidates)){
#
#       # only use variables for calculation of the gains curve that are used on left and right segment.
#       # See curve smoothing section in paper
#       (sum(loglikelihood(x[(start + 1) : split_point, fit_left$inds, drop = F],
#                          fit_global$mu[fit_left$inds, drop = F],
#                          fit_global$wi[fit_left$inds, fit_left$inds, drop = F])) +
#           sum(loglikelihood(x[(split_point + 1) : end, fit_right$inds, drop = F],
#                             fit_global$mu[fit_right$inds, drop = F],
#                             fit_global$wi[fit_right$inds, fit_right$inds, drop = F])) -
#           sum(loglikelihood(x[(start + 1) : split_point, fit_left$inds, drop = F],
#                             fit_left$mu[fit_left$inds, drop = F],
#                             fit_left$wi[fit_left$inds, fit_left$inds, drop = F])) -
#           sum(loglikelihood(x[(split_point + 1) : end, fit_right$inds, drop = F],
#                             fit_right$mu[fit_right$inds, drop = F],
#                             fit_right$wi[fit_right$inds, fit_right$inds, drop = F] ))) / nrow(x)
#
#     } else {
#       gain <- rep(NA, nrow(x))
#       cumsum_right <- cumsum(loglikelihood(x[(start + 1) : end, fit_right$inds, drop = F],
#                                            fit_right$mu[fit_right$inds, drop = F],
#                                            fit_right$wi[fit_right$inds, fit_right$inds, drop = F] ))
#       cumsum_left <- cumsum(loglikelihood(x[(start + 1) : end, fit_left$inds, drop = F],
#                         fit_left$mu[fit_left$inds, drop = F],
#                         fit_left$wi[fit_left$inds, fit_left$inds, drop = F]))
#       loglik_global <- sum(loglikelihood(x[(start + 1) : split_point, fit_left$inds, drop = F],
#                                          fit_global$mu[fit_left$inds, drop = F],
#                                          fit_global$wi[fit_left$inds, fit_left$inds, drop = F])) +
#         sum(loglikelihood(x[(split_point + 1) : end, fit_right$inds, drop = F],
#                           fit_global$mu[fit_right$inds, drop = F],
#                           fit_global$wi[fit_right$inds, fit_right$inds, drop = F]))
#
#       gain[split_candidates] <- (loglik_global - cumsum_right[length(cumsum_right)] - cumsum_left + cumsum_right)[split_candidates - start]/nrow(x)
#
#
#       gain
#
#     }
#   }
# }

#' Glasso Gain Function
#'
#' @inheritParams classifier_best_split_function
#' @inheritParams hdcd
#' @inheritParams classifier_gain_function
#'
#' @return A closure with parameters \code{x}, \code{start} and \code{end}, that when evaluated
#' will itself return a closure with arguments \code{split_point} and \code{split_candidates}.
#' Evaluated this function returns the gain when splitting the segment \code{(start, end]} of \code{x}
#'  at \code{split_point}.
#' @export
glasso_gain_function <-
  function(x, start, end, lambda, control = hdcd_control()) {
    if (is.null(lambda)) {
      stop('Please supply a value for lambda to the glasso gain function')
    }
    
    fit_global <-
      get_glasso_fit(x[(start + 1):end, , drop = F], lambda = lambda, control = control)

    function(split_point, split_candidates = split_point) {
      fit_left <-
        get_glasso_fit(x[(start + 1):split_point, , drop = F], lambda = lambda, control = control)
      fit_right <-
        get_glasso_fit(x[(split_point + 1):end, , drop = F], lambda = lambda, control = control)
      
      gain <- rep(NA, nrow(x))
      
      cumsum_right <-
        cumsum(loglikelihood(x[(start + 1):end, fit_right$inds, drop = F],
                             fit_right$mu[fit_right$inds, drop = F],
                             fit_right$wi[fit_right$inds, fit_right$inds, drop = F]))
      cumsum_left <-
        cumsum(loglikelihood(x[(start + 1):end, fit_left$inds, drop = F],
                             fit_left$mu[fit_left$inds, drop = F],
                             fit_left$wi[fit_left$inds, fit_left$inds, drop = F]))
      loglik_global <-
        sum(loglikelihood(x[(start + 1):split_point, fit_left$inds, drop = F],
                          fit_global$mu[fit_left$inds, drop = F],
                          fit_global$wi[fit_left$inds, fit_left$inds, drop = F])) +
        sum(loglikelihood(x[(split_point + 1):end, fit_right$inds, drop = F],
                          fit_global$mu[fit_right$inds, drop = F],
                          fit_global$wi[fit_right$inds, fit_right$inds, drop = F]))
      
      gain[(start + 1):end] <-
        (loglik_global - cumsum_right[length(cumsum_right)] - cumsum_left + cumsum_right) / nrow(x)
      best_split <-
        split_candidates[which.max(gain[split_candidates])]
      max_gain <- gain[best_split]
      
      list(gain = gain,
           best_split = best_split,
           max_gain = max_gain)
    }
  }

# kNN_best_split_function <- function(x, control = hdcd_control()){
#
#   distance_matrix <- as.matrix(dist(x)) # calculate pairwise distances of all observations
#
#   function(x, start, end, lambda, control){
#     n <-  end - start
#
#     k <- ceiling(sqrt(n))
#     # get matrix of neighbors. the ijs entry is monotonely increasing in the distance of the ith observation to
#     # the jth observation. Entries are unique in the columns
#     neighbors <- apply(distance_matrix[(start + 1) : end, (start + 1) : end], 1, function(x) order(order(x)))
#
#     # The k-th nearest neigbors are (in rows) those entries >=2 and <= k+1.
#     # The prediction for the jth observation is the proportion of how many of its k-nn are before the split,
#     # i.e the jth entry of the rowwise sum of the k-NN matrix up to split, i.e. sum(kNN_matrix[1 : split, j])
#
#     # thus predictions[j, ] are the predictions if the split_point is set at j.
#     predictions <- apply(neighbors >= 2 & neighbors <= k + 1, 2, cumsum) / k
#
#     function(split_point, split_candidates = split_point){
#
#
#       if(split_point > start + 1 & split_point < end){
#         sum(log_eps( (n - 1) / (split_point - start - 1) * predictions[split_point - start, 1 : (split_point - start)])) +
#           sum(log_eps( (n - 1) / (end - split_point) * (1 - predictions[split_point - start, (split_point - start + 1) : (end - start)])))
#       } else {
#         NA
#       }
#     }
#   }
#
# }

#' kNN gain function
#'
#'@importFrom stats dist
#'@inheritParams hdcd
kNN_gain_function <- function(x, control = hdcd_control()) {
  distance_matrix <-
    as.matrix(stats::dist(x)) # calculate pairwise distances of all observations
  
  function(x, start, end, lambda, control) {
    n <-  end - start
    
    k <- ceiling(sqrt(n))
    # get matrix of neighbors. the ijs entry is monotonely increasing in the distance of the ith observation to
    # the jth observation. Entries are unique in the columns
    neighbors <-
      apply(distance_matrix[(start + 1):end, (start + 1):end], 1, function(x)
        order(order(x)))
    
    # The k-th nearest neigbors are (in rows) those entries >=2 and <= k+1.
    # The prediction for the jth observation is the proportion of how many of its k-nn are before the split,
    # i.e the jth entry of the rowwise sum of the k-NN matrix up to split, i.e. sum(kNN_matrix[1 : split, j])
    
    # thus predictions[j, ] are the predictions if the split_point is set at j.
    predictions <-
      apply(neighbors >= 2 & neighbors <= k + 1, 2, cumsum) / k
    
    function(split_point, split_candidates = split_point) {
      if (split_point > start + 1 & split_point < end) {
        sum(log_eps((n - 1) / (split_point - start - 1) * predictions[split_point - start, 1:(split_point - start)])) +
          sum(log_eps((n - 1) / (end - split_point) * (1 - predictions[split_point - start, (split_point - start + 1):(end - start)])))
      } else {
        NA
      }
    }
  }
}
#' get_gain_function_from_loss_function
#'
#' Returns a closure with formal arguments \code{split_point}
#'
#' @inheritParams hdcd
#' @param lambda If \code{loss_function} has argument \code{lambda}, then this is the standard value used
#' for its evaluation as long as no different lambda is supplied
# gain_function_from_loss_function <- function(loss_function, lambda = NULL){
#
#   function(x, start, end, lambda = lambda){
#
#     stopifnot(!is.null(lambda))
#
#     # calculate global loss once to compare to later
#     global_loss <- loss_function(x[(start + 1) : end, , drop = F], lambda)
#
#     function(split_point){
#
#       # series of check to avoid unwanted behaviour
#       stopifnot(start <= split_point)
#       stopifnot(split_point <= end)
#
#       global_loss -
#         loss_function(x[(start + 1) : split_point, , drop = F], lambda = lambda) -
#         loss_function(x[(split_point + 1) : end, , drop = F], lambda = lambda)
#     }
#   }
# }
