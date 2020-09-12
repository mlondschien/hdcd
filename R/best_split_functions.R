#' Get a get_best_split function from a classifier
#'
#' Returns a closure with formal arguments \code{x}, \code{start}, \code{end} and \code{split_candidates} that finds the
#' maximizer of the gain_function given data \code{x}, \code{start} and \code{end} on \code{split_candidates}
#' using the optimizer \code{optimizer}.
#'
#' @param classifier A function with formal arguments \code{x_train}, \code{y_train} and \code{control} that returns class probabilities in
#' a matrix of size nrow(x_train) x k, where k is the number of classes in y_train.
#' @inheritParams hdcd
#' @export
classifier_best_split_function <-
  function(classifier, optimizer, control) {
    stopifnot(optimizer %in% c('section_search', 'line_search', 'two_step_search'))
    
    function(x, start, end, split_candidates, lambda = 0) {
      # Sequence of checks to avoid unwanted behaviour
      stopifnot(start >= 0)
      stopifnot(end <= nrow(x))
      stopifnot(start < min(split_candidates))
      stopifnot(end > max(split_candidates))
      
      gain_function <-
        classifier_gain_function(
          x = x,
          start = start,
          end = end,
          classifier = classifier,
          control = control
        )
      n <- nrow(x)
      # Apply optimizer
      if (optimizer == 'section_search') {
        res <-
          section_search(
            gain_function = function(x)
              gain_function(x)$max_gain,
            split_candidates = split_candidates,
            n = n,
            control = control
          )
        #  } else if(optimizer == 'smooth_section_search'){
        #    res <- smooth_section_search(gain_function = function(x) gain_function(x)$gain[x],
        #                                 split_candidates = split_candidates, n = n, control = control)
      } else if (optimizer == 'line_search') {
        res <-
          line_search(
            gain_function = function(x)
              gain_function(x)$max_gain,
            split_candidates = split_candidates,
            n = n,
            control = control
          )
      } else if (optimizer == 'two_step_search') {
        res <- two_step_search(
          gain_function = gain_function,
          split_candidates = split_candidates,
          n = n,
          control = control
        )
      }
      
      # In some situations (e.g. with missing values) the gain function cannot be evaluated at any of the
      # split_candidates. In such a situation NA is returned as best_split and binary_segmentation is stopped.
      if (is.na(res$best_split)) {
        warning(
          'the gain_function could not be evaluated at any of the split_candidates. NA is returned and binary_segmentation
          stopped for segment (',
          start,
          ', ',
          end,
          ']',
          sep = ''
          )
      }
      
      res
    }
  }

#' Get a get_best_split function from the kNN classifier
#'
#' Returns a closure with formal arguments \code{x}, \code{start}, \code{end} and \code{split_candidates} that finds the
#' maximizer of the gain_function given data \code{x}, \code{start} and \code{end} on \code{split_candidates}
#' using the optimizer \code{optimizer}.
#'
#' @inheritParams hdcd
#' @export
kNN_best_split_function <- function(x, control) {
  distance_matrix <-
    as.matrix(dist(x)) # calculate pairwise distances of all observations
  
  function(x, start, end, split_candidates, lambda = 0) {
    n <-  end - start
    
    k <- min(ceiling(sqrt(n)), floor(n / 2))
    
    # get matrix of neighbors. the ijs entry is monotonely increasing in the distance of the ith observation to
    # the jth observation. Entries are unique in the columns
    # The ith observation is the j's (neighbors[i, j])-nearest neighbor
    neighbors <-
      apply(distance_matrix[(start + 1):end, (start + 1):end], 1, function(x)
        order(order(x)))
    
    neighbors <- neighbors >= 2 & neighbors <= k + 1
    # The k-th nearest neigbors are (in rows) those entries >=2 and <= k+1.
    # The prediction for the jth observation is the proportion of how many of its k-nn are before the split,
    # i.e the jth entry of the rowwise sum of the k-NN matrix up to split, i.e. sum(kNN_matrix[1 : split, j])
    
    # thus predictions[j, ] are the predictions if the split_point is set at j.
    predictions <- apply(neighbors, 2, cumsum) / k
    
    predictions <-
      (n - 1) * ((predictions / c(1, 1:(n - 1)) * !upper.tri(predictions)) + ((1 - predictions) / c((n - 1):1, 1) * upper.tri(predictions)))
    gain <- rep(NA, nrow(x))
    gain[(start + 1):end] <- rowSums(log_eps(predictions)) / nrow(x)
    best_split <-
      split_candidates[which.max(gain[split_candidates])]
    max_gain <- gain[best_split]
    
    n_perm <- control$permutation_test_n
    
    permutation_test <- array(NA, n_perm)
    set.seed(0)
    for (i in 1:n_perm) {
      sigma <- sample(1:nrow(x))
      sigma <- sigma[sigma > start & sigma <= end] - start
      neighbors_perm <- neighbors[sigma, sigma]
      predictions_perm <- apply(neighbors_perm, 2, cumsum) / k
      predictions_perm <-
        (n - 1) * ((
          predictions_perm / c(1, 1:(n - 1)) * !upper.tri(predictions_perm)
        ) + ((1 - predictions_perm) / c((n - 1):1, 1) * upper.tri(predictions_perm)
        ))
      permutation_test[i] <-
        max(rowSums(log_eps(predictions_perm))[split_candidates - start]) / nrow(x)
    }
    
    list(
      gain = gain,
      max_gain = max_gain,
      best_split = best_split,
      permutation_test = permutation_test,
      pval = (1 + sum(permutation_test >= max_gain)) / (1 + length(permutation_test))
    )
  }
}

#' Get a best_split_function from a gain_function
#'
#' @inheritParams hdcd
#' @param gain_function A function with argiument \code{x}, \code{start}, \code{end}, \code{lambda}, \code{control} that returns a closure with arguments
#'  split_point and possibly split_candidates that returns the gain after splitting the segment (\code{start}, \code{end}]
#'  at \code{split_point} given data \code{x} and tuning parameter \code{lambda}.
#' @return A function with formal arguments \code{x}, \code{start}, \code{end}, \code{split_candidates} and \code{lambda} that
#'  uses the optimizer specified to search for a maximum of the gain_function on the \code{split_candidates} given a segment
#'  (\code{start}, end], data \code{x} and a tuning parameter \code{lambda} and returns a list with arguments \code{gain}, \code{max_gain},
#'  \code{best_split} and possible \code{permutation_test} and \code{pval}.
best_split_function_from_gain_function <-
  function(gain_function, optimizer, control) {
    stopifnot(
      optimizer %in% c(
        'section_search',
        'smooth_section_search',
        'line_search',
        'two_step_search',
        'two_step_search_shift_in_mean',
        'two_step_search_adjusted'
      )
    )
    
    # Stop with error if gain_function does not take required arguments
    if (!all(c('x', 'start', 'end') %in% methods::formalArgs(gain_function))) {
      stop(
        'gain_function is not of the required form. Make sure that gain_function takes formal arguments x, start, end, lambda.'
      )
    }
    
    # Return closure that estimates /  calculates the location of the best split (with maximum gain) within (start, end]
    function(x, start, end, split_candidates, lambda = 0) {
      # Sequence of checks to avoid unwanted behaviour
      stopifnot(start >= 0)
      stopifnot(end <= nrow(x))
      stopifnot(start < min(split_candidates))
      stopifnot(end > max(split_candidates))
      
      n <- nrow(x)
      
      gain_function <-
        gain_function(
          x = x,
          start = start,
          end = end,
          lambda = lambda,
          control = control
        )
      
      # Apply optimizer
      if (optimizer == 'section_search') {
        res <-
          section_search(function(x)
            gain_function(x)$max_gain,
            split_candidates,
            n,
            control)
      } else if (optimizer == 'line_search') {
        res <-
          line_search(function(x)
            gain_function(x)$max_gain,
            split_candidates,
            n,
            control)
      } else if (optimizer == 'two_step_search') {
        res <- two_step_search(gain_function, split_candidates, n, control)
      } else {
        stop(
          'Make sure that optimizer is one of \'section_search\', \'line_search\' or \'two_step_search\''
        )
      }
      
      # In some situations (e.g. with missing values) the gain function cannot be evaluated at any of the
      # split_candidates. In such a situation NA is returned as best_split and binary_segmentation is stopped.
      if (is.na(res$best_split)) {
        warning(
          'the gain_function could not be evaluated at any of the split_candidates. NA is returned and binary_segmentation
          stopped for segment (',
          start,
          ', ',
          end,
          ']',
          sep = ''
          )
      }
      
      res
    }
  }