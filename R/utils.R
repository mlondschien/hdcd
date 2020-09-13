#' Sample Folds
#'
#' @param n number of observations
#' @param k number of folds
#' @param randomize should folds be randomized? Equispaced folds are returned for randomize = FALSE
sample_folds <- function(n, k, randomize = FALSE) {
  if (k == 1) {
    as.factor(rep(0, n))
  } else  if (randomize) {
    random_draw <- runif(n)
    k_quantiles <- stats::quantile(random_draw, 0:k / k)
    cut(random_draw,
        k_quantiles,
        labels = 1:k,
        include.lowest = TRUE)
  } else {
    as.factor(rep(1:k, ceiling(n / k))[1:n])
  }
}

#' Bounded approximation to the logarithm
#'
#' @param x numeric
#' @param threshold bound
log_eps <- function(x, threshold = exp(-6)) {
  log((1 - threshold) * x + threshold)
}

#' Logarithmically Scaled Sequence Generation
#'
#' Generates a logarithmically scaled sequence
#'
#' @param from the starting value of the sequence
#' @param to the end value of the sequence
#' @param length.out the length of the sequence
log_space <- function(from, to, length.out) {
  exp(seq(
    from = log(from),
    to = log(to),
    length.out = length.out
  ))
}

#'Rand type performance indices
#'
#' Calculate rand type performance indices for two sets of changepoints. Typically one
#' of them will be the oracle estimate. See clues package for more details.
#'
#' @param cpts_a A sequence of changepoints.
#' @param cpts_b A sequence of changepoints.
#' @param n Total size of dataset from which both changepoint estimates originate.
#' @return Returns a vector of the index values.
#' @export
compare_change_points <- function(cpts_a, cpts_b, n) {
  if (!requireNamespace("clues", quietly = TRUE)) {
    stop("Please install clues: install.packages('clues')")
  } else {
    cpts_a <- setdiff(sort(cpts_a)[!duplicated(sort(cpts_a))], c(0, n))
    cpts_b <-
      setdiff(sort(cpts_b)[!duplicated(sort(cpts_b))], c(0, n))
    
    MarkGroupings <- function(cpts) {
      diffs <- c(cpts, n) - c(0, cpts)
      rep(1:length(diffs), diffs)
    }
    
    out <-
      c(
        clues::adjustedRand(MarkGroupings(cpts_a), MarkGroupings(cpts_b)),
        max_min_1 = max_min_distance(cpts_a, cpts_b, n) / n,
        max_min_2 = max_min_distance(cpts_b, cpts_a, n) / n
      )
    
    c(out, Hausdorff = max(out[6], out[7]))
  }
}

max_min_distance <- function(cpts_a, cpts_b, n) {
  if (length(cpts_b) == 0 & length(cpts_a) > 0) {
    n
  } else if (length(cpts_a) == 0) {
    0
  } else {
    max(sapply(cpts_a, function(i)
      min(sapply(cpts_b, function(j)
        abs(i - j)))))
  }
}


#' Draw segment (-boundaries) for WBS or BS
#'
#' @param n total number of observations available
#' @param n_segments Number of segments to be drawn (only relevant of WBS)
#' @param delta Minimal relative segment length
#' @param segmentation One of BS, WBS or SBS
#' @param alpha decay parameter in [1/2, 1) for SBS
draw_segments <-
  function(n,
           n_segments = n,
           delta = 0.1,
           segmentation = 'BS',
           alpha = 1 / 2) {
    if (segmentation == 'BS') {
      return(data.frame(start = 0, end = n))
      
    } else if (segmentation == 'SBS') {
      k <- ceiling(log(2 * delta) / log(alpha))
      df <- data.frame(start = 0, end = n)
      for (i in 1:k) {
        n_intervals <- ceiling(2 * (1 - alpha ^ i) / alpha ^ i) + 1
        length_interval <- round(n * alpha ^ i)
        df <- rbind(df,
                    data.frame(start = round(
                      seq(
                        from = 0,
                        to = n - length_interval,
                        length.out = n_intervals
                      )
                    ),
                    end = round(
                      seq(
                        from = length_interval,
                        to = n,
                        length.out = n_intervals
                      )
                    )))
      }
      return(df)
      
    } else if (segmentation == 'WBS') {
      df <- data.frame(start = 0, end = n)
      while (nrow(df) < n_segments) {
        ind <- sample(0:n, 2)
        if (abs(diff(ind)) >= 2 * delta * n) {
          df <- rbind(df, data.frame(start = min(ind), end = max(ind)))
        }
      }
    }
  }

#' find best split shift in mean and variance
#'
#' Efficiently calculates the gain in a one-dimensional shift in mean and variance scenario.
#'
#' @param x Array with entries that are assumed to have a shift in mean and variance at some split point.
#' @param alpha array of segment boundaries
#' @param train_fold array containing indices in training fold
#' @return An array on length \code{length(x)} with gains resulting from splitting at that specific split point.
#'
#' The negative gaussian loglikelihood of observations x for estimated mean and variance is \eqn{-n/2 * (log(2 \pi \hat\sigma^2) + 1)}.
#' The gaussian maximum likelihood estimate \eqn{\hat\sigma^2} is \code{(sum(x^2) - sum(x)^2/length(x))/length(x)}
shift_in_mean_and_variance <- function(x, ...){
  y <- x[!is.na(x)]
  n <- length(y)
  cumsum_x <- cumsum(y)
  cumsum_x_2 <- cumsum(y ^ 2)

  sigma_1 <- cumsum_x_2 / (1 : n) - cumsum_x ^ 2 / (1 : n) ^ 2
  sigma_2 <- (cumsum_x_2[n] - cumsum_x_2) / ((n - 1) : 0) - (cumsum_x[n] - cumsum_x) ^ 2 / ((n - 1) : 0)^2
  sigma_1 <- sigma_1 + 0.0001 * max(sigma_1[is.finite(sigma_1)])
  sigma_2 <- sigma_2 + 0.0001 * max(sigma_2[is.finite(sigma_2)])
  sigma_1[1] <- sigma_2[n] <- sigma_2[n-1] <- NA

  gain <- rep(NA, length(x))
  gain[!is.na(x)] <- (log(sigma_1[n]) - (1:n)/n * log(sigma_1) - ((n-1):0)/n * log(sigma_2))
  gain
}

#' Divide data into training + testing data
train_test_split <- function(x, alpha, train_fold = 1:nrow(x)) {
  if (is.logical(train_fold)) {
    train_fold <- which(train_fold)
  }
  if (alpha[1] != 0)
    alpha <- c(0, alpha)
  if (alpha[length(alpha)] != length(train_fold))
    alpha <- c(alpha, length(train_fold))
  
  segment_lengths <-  alpha[-1] - alpha[-length(alpha)]
  
  y_train <-
    rep.int(1:length(segment_lengths), times = segment_lengths)
  y <- rep(NA, nrow(x))
  y[train_fold] <- y_train
  for (i in 2:length(train_fold)) {
    if (y[train_fold[i - 1]] == y[train_fold[i]]) {
      y[train_fold[i - 1]:train_fold[i]] <- y[train_fold[i]]
    }
  }
  y[1:train_fold[1]] <- y[train_fold[1]]
  y[train_fold[length(train_fold)]:nrow(x)] <-
    y[train_fold[length(train_fold)]]
  
  y_test <- y[-train_fold]
  x_test <- x[-train_fold,][!is.na(y_test),]
  y_test <- y_test[!is.na(y_test)]
  
  
  stopifnot(nrow(x_test) == length(y_test))
  
  return(
    list(
      x_train = x[train_fold,],
      y_train = y_train,
      x_test = x_test,
      y_test = y_test,
      alpha = train_fold[alpha[-c(1, length(alpha))]]
    )
  )
}
