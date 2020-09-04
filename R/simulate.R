#' Simulate Observations from a model created by create_model
#'
#' @param model An object as created by \link{create_model}
#'
#' @return A n times p matrix of simulated data.
#' @export
simulate_from_model <- function(model) {
  segment_lengths <- model$segment_lengths
  
  data <-
    matrix(NA,
           nrow = sum(segment_lengths),
           ncol = length(model[["segment_means"]][[1]]))
  
  for (i in seq_along(segment_lengths)) {
    seg_start <- ifelse(i == 1, 1, sum(segment_lengths[(i - 1):1]) + 1)
    seg_end <- seg_start + segment_lengths[i] - 1
    data[seg_start:seg_end,] <-
      MASS::mvrnorm(segment_lengths[[i]], model[["segment_means"]][[i]], model[["cov_mats"]][[i]])
  }
  return(data)
}

#' Create GGM with changepoints
#'
#' Create a model to generate data from a GGM for simulating the detection of changepoints
#'
#' @param n Number of observations
#' @param p Number of dimensions
#' @param changepoints An array with the changepoints for the model to be generated. The change points are the last
#' observation of each segment.
#' @param mean_vecs If NULL, the mean for each segment will be zero.
#' Otherwise mean_vecs should be a list containing a p-dimensional numeric vector with the means for each segment.
#' @param model_function A function that spawns covariance matrices of dimension p.
#' @param ... Addtional arguments to be supplied to modelFUN
#'
#' @return An object to be used by \link{simulate_from_model}
#' @export
create_model <-
  function(n,
           p,
           changepoints,
           model_function,
           mean_vecs = NULL,
           ...) {
    model_args <- list(...)
    
    # make sure changepoints are adapt to setting
    if (!(length(changepoints) == 0 || (changepoints < n))) {
      stop('Changepoints cannot be greater than n')
    }
    
    if (!is.null(mean_vecs) &&
        length(mean_vecs) != length(changepoints) + 1) {
      stop('Please make sure mean_vecs is a list with one element for each segment')
    }
    
    if (is.null(mean_vecs)) {
      mean_vecs <-
        replicate(length(changepoints) + 1, rep(0, p), simplify = F)
    }
    
    segment_lengths <-
      c(changepoints - c(0, changepoints[-length(changepoints)]), n - changepoints[length(changepoints)])
    
    if (length(segment_lengths) == 0) {
      segment_lengths <- n
    }
    cov_mats <-
      replicate(length(changepoints) + 1,
                do.call(model_function, c(list(p = p), model_args)),
                simplify = F)
    
    list(
      segment_lengths = segment_lengths,
      segment_means = mean_vecs,
      cov_mats = cov_mats,
      true_changepoints = changepoints
    )
  }

#' Delete values from a design matrix
#' 
#' \code{delete_values} returns the matrix \code{x} with a proportion of \code{m} values deleted.
#'
#' @param x Design matrix
#' @param m proportion of values to be deleted
#' @param missingness One of \code{'mcar'}, \code{'nmar'} or \code{'blockwise'}. For \code{mcar} the values are
#' deleted uniformly at random. If \code{nmar} is supplied, the probability of deletion of a given value is set 
#' to be proportional to its absolute value. If \code{blockwise} is supplied, repeatedly values \code{l ~ Pois(p / 20)}
#' and \code{k ~ Exp(n/8)} are drawn and a block of size \code{l x k} is deleted.
#' @importFrom stats rpois
#' @importFrom stats rexp
#' @export
delete_values <- function(x, m, missingness = 'mcar'){
  n <- nrow(x)
  p <- ncol(x)
  
  if(!(missingness %in% c('mcar', 'nmar', 'blockwise', 'none', 'blockwise_big'))){
    stop('missingness needs to be one of \'mcar\', \'nmar\', \'blockwise\' or \'none\'.')
  }
  
  if (missingness == 'mcar'){
    # randomly draw indices corresponding to available observations to later delete
    inds <- sample(which(!is.na(x)), floor(m * n * p), replace = F)
  
  } else if (missingness == 'nmar'){
    # delete values with probability proportional to absolute value of entry
    i <- which(!is.na(x))
    inds <- sample(i, floor(m * n * p), replace = F, prob = abs(c(x[i])))
  } else if (missingness == 'blockwise'){
    total <- prod(dim(x))
    missing <- sum(is.na(x))
    m <- m + missing / total
    missing_max <- ceiling(m * total)
    while (missing < missing_max){
      o <- stats::rpois(1, p / 20)
      l <- min(floor(stats::rexp(1, 8/n)), ceiling((missing_max - missing)/o))
      j <- sample(p, o)###lala
      k <- sample((1 - floor(l/2)) : (n + floor(l/2) - 1), 1)
      int <- max((k - floor(l/2)), 1) : min((k + floor(l/2)), n)
      x[int, j] <- NA
      missing <- sum(is.na(x))
    }
    return(x)
} else if (missingness == 'blockwise_big'){
    total <- prod(dim(x))
    missing <- sum(is.na(x))
    m <- m + missing / total
    missing_max <- ceiling(m * total)
    while (missing < missing_max){
      o <- rpois(1, p / 5)
      l <- min(floor(rexp(1, 8/n)), ceiling((missing_max - missing)/o))
      j <- sample(p, o)###lala
      k <- sample((1 - floor(l/2)) : (n + floor(l/2) - 1), 1)
      int <- max((k - floor(l/2)), 1) : min((k + floor(l/2)), n)
      x[int, j] <- NA
      missing <- sum(is.na(x))
    }
    return(x)
  } else if (missingness == 'none'){
    inds <- numeric(0)
  }
  
  x_del <- x
  x_del[inds] <- NA
  
  return(x_del)
}

#' Plot the missingness structure of a design matrix
#' 
#' @param x design matrix
#' @param order Should the variables of x be reoderd to better display blocks?
#' @import data.table
#' @importFrom stats dist
#' @export
plot_missingness_structure <- function(x, order = F){
  index <- variable <- value <- NULL
  if (! requireNamespace("TSP", quietly = TRUE)) {
    stop("Please install TSP: install.packages('TSP')")
  } else  if (! requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install ggplot2: install.packages('ggplot2')")
  } else {
  if (order){
    y <- is.na(x)
    y <- cbind(y, y[, 1])
    y <- abs(y[, -1] - y[, -ncol(y)])
    distance <- 2 * stats::dist(t(y)) + stats::dist(t(is.na(x)))
    tsp <- TSP::TSP(distance)
    tour <- TSP::solve_TSP(tsp)
    x <- x[, c(1, TSP::cut_tour(tour, 1))]
  }
  dt <- data.table::data.table(is.na(x))
  colnames(dt) <-  as.character(1 : ncol(x))
  dt$index <- 1 : nrow(x)
  dt <- data.table::melt(dt, id.vars = 'index')
  dt$value <-  ifelse(dt$value, 'missing', 'not_missing')
  ggplot2::ggplot(dt, ggplot2::aes(x = index, y = variable, col = value)) + ggplot2::geom_point()
  }
}

#' simulate (non-parametric) data sets
#' 
#' @param model One if \code{Iris}, \code{Noise}, \code{Wine_1}, \code{Wine_2}, \code{RN},
#' \code{CIM} or \code{CIC}
#' @param seed seed used for dataset creation
#' @export
simulate_non_parametric <- function(model, seed){
  
  set.seed(seed)
  
  if(model == 'Iris'){
    i <- c(sapply(0 : 2, function(i) i * 50 + sample(1 : 50)))
    x <- as.matrix(datasets::iris[i, -5])
    x[1 : 50, ] <- scale(x[1 : 50, ], scale = F)
    x[51 : 100, ] <- scale(x[51 : 100, ], scale = F)
    x[101 : 150, ] <- scale(x[101 : 150, ], scale = F)
    x <- cbind(x, matrix(c(stats::rnorm(150), runif(150)), nrow = 150))
    list(x = scale(x), true_change_points = c(50, 100))
  # } else  if(model == 'iris_C_6'){
  #   i <- c(sapply(0 : 2, function(i) i * 50 + sample(1 : 50)))
  #   x <- as.matrix(iris[i, -5])
  #   x[1 : 50, ] <- scale(x[1 : 50, ], scale = F)
  #   x[51 : 100, ] <- scale(x[51 : 100, ], scale = F)
  #   x[101 : 150, ] <- scale(x[101 : 150, ], scale = F)
  #   x <- cbind(x, matrix(c(rnorm(3 * 150), runif(3 * 150)), nrow = 150))
  #   list(x = scale(x), true_change_points = c(50, 100))
  # } else  if(model == 'iris_C_18'){
  #   i <- c(sapply(0 : 2, function(i) i * 50 + sample(1 : 50)))
  #   x <- as.matrix(iris[i, -5])
  #   x[1 : 50, ] <- scale(x[1 : 50, ], scale = F)
  #   x[51 : 100, ] <- scale(x[51 : 100, ], scale = F)
  #   x[101 : 150, ] <- scale(x[101 : 150, ], scale = F)
  #   x <- cbind(x, matrix(c(rnorm(9 * 150), runif(9 * 150)), nrow = 150))
  #   list(x = x, true_change_points = c(50, 100))
  # } else if(model == 'iris_CS_2'){
  #   i <- c(sapply(0 : 2, function(i) i * 50 + sample(1 : 50)))
  #   x <- as.matrix(iris[i, -5])
  #   x[1 : 50, ] <- scale(x[1 : 50, ], scale = T)
  #   x[51 : 100, ] <- scale(x[51 : 100, ], scale = T)
  #   x[101 : 150, ] <- scale(x[101 : 150, ], scale = T)
  #   x <- cbind(x, matrix(c(rnorm(150), runif(150)), nrow = 150))
  #   list(x = scale(x), true_change_points = c(50, 100))
  # } else if(model == 'iris_CS_6'){
  #   i <- c(sapply(0 : 2, function(i) i * 50 + sample(1 : 50)))
  #   x <- as.matrix(iris[i, -5])
  #   x[1 : 50, ] <- scale(x[1 : 50, ], scale = T)
  #   x[51 : 100, ] <- scale(x[51 : 100, ], scale = T)
  #   x[101 : 150, ] <- scale(x[101 : 150, ], scale = T)
  #   x <- cbind(x, matrix(c(3 * rnorm(150), 3 * runif(150)), nrow = 150))
  #   list(x = x, true_change_points = c(50, 100))
  # } else if(model == 'iris_CS_18'){
  #   i <- c(sapply(0 : 2, function(i) i * 50 + sample(1 : 50)))
  #   x <- as.matrix(iris[i, -5])
  #   x[1 : 50, ] <- scale(x[1 : 50, ], scale = T)
  #   x[51 : 100, ] <- scale(x[51 : 100, ], scale = T)
  #   x[101 : 150, ] <- scale(x[101 : 150, ], scale = T)
  #   x <- cbind(x, matrix(c(9 * rnorm(150), 9 * runif(150)), nrow = 150))
  #   list(x = scale(x), true_change_points = c(50, 100))
  #   
  #   
  # } else if (model == 'shift_in_covariance_2'){
  #   Sigma <- diag(2)
  #   Sigma[1,2] <- Sigma[2,1] <- 0.9
  #   x <- rbind(
  #     MASS::mvrnorm(n = 100, mu = rep(0, 2), Sigma = diag(2)),
  #     MASS::mvrnorm(n = 100, mu = rep(0, 2), Sigma = Sigma),
  #     MASS::mvrnorm(n = 100, mu = rep(0, 2), Sigma = diag(2))
  #   )
  #   list(x = x, true_change_points = c(100, 200))
  } else if (model == 'Noise'){
    Sigma <- diag(5)
    Sigma[1,2] <- Sigma[2,1] <- 0.9
    x <- rbind(
      MASS::mvrnorm(n = 100, mu = rep(0, 5), Sigma = diag(5)),
      MASS::mvrnorm(n = 100, mu = rep(0, 5), Sigma = Sigma),
      MASS::mvrnorm(n = 100, mu = rep(0, 5), Sigma = diag(5))
    )
    list(x = x, true_change_points = c(100, 200))
  # } else if (model == 'shift_in_covariance_9'){
  #   Sigma <- diag(9)
  #   Sigma[1,2] <- Sigma[2,1] <- 0.9
  #   x <- rbind(
  #     MASS::mvrnorm(n = 100, mu = rep(0, 9), Sigma = diag(9)),
  #     MASS::mvrnorm(n = 100, mu = rep(0, 9), Sigma = Sigma),
  #     MASS::mvrnorm(n = 100, mu = rep(0, 9), Sigma = diag(9))
  #   )
  #   list(x = x, true_change_points = c(100, 200))
  #   
  #   
  } else if (model == 'Wine_1'){
    n <- 500
    segment_lengths <- sample(c(120, 120, 70, 190))
    true_change_points <- c(0,cumsum(segment_lengths))
    qualities <- sample(4 : 8, 4, replace = F)
    
    wine_dataset$type <- ifelse(wine_dataset$type == 'red', 0, 1)
    quality <- wine_dataset$quality
    id <- 1 : nrow(wine_dataset)
    
    ids <- c(sample(id[quality == qualities[1]], segment_lengths[1]),
             sample(id[quality == qualities[2]], segment_lengths[2]),
             sample(id[quality == qualities[3]], segment_lengths[3]),
             sample(id[quality == qualities[4]], segment_lengths[4])
    )
    stopifnot(all(unique(ids)))
    x <- as.matrix(wine_dataset[ids, -13])
    list(x = scale(x), true_change_points = true_change_points[c(-1, -length(true_change_points))])
  # } else if (model == 'wine_500_12'){
  #   n <- 500
  #   segment_lengths <- sample(c(120, 120, 70, 190))
  #   true_change_points <- c(0,cumsum(segment_lengths))
  #   qualities <- sample(4 : 8, 4, replace = F)
  #   
  #   load('./data/wine_dataset.rda')
  #   
  #   wine_dataset$type <- ifelse(wine_dataset$type == 'red', 0, 1)
  #   quality <- wine_dataset$quality
  #   id <- 1 : nrow(wine_dataset)
  #   
  #   ids <- c(sample(id[quality == qualities[1]], segment_lengths[1]),
  #            sample(id[quality == qualities[2]], segment_lengths[2]),
  #            sample(id[quality == qualities[3]], segment_lengths[3]),
  #            sample(id[quality == qualities[4]], segment_lengths[4])
  #   )
  #   stopifnot(all(unique(ids)))
  #   x <- as.matrix(wine_dataset[ids, -13])
  #   x <- cbind(x, matrix(c(runif(500 * 6), rnorm(500 * 6)), nrow = 500))
  #   list(x = scale(x), true_change_points = true_change_points[c(-1, -length(true_change_points))])
  # } else if (model == 'wine_500_C_0'){
  #   n <- 500
  #   segment_lengths <- sample(c(120, 120, 70, 190))
  #   true_change_points <- c(0,cumsum(segment_lengths))
  #   qualities <- sample(4 : 8, 4, replace = F)
  #   
  #   load('./data/wine_dataset.rda')
  #   
  #   wine_dataset$type <- ifelse(wine_dataset$type == 'red', 0, 1)
  #   quality <- wine_dataset$quality
  #   id <- 1 : nrow(wine_dataset)
  #   
  #   x <- rbind(
  #     as.matrix(center(wine_dataset[sample(id[quality == qualities[1]], segment_lengths[1])]), scale = F),
  #     as.matrix(center(wine_dataset[sample(id[quality == qualities[2]], segment_lengths[2])]), scale = F),
  #     as.matrix(center(wine_dataset[sample(id[quality == qualities[3]], segment_lengths[3])]), scale = F),
  #     as.matrix(center(wine_dataset[sample(id[quality == qualities[4]], segment_lengths[4])]), scale = F)
  #   )
  #   
  #   list(x = x, true_change_points = true_change_points[c(-1, -length(true_change_points))])
  # } else if (model == 'wine_500_C_12'){
  #   n <- 500
  #   segment_lengths <- sample(c(120, 120, 70, 190))
  #   true_change_points <- c(0,cumsum(segment_lengths))
  #   qualities <- sample(4 : 8, 4, replace = F)
  #   
  #   load('./data/wine_dataset.rda')
  #   
  #   wine_dataset$type <- ifelse(wine_dataset$type == 'red', 0, 1)
  #   quality <- wine_dataset$quality
  #   id <- 1 : nrow(wine_dataset)
  #   
  #   x <- rbind(
  #     as.matrix(center(wine_dataset[sample(id[quality == qualities[1]], segment_lengths[1])]), scale = F),
  #     as.matrix(center(wine_dataset[sample(id[quality == qualities[2]], segment_lengths[2])]), scale = F),
  #     as.matrix(center(wine_dataset[sample(id[quality == qualities[3]], segment_lengths[3])]), scale = F),
  #     as.matrix(center(wine_dataset[sample(id[quality == qualities[4]], segment_lengths[4])]), scale = F)
  #   )[, -13]
  #   
  #   x <- cbind(x, matrix(c(runif(500 * 6), rnorm(500 * 6)), nrow = 500))
  #   
  #   list(x = x, true_change_points = true_change_points[c(-1, -length(true_change_points))])
  #   
  } else if (model == 'Wine_2'){
    n <- 1000
    segment_lengths <- sample(c(120, 120, 120, 310, 70, 70, 190))
    true_change_points <- c(0,cumsum(segment_lengths))
    
    quality_0 <- sample(0 : 2, 1)
    qualities <- sample(1 : 2, 6, replace = T)
    qualities <- ( quality_0 + c(0, cumsum(qualities)) ) %% 3 + 5
    
    wine_dataset$type <- ifelse(wine_dataset$type == 'red', 0, 1)
    quality <- wine_dataset$quality
    id <- 1 : nrow(wine_dataset)
    
    ids <- c(sample(id[quality == qualities[1]], segment_lengths[1]),
             sample(id[quality == qualities[2]], segment_lengths[2]),
             sample(id[quality == qualities[3]], segment_lengths[3]),
             sample(id[quality == qualities[4]], segment_lengths[4]),
             sample(id[quality == qualities[5]], segment_lengths[5]),
             sample(id[quality == qualities[6]], segment_lengths[6]),
             sample(id[quality == qualities[7]], segment_lengths[7])
    )
    
    x <- as.matrix(wine_dataset[ids, -13])
    
    list(x = scale(x), true_change_points = true_change_points[c(-1, -length(true_change_points))])
  # } else if (model == 'wine_1000_12'){
  #   n <- 1000
  #   segment_lengths <- sample(c(120, 120, 120, 310, 70, 70, 190))
  #   true_change_points <- c(0,cumsum(segment_lengths))
  #   
  #   quality_0 <- sample(0 : 2, 1)
  #   qualities <- sample(1 : 2, 6, replace = T)
  #   qualities <- ( quality_0 + c(0, cumsum(qualities)) ) %% 3 + 5
  #   
  #   load('./data/wine_dataset.rda')
  #   
  #   wine_dataset$type <- ifelse(wine_dataset$type == 'red', 0, 1)
  #   quality <- wine_dataset$quality
  #   id <- 1 : nrow(wine_dataset)
  #   
  #   ids <- c(sample(id[quality == qualities[1]], segment_lengths[1]),
  #            sample(id[quality == qualities[2]], segment_lengths[2]),
  #            sample(id[quality == qualities[3]], segment_lengths[3]),
  #            sample(id[quality == qualities[4]], segment_lengths[4]),
  #            sample(id[quality == qualities[5]], segment_lengths[5]),
  #            sample(id[quality == qualities[6]], segment_lengths[6]),
  #            sample(id[quality == qualities[7]], segment_lengths[7])
  #   )
  #   
  #   x <- as.matrix(wine_dataset[ids, -13])
  #   x <- cbind(x, matrix(c(runif(1000 * 6), rnorm(1000 * 6)), nrow = 1000))
  #   list(x = x, true_change_points = true_change_points[c(-1, -length(true_change_points))])
  # 
  } else if (model == 'chain_network'){
    segment_lengths <- sample(c(120, 120, 70, 190))
    true_change_points <- c(0,cumsum(segment_lengths))
    alpha <- true_change_points[-c(1,5)]
    x <- simulate_from_model(create_model(500, 100, alpha, ChainNetwork))
    list(x = x, true_change_points = alpha)
  } else if (model == 'RN'){
    segment_lengths <- sample(c(120, 120, 70, 190))
    true_change_points <- c(0,cumsum(segment_lengths))
    alpha <- true_change_points[-c(1,5)]
    x <- simulate_from_model(create_model(500, 100, alpha, RandomNetwork))
    list(x = x, true_change_points = alpha)
  } else if (model == 'CIM'){
    
    x <- rbind(
      MASS::mvrnorm(n = 100, mu = rep(0, 2), Sigma = diag(2)),
      MASS::mvrnorm(n = 100, mu = rep(1, 2), Sigma = diag(2)),
      MASS::mvrnorm(n = 100, mu = rep(0, 2), Sigma = diag(2))
    )
    list(x = x, true_change_points = c(100, 200))
    
  } else if (model == 'CIC'){
    Sigma <- diag(2)
    Sigma[1,2] <- Sigma[2,1] <- 0.7
    
    x <- rbind(
      MASS::mvrnorm(n = 100, mu = rep(0, 2), Sigma = diag(2)),
      MASS::mvrnorm(n = 100, mu = rep(0, 2), Sigma = Sigma),
      MASS::mvrnorm(n = 100, mu = rep(0, 2), Sigma = diag(2))
    )
    list(x = x, true_change_points = c(100, 200))
  } else if (model == 'Gauss_1cp'){
    d <- 20
    Sigma <- diag(d)
    Sigma[1, 1] <- d^(1 / 3)
    delta <- 0.2
    x <- rbind(
      MASS::mvrnorm(80, rep(0, d), diag(d)),
      MASS::mvrnorm(120, rep(delta, d), Sigma)
    )
    list(x = x, true_change_points = 80)
  } else if (model == 'Wine_1cp'){
    load('./data/wine_dataset.rda')

    type <- wine_dataset$type
    id <- 1 : nrow(wine_dataset)
    x <- rbind(
      as.matrix(scale(wine_dataset[sample(id[type == 'red'], 80), -1], scale = T)),
      as.matrix(scale(wine_dataset[sample(id[type == 'white'], 120), -1], scale = T))
    )
    list(x = scale(x), true_change_points = 80)
  } else if (model == 'Gauss_null'){
    x <- MASS::mvrnorm(200, rep(0, 20), diag(20))
    list(x = x, true_change_points = NULL)
  } else if (model == 'Wine_null'){
    x <- load('./data/wine_dataset.rda')
    quality <- sample(5 : 7, 1)
    id <- (1 : nrow(wine_dataset))[wine_dataset$quality == quality]
    wine_dataset$type <- ifelse(wine_dataset$type == 'red', 0, 1)
    list(x = scale(as.matrix(wine_dataset[sample(id, 200, replace = F), -13])), true_change_points = NULL)
  }
}


#' ChainNetwork
#'
#' Spawn a chain network covariance matrix
#'
#' @param p Positive integer.The desired number of dimensions.
#' @param n_perm Positive integer. The first n_perm dimensions will be permuted randomly.
#' @param a Positive float between 0 and 1. Scale parameter for the elements of the covariance matrix
#' @param prec_mat Should the precision matrix be returned? If false the covariance matrix will be returned (default).
#' @param scaled If TRUE the created precision matrix will be inverted and scaled to a correlation matrix.
#'
#' @return A covariance or precision matrix.
#' @export
#'
#' @importFrom stats runif cov2cor
#'
#' @examples
#' ChainNetwork(50)
ChainNetwork <- function(p, n_perm = p, a = 0.5, prec_mat = F, scaled = T) {
  stopifnot(p >= n_perm)
  s_vec <- cumsum(runif(p, 0.5, 1))
  if (!is.null(n_perm) && n_perm >= 0) {
    perm_inds <- sample(1:n_perm, n_perm, replace = FALSE)
    
    if (n_perm < p) perm_inds <- c(perm_inds, (n_perm + 1):p)
    
    s_vec <- s_vec[perm_inds]
  }
  
  omega <- matrix(0, nrow = p, ncol = p)
  
  for (i in seq_len(p)) {
    for (j in seq_len(p)) {
      omega[i, j] <- exp(-a * abs(s_vec[i] - s_vec[j]))
    }
  }
  
  if (scaled) {
    sigma <- cov2cor(solve(omega))
    if (prec_mat) {
      solve(sigma)
    } else {
      sigma
    }
  } else {
    if (prec_mat) {
      omega
    } else {
      solve(omega)
    }
  }
}

#' ScaleNetwork
#'
#' Spawn a hub network covariance matrix
#'
#' @inheritParams ChainNetwork
#' @inheritParams RandomNetwork
#' @param preferential_power Power coefficient alpha for weighting of degree number as alpha in prefential attachment mechanism.
#'
#' @return A covariance or precision matrix.
#' @export
#'
#' @examples
#' ScaleNetwork(50)
ScaleNetwork <- function(p, preferential_power = 1, u = 0.1, v = 0.3, prec_mat = T, scaled = F) {
  theta <- matrix(0, p, p)
  probs <- numeric(p)
  
  theta[1, 2] <- theta[2, 1] <- TRUE
  
  for (i in seq(3, p)) {
    probs <- colSums(theta)^preferential_power
    probs <- probs / sum(probs)
    
    edges <- sample.int(i - 1, 1, prob = probs[1:(i - 1)])
    
    theta[edges, i] <- theta[i, edges] <- TRUE
  }
  
  diag(theta) <- 0
  
  omega <- theta * v
  
  diag(omega) <- abs(min(eigen(omega)$values)) + u
  
  if (scaled) {
    sigma <- cov2cor(solve(omega))
    if (prec_mat) {
      solve(sigma)
    } else {
      sigma
    }
  } else {
    if (prec_mat) {
      omega
    } else {
      solve(omega)
    }
  }
}

#' RandomNetwork
#'
#' Spawn a erdos renyi type random network covariance matrix
#'
#' @inheritParams ChainNetwork
#' @param prob Probabilty that a pair of nodes have a common edge.
#' @param u Constant added to the diagonal elements of the precision matrix for controlling the magnitude of partial correlations.
#' @param v Constant added to the off diagonal of the precision matrix for controlling the magnitude of partial correlations.
#'
#' @return A covariance matrix.
#' @export
#'
#' @examples
#' RandomNetwork(50)
RandomNetwork <- function(p, prob = min(1, 5 / p), u = 0.1, v = 0.3, prec_mat = F, scaled = F) {
  theta <- matrix(0, p, p)
  
  tmp <- matrix(runif(p^2, 0, 0.5), p, p)
  tmp <- tmp + t(tmp)
  theta[tmp < prob] <- 1
  diag(theta) <- 0
  
  omega <- theta * v
  diag(omega) <- abs(min(eigen(omega)$values)) + u
  
  if (scaled) {
    sigma <- cov2cor(solve(omega))
    if (prec_mat) {
      solve(sigma)
    } else {
      sigma
    }
  } else {
    if (prec_mat) {
      omega
    } else {
      solve(omega)
    }
  }
}

#' DiagMatrix
#'
#' Spawn a p-dimensional identity matrix.
#'
#' @inheritParams ChainNetwork
#'
#' @return A covariance matrix.
#' @export
#'
#' @examples
#' DiagMatrix(50)
DiagMatrix <- function(p) {
  diag(p)
}


#' MoveEdges
#'
#' Radomly move a share of the edges in a random graph
#'
#' In order to create a slightly different graph with the same level of sparsity
#' the selected share of edges will be randomly moved to positions where no edge existed
#' before. Make sure to choose share_moves low for dense graphs.
#'
#' @param x A precision matrix
#' @param share_moves Share of edges to be moved.
#' @param tol Tolerance for zero values.
#'
#' @export
MoveEdges <- function(x, share_moves = 0.1, tol = 1e-16) {
  if (share_moves == 0) return(x)
  
  edges <- which(upper.tri(x) & abs(x) >= tol)
  not_edges <- which(upper.tri(x) & x == 0)
  
  n_moves <- floor(share_moves * length(edges))
  
  if (length(not_edges) <= n_moves) {
    n_moves <- length(not_edges)
    warning("Cannot move edges because graph is not sparse enough. Will move as many as possible.")
  }
  
  d <- diag(x)
  
  # ensure that matrix is not pd to begin with
  diag(x) <- diag(x) - abs(min(eigen(x)$values)) - 0.1
  
  # sample until matrix is pd again
  while (any(eigen(x)$values <= 0)) {
    sel_edges <- sample(edges, n_moves)
    x[sample(not_edges, n_moves)] <- x[sel_edges]
    x[sel_edges] <- 0
    
    # make symmetric again
    x[lower.tri(x)] <- t(x)[lower.tri(x)]
    diag(x) <- d
  }
  x
}


#' RegrowNetwork
#'
#' Prune graph and regrow edges according of scale-free network
#'
#' Note that v and preferential_power agruments need to be equal to the ones which
#' initially created omega.
#'
#' @param omega A precision matrix as created by ScaleNetwork
#' @param n_nodes Number of nodes to prune and regrow. Default is 0.1 of all nodes.
#' @inheritParams ScaleNetwork
#' @inheritParams RandomNetwork
#'
#' @export
#'
#' @examples
#' omega <- ScaleNetwork(20, v = 1)
#' omega_hat <- RegrowNetwork(omega, v = 1)
#' omega
#' omega_hat
RegrowNetwork <- function(omega, n_nodes = ncol(omega) * 0.1, preferential_power = 1, v = 0.3) {
  d <- diag(omega)
  diag(omega) <- 0
  
  n_nodes <- floor(n_nodes)
  
  p <- ncol(omega)
  
  # prune
  pruned <- numeric(n_nodes)
  for (i in seq(1, n_nodes)) {
    candidates <- which(colSums(omega != 0) == 1)
    if (length(candidates) == 0) stop("Not enough nodes to prune from graph!")
    pruned[i] <- sample(candidates, 1)
    omega[, pruned[i]] <- omega[pruned[i], ] <- 0
  }
  
  # regrow
  pruned <- rev(pruned)
  probs <- numeric(p)
  
  omega_temp <- omega
  
  # Regrow network and discard permutation if it leads to non pd matrix
  while (any(eigen(omega)$values <= 0)) {
    omega <- omega_temp
    for (i in seq(1, n_nodes)) {
      probs <- colSums(omega != 0)^preferential_power
      if (sum(probs) == 0 | any(is.na(probs))) probs <- rep(1, p)
      probs <- probs / sum(probs)
      
      edges <- sample.int(p, 1, prob = probs)
      
      omega[edges, pruned[i]] <- omega[pruned[i], edges] <- v
    }
    diag(omega) <- d
  }
  omega
}

