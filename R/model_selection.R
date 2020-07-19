#' Loglikelihood from classifier
#' 
#' @param y_train training labels
#' @param predictions matrix of predicted class probabilities
#' @param y_test test labels
#' @param oob boolean whether predictions are oob
classifier_loglikelihood <- function(y_train, predictions, y_test = y_train, oob = F){
  
  n_train <- length(y_train)
  n_test <- length(y_test)
  
  k <- ncol(predictions)
  
  stopifnot(k == length(unique(y_train)))
  stopifnot(nrow(predictions) == n_test)
  
  indices <- matrix(rep(1 : k, each = n_test), nrow = n_test) == matrix(rep(y_test, times = k), nrow = n_test)
  
  if(oob){
    stopifnot(all(y_train == y_test))
    pi <- (n_train - 1) / (matrix(rep(table(y_train), each = n_train), nrow = n_train) - indices)
  } else {
    pi <- n_train / matrix(rep(table(y_train), each = n_test), nrow = n_test)
  }
  
  likelihood_matrix <- log_eps(predictions * pi)
  
  gain <- sum(likelihood_matrix[indices])
  
  list(likelihood_matrix = likelihood_matrix, gain = gain)
}
# 
# classifier_model_select <- function(tree, x, classifier, control, fold = 1 : nrow(x)){
#   
#   alpha <- get_splits(tree)
#   
#   gain_frame <- data.frame(id = 1 : length(alpha))
#   gain_frame$k <- sapply(alpha, length)
#   
#   f <- function(alpha){
#     
#    # if(is.null(fold)){
#    #   
#    #    if(length(alpha) == 0){
#    #      return(list(train_gain = 0))
#    #    }
#    #    alpha <- c(0, alpha, nrow(x))
#    #    segment_lengths <- alpha[-1] - alpha[-length(alpha)]
#    #    x_train <- x
#    #    y_train <- rep.int(1 : length(segment_lengths), times = segment_lengths)
#    #    
#    #    predictions <- classifier(x_train = x_train, y_train = y_train, control = control)
#    #    
#    #    list(train_gain = classifier_loglikelihood(y_train = y_train, predictions = predictions$train, oob = control$classifier_oob)$gain)
#    #         
#    #  } else {
#    #    
#       if(length(alpha) == 0){
#         return(list(train_gain = 0, test_gain = 0))
#       }
#       
#       train_test_data <- train_test_split(x, alpha, fold)
#       x_train <- train_test_data$x_train
#       x_test <- train_test_data$x_test
#       y_train <- train_test_data$y_train
#       y_test <- train_test_data$y_test
#       
#       if(nrow(x_test) > 0){
#         
#         predictions <- classifier(x_train = x_train, y_train = y_train, x_test = x_test, control = control)
#           
#         list(train_gain = classifier_loglikelihood(y_train = y_train, predictions = predictions$train, oob = control$classifier_oob)$gain,
#              test_gain = classifier_loglikelihood(y_train = y_train, predictions = predictions$test, y_test = y_test)$gain * (nrow(x) - nrow(x_train)) / nrow(x_test))
#       } else {
#         predictions <- classifier(x_train = x_train, y_train = y_train, control = control)
#         list(train_gain = classifier_loglikelihood(y_train = y_train, predictions = predictions$train, oob = control$classifier_oob)$gain,
#              test_gain = NULL)
#       }
#     # }
#   }
#   gain_frame$alpha <- alpha
#   gain_frame <- cbind(gain_frame, apply(sapply(alpha, f), 1, unlist))
#   
#   gain_frame
# }

#' Glasso cross validation function
#'
#' @importFrom stats sd
#' @inheritParams glasso_gain_function
#' @param folds folds
glasso_cross_validation_function <-  function(x, start, end, lambda, folds, control){
  
  # get parameters
  n_folds <- control$cv_inner_n_folds
  search_lambda_inner <- control$cv_inner_search_lambda
  lambda_inner_step <- control$cv_inner_lambda_step
  lambda_inner <- control$cv_inner_lambda
  
  stopifnot(length(folds) == nrow(x))
  
  folds_current <- folds[(start + 1) : end]
  
  # function that returns a matrix of dimension (n_folds, length(lambda)) with respective losses
  evaluate_fit <- function(lambda){
    
    # prepare matrix to store values in later
    out <- array(NA, dim = c(length(unique(folds)), length(lambda)) )
    
    # iterate over folds and lambdas
    for (l in 1 : length(lambda)){
      for (i in as.integer(unique(folds_current))){
        glasso_fit <- get_glasso_fit(x[(start + 1) : end, , drop = F][folds_current != i, , drop = F], lambda[l], control)
        out[i, l] <- sum(loglikelihood(x[(start + 1) : end, glasso_fit$inds, drop = F][folds_current == i, ,drop = F],
                                       glasso_fit$mu[glasso_fit$inds],
                                       glasso_fit$wi[glasso_fit$inds, glasso_fit$inds, drop = F]), na.rm = T)
        
      }
    }
    
    # return results
    out
  }
  
  
  if(search_lambda_inner){
    
    # if search_lambda_inner is true, we search for an optimal lambda by adjusting lambda0 in steps of size
    # lambda_inner_step until a valley is reached
    lambda <- c(1 / lambda_inner_step, 1, lambda_inner_step) * lambda
    
    loss <- evaluate_fit(lambda)
    loss_sum <- apply(loss, 2, sum, na.rm = T)
    loss_sd <- apply(loss, 2, stats::sd, na.rm = T)
    i <- which.min(loss_sum)
    
    if(i == 2){ # lambda already valley
      list(loss_array = loss[, i], lambda_opt = lambda[2], cv_loss = loss_sum[i])
    } else if (i == 1){
      while(TRUE){ # bad form, change this
        lambda <- c(lambda[1] / lambda_inner_step, lambda)
        loss_cur <- evaluate_fit(lambda[1])
        loss_sum <- c(sum(loss_cur, na.rm = T), loss_sum)
        loss_sd <- c(stats::sd(loss_cur, na.rm = T), loss_sd)
        loss <- cbind(loss_cur, loss)
        if(loss_sum[1] >= loss_sum[2]){
          return(list(loss_array = loss[, 2], lambda_opt = lambda[2], cv_loss = loss_sum[2]))
        }
      }
    } else {
      while (TRUE){
        k <- length(lambda)
        lambda <- c(lambda, lambda[k] * lambda_inner_step)
        loss_cur <- evaluate_fit(lambda[k + 1])
        loss_sum <- c(loss_sum, sum(loss_cur, na.rm = T))
        loss_sd <- c(loss_sd, stats::sd(loss_cur, na.rm = T))
        loss <- cbind(loss, loss_cur)
        if(loss_sum[k + 1] >= loss_sum[k]){
          return(list(loss_array = loss[, k], lambda_opt = lambda[k], cv_loss = loss_sum[k]))
        }
      }
    }
  } else {
    stopifnot(!is.null(lambda_inner))
    loss <- evaluate_fit(lambda_inner)
    loss_sum <- apply(loss, 2, sum)
    loss_sd <- apply(loss, 2, stats::sd)
    i <- which.min(loss_sum)
    return(list(loss_array = loss[, i], lambda_opt = lambda_inner[i], cv_loss = loss_sum[i]))
  }
}

# 
# cross_validate <- function(x, classifier, get_best_split, delta, lambda, segmentation, control){
#   
#   n <- nrow(x)
#   
#   # sample folds for cross validation
#   folds <- sample_folds(n, control$cv_classifier_n_folds, randomize = control$cv_classifier_randomize_folds)
#   
#   result <- data.frame(id = NULL, fold = NULL, gain = NULL, cv_loss = NULL, k = NULL, alpha = NULL)
#   trees <- list()
#   
#   for(i in 1 : control$cv_classifier_n_folds){
#     
#     # calculate a wbs tree for the current fold / trianing_data
#     trees[[i]] <- binary_segmentation(x[folds != i, ], get_best_split, delta, lambda, control = control, segmentation = segmentation)
#     
#     alpha <- get_splits(trees[[i]])
#     
#     for(j in 1 : length(alpha)){
#       
#       if(is.null(alpha[[j]])){
#         result <- rbind(result,
#                         list(fold = i, j = j, gain = 0, cv_loss = 0, k = 0)
#         )
#       } else {
#         train_test_data <- train_test_split(x, alpha[[j]], folds != i)
#         predictions <- classifier(x_train = train_test_data$x_train, y_train = train_test_data$y_train,
#                                   x_test = train_test_data$x_test, control = control)
#         loglikelihoods <- loglikelihood_estimate(train_test_data$y_train, predictions$train, train_test_data$y_test, predictions$test, oob = T)
#         
#         result <- rbind(result,
#                         data.frame(fold = i, j = j, gain = loglikelihoods$train_loss, cv_loss = loglikelihoods$test_loss,
#                                    k = length(alpha[[j]]))
#         )
#       }
#     }
#   }
#   lapply()
#   # get a sequence of change points, with the corresponding gamma when they were added
#   alpha <- as.matrix(trees[[i]]$Get(function(node) c(split_point = node$split_point, max_gain = node$max_gain), filterFun = function(x){!is.na(x[['max_gain']]) && !is.null(x[['max_gain']])}))# &&  x[['max_gain']] > 0}))
#   
#   # Extract the locations of change points and the corresponding max_gains. This is coded like this
#   # to work with trivial trees / single splits with positive gain
#   gamma <- alpha[2, ]
#   alpha <- alpha[1, ]
#   
#   alpha <- alpha[order(-gamma)]
#   gamma <- gamma[order(-gamma)]
#   
#   alpha <- alpha[gamma > 0]
#   gamma <- gamma[gamma > 0]
#   
#   
#   # loop over different change point configurations
#   for(j in seq_along(alpha)){
#     
#     change_points <- sort(alpha[1 : j])
#     k <- length(change_points) + 1 # number of segments
#     n_train <- sum(folds != i) # number of training observations in this fold
#     n_test <- sum(folds == i) # number of test observations in this fold
#     
#     # get vector of segment lengths and prepare target variable
#     segment_lengths <- c(change_points, n_train) - c(0, change_points)
#     y <- as.factor(rep.int(1 : k, times = segment_lengths))
#     
#     # get y_test TODO rewrite this more efficient
#     y_full <- rep(NA,n_test + n_train)
#     y_full[folds != i] <- y
#     y_full[1] <- 1
#     while(any(is.na(y_full))){
#       y_full[which(is.na(y_full))] <- y_full[which(is.na(y_full)) - 1]
#     }
#     y_test <- y_full[folds == i]
#     
#     # fit a probability machine
#     predictions <- classifier(x_train = x[folds != i, , drop = F], y_train = y, x_test = x[folds == i, , drop = F], control = control)
#     
#     if(ncol(predictions) == 1){
#       predictions <- cbind(predictions, 1 - predictions)
#     }
#     # rescale by dividing by the expected predictions under H0, i.e. segment_length / n & calculate log-likelihood
#     predictions <- matrix(n_train / segment_lengths, nrow = n_test, ncol = k, byrow = T) * predictions
#     predictions <- log(predictions[matrix(1 : k, nrow = n_test, ncol = k, byrow = TRUE) == y_test])# / rowMeans(predictions))
#     
#     loss <- rbind(loss, list(gamma = gamma[j], fold = i, loss = sum(predictions)))
#   }
# }
# #loss <- rbind(loss, list(gamma = max(loss$gamma) + 1, fold = 1 : n_folds, loss = -Inf))
# temp <- dcast(loss, gamma ~ fold, value.var = 'loss')
# temp[, -1] <- data.frame(lapply(temp[, -1], function(y) na.omit(y)[cumsum(!is.na(c(0, y[-length(y)])))]))
# temp_2 <- cbind(temp[, 1], rowMeans(temp[, -1]))
# i <- which.max(temp_2[,2])
# gamma <- mean(temp_2[, 1][pmax(1, c(i, i - 1))])
# 
# tree <- binary_segmentation(x, get_best_split, delta, lambda, control = control, segmentation = segmentation, gamma = gamma)
# tree$gamma <- gamma
# tree$cv_results <- list(reduced = temp_2, full = temp, trees = trees)
# class(tree) <- c(class(tree), 'cross_validated_binary_segmentation_tree')
# tree
# }
