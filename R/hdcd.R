#' Random Forest change point detection
#'
#' Find (non-parametric) breaks in the distribution of a time series using Random Forest Classifiers
#' 
#'@inheritParams hdcd
RFcd <- function(x, delta = 0.1, control = hdcd_control()){
  
  tree <- hdcd(x = x,
               delta = delta,
               gamma = 0,
               method = 'RFcd',
               optimizer = 'two_step_search',
               segmentation = 'SBS',
               control = control)
  
  tree$change_points <- unname(tree$Get('split_point', filterFun = function(x) !is.null(x[['pvalue']]) && x[['pvalue']] < 0.02))
  
  tree
}

#' k-Nearest Neighbor based change point detection
#' 
#' Find (non-parametric) breaks in the distribution of a time series using k-Nearest Neigbor Classifiers
#' 
#' @inheritParams RFcd
kNNcd <- function(x, delta = 0.1, control = hdcd_control()){
  
  tree <- hdcd(x = x,
               delta = delta,
               gamma = 0,
               method = 'kNNcd',
               optimizer = 'line_search',
               segmentation = 'SBS',
               control = control)
  
  tree$change_points <- unname(tree$Get('split_point', filterFun = function(x) !is.null(x[['pvalue']]) && x[['pvalue']] < 0.05))
  
  tree
}

#' High-dimensional change point detection
#' 
#' Find breaks in distribution in a possiblity high-dimensional time series
#' 
#' @param x A matrix with observations in rows
#' @param delta Minimal relative segment length, defaults to 0.1
#' @param gamma Minimal required gain to split in the Binary Segmentation Procedure. Defaults to 0.
#' @param lambda Tuning parameter passed to the method.
#' @param method Method to be applied for change point detection. One of \code{glassocd}, \code{RFcd}, \code{kNNcd},
#'  \code{custom_classifier}, \code{custom_gain_function}, \code{custom_best_split_function}.
#' @param optimizer Optimization method to be used to find an (approximate) of the maximum of the gain function.
#' One of \code{line_search} for BS or similar, \code{section_search} for OBS or similar or \code{two_step_search}.
#' @param segmentation Segmentation method to be used. One of \code{BS}, \code{SBS} or \code{WBS}.
#' @param control Control parameter as returned by \link{hdcd_control}
#' @export
hdcd <- function(x,
                 delta = 0.1, 
                 gamma = 0,
                 lambda = NULL,
                 method = NULL,
                 optimizer = 'line_search',
                 segmentation = 'BS',
                 control = hdcd_control()
                 ){
  
  # First make sure x is of matrix form
  if(!is.matrix(x)){
    x <- as.matrix(x)
    warning('x has been coerced to matrix by hdcd')
  }

  control$n = length(x)
  
  if(method == 'RFcd'){
    get_best_split <- classifier_best_split_function(classifier = random_forest, optimizer = optimizer, control = control)
    cross_validation_function <- NULL
  } else if (method == 'kNNcd'){
    get_best_split <- best_split_function_from_gain_function(kNN_gain_function(x, control), optimizer, control)
    cross_validation_function <- NULL
  } else if (method == 'glasso'){
    if(is.null(lambda)){
      cov_mat <- get_cov_mat(x, control$glasso_NA_method)$mat
      lambda <- max(abs(cov_mat[upper.tri(cov_mat)])) / 10
      warning('Lambda for glassocd set by asymptotic theory to ', lambda)
    }
    if(control$glasso_cv_inner){
      cross_validation_function <- glasso_cross_validation_function
    } else { 
      cross_validation_function <- NULL
    }
    get_best_split <- best_split_function_from_gain_function(glasso_gain_function, optimizer, control)
  } else {
    stop("method should be one of 'glasso', 'RFcd' or 'kNNcd'. Got %s", method)
  }
  
  tree <- binary_segmentation(x = x, gamma = gamma,
                              get_best_split = get_best_split, 
                              delta = delta, 
                              lambda = lambda,
                              segmentation = segmentation, 
                              cross_validation_function = cross_validation_function,
                              control = control)
  
  tree
}

#add glassocd

#' hdcd cv
# cv_hdcd <- function(x,
#                     method,
#                     optimizer,
#                     segmentation,
#                     delta = 0.1,
#                     lambda = 0,
#                     control = hdcd_control()){
#   
#   n <- nrow(x)
#   
#   classifier <- switch(method,
#                        inbag_random_forest = inbag_random_forest,
#                        random_forest = random_forest,
#                        extra_trees = extra_trees,
#                        logistic_regression = logistic_regression,
#                        bagged_trees = bagged_trees,
#                        nnet = nnet,
#                        kNN = kNN
#   )
#   
#   if(method == 'kNN'){
#     get_best_split <- best_split_function_from_gain_function(kNN_gain_function(x, control), optimizer, control)
#   } else {
#     get_best_split <- classifier_best_split_function(classifier = classifier, optimizer = optimizer, control = control)
#   }
#   
#   tree <- binary_segmentation(x = x, gamma = 0,
#                               get_best_split = get_best_split, 
#                               delta = delta, 
#                               lambda = lambda,
#                               segmentation = segmentation,
#                               control = control)
#   tree$alpha <- classifier_model_select(tree = tree, x = x, classifier = random_forest, control = control)
#   tree$alpha$train_gain <- tree$alpha$train_gain / nrow(x)
#   tree$alpha_reduced <- data.table::as.data.table(tree$alpha)[, .SD[train_gain == max(train_gain)], by = k]
#   
#   # sample folds for cross validation
#   folds <- sample_folds(n, control$cv_classifier_n_folds, randomize = control$cv_classifier_randomize_folds)
#   
#   result <- data.table::data.table(id = 0, fold = 0, train_gain = 0, test_gain = 0, k = 0, alpha = 0)
#   trees <- list()
#   
#   for(i in 1 : control$cv_classifier_n_folds){
#     
#     # calculate a tree for the current fold / trianing_data
#     trees[[i]] <- binary_segmentation(x = x[folds != i, ], get_best_split = get_best_split, delta = delta, lambda = lambda, control = control, segmentation = segmentation)
#     
#     temp <- classifier_model_select(x = x, tree = trees[[i]], classifier =  classifier, control = control, fold = (folds != i))
#     temp$train_gain <- temp$train_gain / sum(folds != i)
#     temp$fold <- i
#     result <- rbind(result, temp)
#   }
#   
#   result <- result[fold != 0, .SD[train_gain == max(train_gain)], by = .(fold, k)]
#   # recover()
#   result[, gamma := 0][, keep := FALSE]
#   # for(i in 1 : control$cv_classifier_n_folds){
#   #   j <- result[fold == i & k == 0, id]
#   #   result[fold == i & k == 0, gamma := Inf]
#   #   while(any(is.na(result[fold == i, gamma]))){
#   #     gain_0 <- result[fold == i & id == j, train_gain]
#   #     k_0 <- result[fold == i & id == j, k]
#   #     
#   #     j <- result[fold == i & is.na(gamma)][ (gain_0 - train_gain) / (k_0 - k) == max((gain_0 - train_gain) / (k_0 - k), na.rm = T), id]
#   #     stopifnot(length(j) == 1)
#   #     result[fold == i & id == j, gamma := (gain_0 - train_gain) / (k_0 - k)]
#   #   }
#   # }
#   #recover()
#   for(i in 1 : control$cv_classifier_n_folds){
#     gain_0 <- k_0 <- 0
#     j_0 <- result[fold == i & k == 0, id]
#     #result[fold == i & id == j_0, keep := T]
#     gamma_0 <- NULL
#     while(nrow(result[fold == i & k > k_0 & train_gain > gain_0])>0){
#       j <- result[fold == i & k > k_0 & train_gain > gain_0][ (gain_0 - train_gain) / (k_0 - k) == max((gain_0 - train_gain) / (k_0 - k)), id]
#       result[fold == i & id == j, keep := T]
#       gamma_0 <- result[fold == i & id == j, (gain_0 - train_gain) / (k_0 - k)]
#       result[fold == i & id == j, gamma := gamma_0]#mean(c(gamma_0, gamma_1))]
#       #gamma_0 <- gamma_1
#       
#       k_0 <- result[fold == i & id == j, k]
#       gain_0 <- result[fold == i & id == j, train_gain]
#       j_0 <- j
#     }
#   }
#   tree$cv_result <- result
#   result <- result[keep == T, ]
#   #recover()
#   #stopifnot(!all(result$gamma == 0))
#   if(nrow(result) > 0){
#     casted_results <- data.frame(data.table::dcast(result, gamma ~ fold, value.var = 'test_gain'))#, fun.aggregate = mean)) ## REPAIR THIS
#     #casted_results <- casted_results[order(-casted_results$gamma), ]
#     casted_results[, -1] <- lapply(casted_results[, -1], function(y) na.omit(y)[cumsum(!is.na(c(0, y[-length(y)])))])
#   
#     casted_results_reduced <- data.frame(gamma = casted_results[, 1], gain = apply(casted_results[, -1, drop = F], 1, mean))
#   
#     index_opt <- which.max(casted_results_reduced$gain)# == max(casted_results_reduced$gain))
#   #gamma_1sd <- max(0, casted_results_reduced$gamma[ casted_results_reduced$gain + casted_results_reduced$sd[index_opt] > casted_results_reduced$gain[index_opt]])
#   #gamma_05sd <- max(0, casted_results_reduced$gamma[ casted_results_reduced$gain + 0.5 * casted_results_reduced$sd[index_opt] > casted_results_reduced$gain[index_opt]])
#     gamma_opt <- mean(casted_results_reduced$gamma[index_opt])
#     tree$casted_results <- casted_results
#     tree$casted_results_reduced <- casted_results_reduced
#   } else {
#     gamma_opt <- 0
#   } 
#   
#   tree$gamma_opt <- gamma_opt
#   #tree$gamma_1sd <- gamma_1sd
#   tree$cv_trees <- trees
#   
#   tree$result <- tree$alpha_reduced[, gain_gamma_opt :=  train_gain - k * gamma_opt]#[, gain_gamma_1sd := train_gain - k * gamma_1sd][, gain_gamma_05sd := train_gain - k* gamma_05sd]
#   #tree$cv_result <- result
#   
#   tree$change_points <- unlist(tree$result[gain_gamma_opt == max(gain_gamma_opt), alpha])
#   if(is.list(tree$change_points)) tree$change_points <- unlist(tree$change_points)
#   
#   tree
# }

