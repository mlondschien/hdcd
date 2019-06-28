#' Sample Folds
#' 
#' @param n number of observations
#' @param k number of folds
#' @param randomize should folds be randomized? equispaced folds are returned if randomize = F
#' @export
sample_folds <- function(n, k, randomize = FALSE) {
    if (k == 1) {
        as.factor(rep(0, n))
    } else if (randomize) {
        random_draw <- runif(n)
        k_quantiles <- quantile(random_draw, 0:k/k)
        cut(random_draw, k_quantiles, labels = 1:k, include.lowest = TRUE)
    } else {
        as.factor(rep(1:k, ceiling(n/k))[1:n])
    }
}

#' print.binary_segmentation_tree
#'
#' S3 method for printing a binary_segmentation_tree object.
#'
#' Decorate the print method of the data.tree package to see more details at each node.
#'
#' @param x A data.tree node.
#' @param ... Further arguments passed to print generic.
#' @export
print.binary_segmentation_tree <- function(x, ...) {
    
    # temp <- list('max_gain') if(!is.null(x$model_selection_statistic)){
    # temp[[length(temp) + 1]] <- 'model_selection_statistic' }
    NextMethod(generic = NULL, object = NULL, "split_point", "model_selection_statistic", 
        "cv_loss", "cv_improvement", "relative_cv_improvement", "lambda", 
        ...)
}

#' Logarithmically Scaled Sequence Generation
#' 
#' Generates a logarithmically scaled sequence
#' 
#' @param from the starting value of the sequence
#' @param to the end value of the sequence
#' @param length.out the length of the sequence
#' @export
log_space <- function(from, to, length.out) {
    exp(seq(from = log(from), to = log(to), length.out = length.out))
}

#' Get Change Points from a binary_segmentation_tree
#'
#' Utility function to get the change points with high enough value for some variable from a binary_segmentation_tree
#'
#' @param tree An object of class \strong{binary_segmentation_tree}
#' @param variable Name of the variable with respect to which should be pruned. Usually one of \code{cv_improvement} or
#' \code{max_gain}.
#' @param value Value until which the tree will be pruned
#' @export
#' @return A vector with the sorted change points.
get_change_points_from_tree <- function(tree, variable = "cv_improvement", value = 0) {
    
    alpha <- tree$Get("split_point", filterFun = function(x) {
        !is.na(x[[variable]]) && !is.null(x[[variable]]) && x[[variable]] > 
            value
    })
  
    data.frame(change_points = sort(unname(alpha)))
}

#'Rand type performance indices
#'
#' Calculate Rand type performance indices for two sets of change points. Typically one
#' of them will be the oracle estimate. See clues package for more details.
#'
#' @param cpts_a A sequence of changepoints.
#' @param cpts_b A sequence of changepoints.
#' @param n Total size of dataset from which both change point estimates originate.
#' @importFrom clues adjustedRand
#' @return Returns a vector of the index values.
#' @export
#'
#' @examples
#' compare_change_points(c(20, 50), c(30, 70), 100)
compare_change_points <- function(cpts_a, cpts_b, n) {
    cpts_a <- sort(cpts_a)[!duplicated(sort(cpts_a))]
    cpts_b <- sort(cpts_b)[!duplicated(sort(cpts_b))]
    
    MarkGroupings <- function(cpts) {
        diffs <- c(cpts, n) - c(0, cpts)
        rep(1:length(diffs), diffs)
    }
    
    clues::adjustedRand(MarkGroupings(cpts_a), MarkGroupings(cpts_b))
}
