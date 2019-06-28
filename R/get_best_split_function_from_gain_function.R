#' Get a best split function from a gain function
#' 
#' Returns a closure with formal arguments \code{x}, \code{start}, \code{end} and \code{split_candidates} that finds the
#' maximizer of the \code{gain_function} given data \code{x}, \code{start} and \code{end} on \code{split_candidates}
#' using the optimizer \code{optimizer}.
#' 
#' @param gain_function A function with formal arguments \code{x}, \code{start}, \code{end} and \code{lambda} that returns
#'  a closure with argument \code{split_point}, that returns the gain after splitting the segment (\code{start}, \code{end}]
#'  at \code{split_point} given data \code{x} and tuning parameter \code{lambda}.
#' @return A function with formal arguments \code{x}, \code{start}, \code{end}, \code{split_candidates} and \code{lambda} that
#'  uses the optimizer specified to search for a maximum of the gain_function on the \code{split_candidates} given a segment
#'  (\code{start}, end], data \code{x} and a tuning parameter \code{lambda} and returns a list with arguments \code{gain}, an
#'  array of length end - start with evaluations of the gain function and \code{best_split}.
get_best_split_function_from_gain_function <- function(gain_function, optimizer, 
    control) {
    
    stopifnot(optimizer %in% c("section_search", "line_search"))
    
    # Stop with error if gain_function does not take required arguments
    if (!all(c("x", "start", "end") %in% formalArgs(gain_function))) {
        stop("gain_function is not of the required form. Make sure that gain_function takes formal arguments x, start, end, lambda.")
    }
    
    # Return closure that estimates / calculates the location of the best
    # split (with maximum gain) within (start, end]
    function(x, start, end, split_candidates, lambda) {
        
        # Sequence of checks to avoid unwanted behaviour
        stopifnot(start >= 0)
        stopifnot(end <= nrow(x))
        stopifnot(start < min(split_candidates))
        stopifnot(end > max(split_candidates))
        
        # Apply optimizer
        if (optimizer == "section_search") {
            res <- section_search(gain_function(x = x, start = start, end = end, 
                lambda = lambda), split_candidates, control)
        } else if (optimizer == "line_search") {
            res <- line_search(gain_function(x = x, start = start, end = end, 
                lambda = lambda), split_candidates, control)
        } else {
            stop("Make sure that optimizer is one of 'section_search', 'line_search' ")
        }
        
        # In some situations (e.g. with missing values) the gain function
        # cannot be evaluated at any of the split_candidates. In such a
        # situation NA is returned as best_split and binary_segmentation is
        # stopped.
        if (is.na(res$best_split)) {
            warning("the gain_function could not be evaluated at any of the split_candidates. NA is returned and binary_segmentation
              stopped for segment (", 
                start, ", ", end, "]", sep = "")
        }
        
        res
    }
}

