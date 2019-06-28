#' get_gain_function_from_loss_function
#'
#' Returns a closure with formal arguments \code{split_point}
#'
#' @inheritParams hdcd
#' @param lambda If \code{loss_function} has argument \code{lambda}, then this is the standard value used
#' for its evaluation as long as no different lambda is supplied
get_gain_function_from_loss_function <- function(loss_function, lambda = NULL) {
    
    function(x, start, end, lambda = lambda) {
        
        stopifnot(!is.null(lambda))
        
        # calculate global loss once to compare to later
        global_loss <- loss_function(x[(start + 1):end, , drop = F], lambda)
        
        function(split_point) {
            
            # series of check to avoid unwanted behaviour
            stopifnot(start <= split_point)
            stopifnot(split_point <= end)
            
            global_loss - loss_function(x[(start + 1):split_point, , drop = F], 
                lambda = lambda) - loss_function(x[(split_point + 1):end, 
                , drop = F], lambda = lambda)
        }
    }
}
