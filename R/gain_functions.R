#' Closure generating function to calculate gains when splitting and learning a glasso model
#' 
#' @param control an object of type \code{hdcd_control} returned by \link{hdcd_control}
#' 
#' @return A closure with parameters \code{x}, \code{start} and \code{end}, that when evaluated
#' will itself return a closure with parameter \code{split_point}. This calculates the gain when
#' splitting the segment \code{(start, end]} of \code{x} at \code{split_point}.
#' @export
get_glasso_gain_function <- function(control){
  
  function(x, start, end, lambda) {
    
    if(is.null(lambda)){
      stop('Please supply a value for lambda to the glasso gain function')
    }
    
    # get a global glasso fit to compare gains to.
    fit_global <- get_glasso_fit(x[(start + 1) : end, , drop = F], lambda = lambda, control = control)
    
    function(split_point){
      
      fit_left <- get_glasso_fit(x[(start + 1) : split_point, , drop = F], lambda = lambda, control = control)
      fit_right <- get_glasso_fit(x[(split_point + 1) : end, , drop = F], lambda = lambda, control = control)
    
      # only use variables for calculation of the gains curve that are used on left and right segment.
      # See curve smoothing section in paper
      (sum(loglikelihood(x[(start + 1) : split_point, fit_left$inds, drop = F],
                        fit_global$mu[fit_left$inds, drop = F],
                        fit_global$wi[fit_left$inds, fit_left$inds, drop = F])) +
            sum(loglikelihood(x[(split_point + 1) : end, fit_right$inds, drop = F],
                          fit_global$mu[fit_right$inds, drop = F],
                          fit_global$wi[fit_right$inds, fit_right$inds, drop = F])) -
            sum(loglikelihood(x[(start + 1) : split_point, fit_left$inds, drop = F],
                          fit_left$mu[fit_left$inds, drop = F],
                          fit_left$wi[fit_left$inds, fit_left$inds, drop = F])) -
            sum(loglikelihood(x[(split_point + 1) : end, fit_right$inds, drop = F],
                          fit_right$mu[fit_right$inds, drop = F],
                          fit_right$wi[fit_right$inds, fit_right$inds, drop = F] ))) / nrow(x)
    }
  }
}