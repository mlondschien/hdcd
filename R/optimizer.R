#' Section search optimisation algorithm (OBS, Optimistic Binary Segmentation)
#' 
#' Finds an (approximate) maximum of \code{gain_function} by intelligently evaluating at some of the \code{split_candidates}.
#' Works best if \code{split_candidates} is an interval.
#' 
#' @inheritParams hdcd
#' @param split_candidates An array of indices where gain function is (possibly) evaluated. 
#' @return An array of length \code{control$n_obs} with entries for gains at different splits
section_search <- function(gain_function, split_candidates, control) {
    
    # Check that gain_function is compatible
    if (!("split_point" %in% formalArgs(gain_function))) {
        stop("gain_function is not of the required form. Make sure that gain_function evaluated at x, start, end and lambda returns
         a closure with a single argument split_point.")
    }
    
    min_points <- control$section_search_min_points
    if (min_points < 5) {
        stop("section_search_min_points needs to be chosen >= 5")
    }
    n_obs <- control$n_obs
    stepsize <- control$section_search_stepsize
    
    stopifnot(!any(c(is.null(min_points), is.null(n_obs), is.null(stepsize))))
    stopifnot(!is.null(n_obs))
    
    gain <- rep(NA, n_obs)
    
    # helper function that selects a new midpoint in the smaller of the two
    # segments (cur_left, cur_middle] and (cur_middle, cur_right] and then
    # calls section_search_recursive
    select_next_midpoint <- function(cur_left, cur_middle, cur_right) {
        
        # If the current segment (cur_left, cur_right] is small, find the
        # maximum of the gain_function on that segment with a line_search
        if (cur_right - cur_left < min_points) {
            
            # select indices for which to calculate gains
            inds <- intersect(intersect((cur_left):(cur_right), split_candidates), 
                which(is.na(gain)))
            
            # Evaluate gain on all of the split_candidates within w_left : w_right
            # to search for maximum
            gain[inds] <<- sapply(inds, gain_function)
            return(gain)
        }
        
        stopifnot(cur_left < cur_middle)
        stopifnot(cur_middle < cur_right)
        
        # select new midpoint in smaller segment
        if (cur_right - cur_middle > cur_middle - cur_left) {
            w <- cur_right - ceiling((cur_right - cur_middle) * stepsize)
            section_search_recursive(cur_left, cur_middle, w, cur_right)
        } else {
            w <- cur_left + ceiling((cur_middle - cur_left) * stepsize)
            section_search_recursive(cur_left, w, cur_middle, cur_right)
        }
    }
    
    # Recursive function that finds the maximum of the gain_function on the
    # segment (cur_left, cur_right] via a binary search style algorithm.
    section_search_recursive <- function(cur_left, w_left, w_right, cur_right) {
        
        # When called, at least one of gain[w_left] and gain[w_right] is NA.
        # We then evaluate gain_function at this w (if it is a
        # split_candidate), in/decrease it by one if doing so is not possible
        # and try again. When w_right - w_left < min_points, the maximum is
        # found via a line_search.
        while (is.na(gain[w_left]) & w_right - w_left > 0) {
            if (w_left %in% split_candidates) {
                gain[w_left] <<- gain_function(w_left)
            }
            if (is.na(gain[w_left])) {
                w_left <- w_left + 1
            }
        }
        
        while (is.na(gain[w_right]) & w_right - w_left > 0) {
            if (w_right %in% split_candidates) {
                gain[w_right] <<- gain_function(w_right)
            }
            if (is.na(gain[w_right])) {
                w_right <- w_right - 1
            }
        }
        
        # if the gain function fails to evaluate at all w_left, w_right, just
        # evaluate w on the whole segment (ie do a line search) and return gain
        if (w_right <= w_left) {
            gain[cur_left:cur_right] <- sapply(cur_left:cur_right, gain_function)
            return(gain)
        }
        
        # if the gain is higher at w_left than at w_right, discard the segment
        # (w_right, cur_right]. Else discard the segment (cur_left, w_right]
        if (gain[w_left] >= gain[w_right]) {
            select_next_midpoint(cur_left, w_left, w_right)
        } else {
            select_next_midpoint(w_left, w_right, cur_right)
        }
    }
    
    cur_left <- min(split_candidates)
    cur_right <- max(split_candidates)
    
    w_left <- ceiling((cur_left + stepsize * cur_right)/(1 + stepsize))
    w_right <- floor((stepsize * cur_left + cur_right)/(1 + stepsize))
    
    gain <- section_search_recursive(cur_left, w_left, w_right, cur_right)
    
    best_split = ifelse(all(is.na(gain)), NA, which.max(gain))
    max_gain = ifelse(all(is.na(gain)), NA, max(gain, na.rm = T))
    
    list(gain = gain, best_split = best_split, max_gain = max_gain)
}

#' Line search optimisation algorithm
#' 
#' Finds an (exact) maximum of \code{gain_function} by evaluating it at each of \code{split_candidates}.
#' 
#' @inheritParams hdcd
#' @param split_candidates array of indices where gain function is evaluated. 
#' @return An array of length \code{control$n_obs} with entries for gains at all of \code{split_candidates}.
line_search <- function(gain_function, split_candidates, control) {
    
    # Check that gain_function is compatible
    if (!("split_point" %in% formalArgs(gain_function))) {
        stop("gain_function is not of the required form. Make sure that gain_function evaluated at x, start, end and lambda returns
         a closure with a single argument split_point.")
    }
    
    stopifnot(!is.null(control$n_obs))
    
    # Create array to save evaluated gains
    gain <- rep(NA, control$n_obs)
    
    # evaluate the gain_function at each of split_candidates and save
    # corresponding gains in gain
    gain[split_candidates] <- sapply(split_candidates, gain_function)
    
    best_split = ifelse(all(is.na(gain)), NA, which.max(gain))
    max_gain = ifelse(all(is.na(gain)), NA, max(gain, na.rm = T))
    
    list(gain = gain, best_split = best_split, max_gain = max_gain)
}
