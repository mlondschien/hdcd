#' Section search optimisation algorithm
#' 
#' Finds an (approximate) maximum of \code{gain_function} by intelligently evaluating at some of the \code{split_candidates}.
#' Works best if \code{split_candidates} is an interval.
#' 
#' @inheritParams hdcd
#' @param split_candidates An array of indices where gain function is (possibly) evaluated. 
#' @param n total number of observations
#' @param gain_function closure
#' @return A list with entries gain, max_gain and best_split
section_search <- function(gain_function, split_candidates, n, control){
  
  # get parameters
  min_points <- control$section_search_min_points
  
  if(min_points < 5){
    stop('section_search_min_points needs to be chosen >= 5')
  }
  
  stepsize <- control$section_search_stepsize
  
  stopifnot(!any(is.null(min_points), is.null(stepsize)))
  
  # initiate vector to store gains
  gain <- rep(NA, n)
  
  # helper function that selects a new midpoint in the smaller of the two segments (cur_left, cur_middle] 
  # and (cur_middle, cur_right] and then calls section_search_recursive
  select_next_midpoint <- function(cur_left, cur_middle, cur_right){
    
    # If the current segment (cur_left, cur_right] is small, find the maximum of the gain_function on that
    # segment with a line_search
    if(cur_right - cur_left < min_points){
      
      # select indices for which to calculate gains
      inds <- intersect(intersect((cur_left) : (cur_right), split_candidates), which(is.na(gain)))
      
      gain[inds] <<- sapply(inds, gain_function)
      return(gain)
    }
    
    stopifnot(cur_left < cur_middle)
    stopifnot(cur_middle < cur_right)
    
    # select new midpoint in smaller segment
    if (cur_right - cur_middle > cur_middle - cur_left){
      w <- cur_right - ceiling((cur_right - cur_middle) * stepsize)
      section_search_recursive(cur_left, cur_middle, w, cur_right)
    } else {
      w <- cur_left + ceiling((cur_middle - cur_left) * stepsize)
      section_search_recursive(cur_left, w, cur_middle, cur_right)
    }
  }
  
  # Recursive function that finds the maximum of the gain_function on the segment (cur_left, cur_right] via a
  # binary search style algorithm.
  section_search_recursive <- function(cur_left, w_left, w_right, cur_right){
    # When called, at least one of gain[w_left] and gain[w_right] is NA. We then evaluate gain_function at this w
    # (if it is a split_candidate), in/decrease it by one if doing so is not possible and try again.
    while(is.na(gain[w_left]) & w_right - w_left > 0){
      if(w_left %in% split_candidates){
        gain[w_left] <<- gain_function(w_left)
      }
      if(is.na(gain[w_left])){
        w_left <- w_left + 1
      }
    }
    
    while(is.na(gain[w_right]) & w_right - w_left > 0){
      if(w_right %in% split_candidates){
        gain[w_right] <<- gain_function(w_right)
      }
      if(is.na(gain[w_right])){
        w_right <- w_right - 1
      }
    }
    
    # if the gain function fails to evaluate at all w_left, w_right, just evaluate w on the whole segment (ie do a line search)
    # and return gain
    if(w_right <= w_left){
      # select indices for which to calculate gains
      inds <- intersect(intersect((cur_left) : (cur_right), split_candidates), which(is.na(gain)))
      
      # Evaluate gain on all of the split_candidates within w_left : w_right to search for maximum
      gain[inds] <<- sapply(inds, gain_function)
      
      return(gain)
    }
    
    # if the gain is higher at w_left than at w_right, discard the segment (w_right, cur_right]. Else discard the
    # segment (cur_left, w_right]
    if(gain[w_left] >= gain[w_right]){
      select_next_midpoint(cur_left, w_left, w_right)
    } else {
      select_next_midpoint(w_left, w_right, cur_right)
    }
  }
  
  cur_left <- min(split_candidates)
  cur_right <- max(split_candidates)
  
  w_left <- ceiling((cur_left + stepsize * cur_right) / (1 + stepsize))
  w_right <- floor((stepsize * cur_left + cur_right) / (1 + stepsize))
  
  gain <- section_search_recursive(cur_left, w_left, w_right, cur_right)
  
  best_split = ifelse(all(is.na(gain)), NA, which.max(gain))
  max_gain = ifelse(all(is.na(gain)), NA, max(gain, na.rm = T))
  
  list(gain = gain, best_split = best_split, max_gain = max_gain)
}

#' smooth section search optimisation algorithm
#' 
#' Finds an (approximate) maximum of \code{gain_function} by intelligently evaluating at some of the \code{split_candidates}.
#' Works best if \code{split_candidates} is an interval.
#' 
#' @inheritParams hdcd
#' @param split_candidates An array of indices where gain function is (possibly) evaluated. 
#' @param n total number of observations
#' @param gain_function closure
#' @return A list with entries gain, max_gain and best_split
smooth_section_search <- function(gain_function, split_candidates, n, control){
  
  # get parameters
  min_points <- control$section_search_min_points
  if(min_points < 5){
    stop('section_search_min_points needs to be chosen >= 5')
  }
  stepsize <- control$section_search_stepsize
  
  stopifnot(!any(is.null(min_points), is.null(stepsize)))
  
  # initiate vector to store gains
  gain <- rep(NA, n)
  
  # helper function that selects a new midpoint in the smaller of the two segments (cur_left, cur_middle] 
  # and (cur_middle, cur_right] and then calls section_search_recursive
  select_next_midpoint <- function(cur_left, cur_middle, cur_right){
    
    # If the current segment (cur_left, cur_right] is small, find the maximum of the gain_function on that
    # segment with a line_search
    if(cur_right - cur_left < min_points){
      
      # select indices for which to calculate gains
      inds <- intersect(intersect((cur_left) : (cur_right), split_candidates), which(is.na(gain)))
      
      if(length(inds) > 0){
        gain[inds] <<- gain_function(split_point = ceiling(stats::median(inds)), split_candidates = inds)[inds]
      }
      return(gain)
    }
    
    stopifnot(cur_left < cur_middle)
    stopifnot(cur_middle < cur_right)
    
    # select new midpoint in smaller segment
    if (cur_right - cur_middle > cur_middle - cur_left){
      w <- cur_right - ceiling((cur_right - cur_middle) * stepsize)
      section_search_recursive(cur_left, cur_middle, w, cur_right)
    } else {
      w <- cur_left + ceiling((cur_middle - cur_left) * stepsize)
      section_search_recursive(cur_left, w, cur_middle, cur_right)
    }
  }
  
  # Recursive function that finds the maximum of the gain_function on the segment (cur_left, cur_right] via a
  # binary search style algorithm.
  section_search_recursive <- function(cur_left, w_left, w_right, cur_right){
    # When called, at least one of gain[w_left] and gain[w_right] is NA. We then evaluate gain_function at this w
    # (if it is a split_candidate), in/decrease it by one if doing so is not possible and try again.
    while(is.na(gain[w_left]) & w_right - w_left > 0){
      if(w_left %in% split_candidates){
        gain[w_left] <<- gain_function(w_left)
      }
      if(is.na(gain[w_left])){
        w_left <- w_left + 1
      }
    }
    
    while(is.na(gain[w_right]) & w_right - w_left > 0){
      if(w_right %in% split_candidates){
        gain[w_right] <<- gain_function(w_right)
      }
      if(is.na(gain[w_right])){
        w_right <- w_right - 1
      }
    }
    
    # if the gain function fails to evaluate at all w_left, w_right, just evaluate w on the whole segment (ie do a line search)
    # and return gain
    if(w_right <= w_left){
      # select indices for which to calculate gains
      inds <- intersect(intersect((cur_left) : (cur_right), split_candidates), which(is.na(gain)))
      
      if(length(inds) > 0){
        #cat('evaluate gain function on', inds, '\n')
        gain[inds] <<- gain_function(split_point = ceiling(stats::median(inds)),  split_candidates = inds)[inds]
      }
      return(gain)
    }
    
    # if the gain is higher at w_left than at w_right, discard the segment (w_right, cur_right]. Else discard the
    # segment (cur_left, w_right]
    if(gain[w_left] >= gain[w_right]){
      select_next_midpoint(cur_left, w_left, w_right)
    } else {
      select_next_midpoint(w_left, w_right, cur_right)
    }
  }
  
  cur_left <- min(split_candidates)
  cur_right <- max(split_candidates)
  
  w_left <- ceiling((cur_left + stepsize * cur_right) / (1 + stepsize))
  w_right <- floor((stepsize * cur_left + cur_right) / (1 + stepsize))
  
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
#' @param split_candidates An array of indices where gain function is (possibly) evaluated. 
#' @param n total number of observations
#' @param gain_function closure
#' @return A list with entries gain, max_gain and best_split
line_search <- function(gain_function, split_candidates, n, control){
  
  stopifnot(!is.null(control$n))
  
  # Create array to save evaluated gains
  gain <- rep(NA, n)
  
  # evaluate the gain_function at each of split_candidates and save corresponding gains in gain
  gain[split_candidates] <- sapply(split_candidates, gain_function)
  
  best_split = ifelse(all(is.na(gain)), NA, which.max(gain))
  max_gain = ifelse(all(is.na(gain)), NA, max(gain, na.rm = T))
  
  list(gain = gain, best_split = best_split, max_gain = max_gain)
}


#' Two step search optimisation algorithm
#' 
#' Find an (approximate) maximum of \code{gain_function} by repeatedly fitting and evaluating.
#' 
#' @inheritParams hdcd
#' @param split_candidates An array of indices where gain function is (possibly) evaluated. 
#' @param n total number of observations
#' @param gain_function closure
#' @return A list with entries gain, max_gain, best_split, permutation_test and pval
two_step_search <- function(gain_function, split_candidates, n, control){
  
  #prepare output
  #res <- list(gain_left = rep(NA, n),
  #            gain_mid = rep(NA, n),
  #            gain_right = rep(NA, n),
  #            gain = rep(NA, n), best_split = NA, max_gain = NA)
  res <- list()
  #split <- round(median(split_candidates)) # first first split
  split <- round(stats::quantile(split_candidates, c(0.25, 0.5, 0.75)))
  
  temp_left <- gain_function(split_point = split[1], split_candidates = split_candidates)
  temp_mid <- gain_function(split_point = split[2], split_candidates = split_candidates)
  temp_right <- gain_function(split_point = split[3], split_candidates = split_candidates)
  
  res$gain_left <- temp_left$gain
  res$gain_mid <- temp_mid$gain
  res$gain_right <- temp_right$gain
  
  if(all(is.na(c(temp_left$max_gain, temp_mid$max_gain, temp_right$max_gain)))){
    res$best_spit <- NA
    res$max_gain <- NA
    return(res)
  }
  
  split <- c(temp_left$best_split, temp_mid$best_split, temp_right$best_split)[which.max(c(temp_left$max_gain, temp_mid$max_gain, temp_right$max_gain))]
  
  temp <- gain_function(split_point = split, split_candidates = split_candidates)
  res$gain <- temp$gain
  res$max_gain <- temp$max_gain
  res$best_split <- temp$best_split
  
  permutation_test <- cbind(temp_left$permutation_test,
                            temp_mid$permutation_test,
                            temp_right$permutation_test,
                            temp$permutation_test)
  res$permutation_test <- apply(permutation_test, 1, max)
  
  res$pval <- (1 + sum(res$permutation_test >= res$max_gain)) / (1 + length(res$permutation_test))
  
  if(is.na(temp$max_gain)){
    res$best_spit <- NA
    res$max_gain <- NA
    return(res)
  }
  
  res$best_split <- temp$best_split
  res$max_gain <- temp$max_gain
  
  res
}