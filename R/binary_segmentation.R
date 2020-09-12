#' Binary Segmentation
#'
#' Search for change points with Binary Segmentation, Wild Binary Segmentation or Seeded Binary Segmentation
#'
#' @param x An n x p data matrix
#' @param get_best_split A function with arguments \code{x}, \code{start}, \code{end}, \code{split_candidates} and \code{lambda}, that returns a list with
#' entries \code{gain}, \cope{max_gain}, \code{best_split} and possibly \code{permutation_test} and \code{pval}.
#' @param cross_validation_function A function with arguments \code{x}, \code{start}, \code{end}, \code{lambda}, \code{folds} and \code{control}
#' @inheritParams hdcd
#' @import data.table
#' @import data.tree
binary_segmentation <-
  function(x,
           get_best_split,
           delta,
           lambda = 0,
           gamma = 0,
           segmentation = 'BS',
           cross_validation_function = NULL,
           control = hdcd_control()) {
    # to remove warnings in CMD check
    start <- end <- max_gain <- permutation_test <- NULL
    
    n <- nrow(x)
    minimal_segment_length <- ceiling(delta * n)

    # get a data.frame of start - end pairs for which gains and optimal splits will be calculated.
    # Only relevant for SBS or WBS, will contain one segment (0, n] for BS
    segments <- data.table::as.data.table(
      draw_segments(
        n = n,
        n_segments = control$wbs_n_segments,
        delta = delta,
        segmentation = segmentation,
        alpha = control$sbs_alpha
      )
    )
    
    # add column with lists of split candidates
    segments$split_candidates <-
      list(mapply(function(start, end) {
        list((start + minimal_segment_length):(end - minimal_segment_length))
      }, segments$start, segments$end))

    # calculate gain, best split, ... for all segments and add to segments
    gains <-
      mapply(
        get_best_split,
        start = segments$start,
        end = segments$end,
        split_candidates = segments$split_candidates,
        MoreArgs = list(x = x, lambda = lambda)
      )
    segments <- cbind(segments, t(gains))
    
    # Create tree structure to save binary segmentation output. One node for each segment,
    # with subsegments due to splitting as children.
    tree <-
      data.tree::Node$new(paste('(', 0, ' ', n, ']', sep = ''),
                          start = 0,
                          end = n)
    class(tree) <- c("binary_segmentation_tree", class(tree))
    # Store estimation parameters in root node.
    tree$lambda <- lambda
    tree$delta <- delta
    tree$n <- n
    
    # In case a cross_validation_function that returns some segment loss, evaluate on the
    # root node.
    if (!is.null(cross_validation_function)) {
      tree$folds <-
        sample_folds(
          tree$n,
          control$glasso_cv_inner_n_folds,
          randomize = control$glasso_cv_inner_randomize_folds
        )
      temp <-
        cross_validation_function(
          x,
          start = 0,
          end = n,
          lambda = lambda,
          folds = tree$folds,
          control = control
        )
      if (is.null(temp$cv_loss) | is.na(temp$cv_loss) | is.null(temp$lambda_opt) | is.na(temp$lambda_opt)) {
        stop(
          "cross_validation_function is not of the required form. Make sure
          cross_validation_function returns a list with
          attributes cv_loss and lambda_opt"
        )
      }
      #Store results of cross_validation_function in tree
      tree$cv_loss <- temp$cv_loss
      tree$lambda <- temp$lambda_opt
      rm(temp)
      }
    
    # Recursive function that creates the binary segmentation tree.
    # If applied on a node where some stopping condition is reached,
    # returns NA
    binary_segmentation_recursive <- function(node) {
      # stop if segment is not long enough to allow for splitting
      segment_length <- node$end - node$start
      if (segment_length < 2 * minimal_segment_length) {
        return(NA)
      }

      # If there is no entry in segments corresponding to the current segment, add a row and calculate
      # the gain and best split. In case of binary segmentation, the resulting table segments will
      # comprise of exactly those segments found while splitting.
      if (!any(segments$start == node$start &
               segments$end == node$end)) {
        split_candidates <-
          (node$start + minimal_segment_length):(node$end - minimal_segment_length)

        temp <-
          get_best_split(
            x = x,
            start = node$start,
            end = node$end,
            split_candidates = split_candidates,
            lambda = node$lambda
          )

        segments <<- rbind(segments,
                           c(
                             list(
                               start = node$start,
                               end = node$end,
                               split_candidates = list(split_candidates)
                             ),
                             lapply(temp, function(i) {
                               if (length(i) == 1) {
                                 i
                               } else{
                                 list(i)
                               }
                             })
                           ))
      }
  
      best_segment <-
        segments[start >= node$start &
                   end <= node$end, ][unlist(max_gain) == max(unlist(max_gain))][1]
      
      node$gain <- unlist(best_segment$gain)
      node$split_point <- unlist(best_segment$best_split)
      node$max_gain <- unlist(best_segment$max_gain)

      if (control$permutation_test) {
        node$permutation_test <-
          apply(matrix(unlist(segments[start >= node$start &
                                         end <= node$end, permutation_test]), nrow = control$permutation_test_n),
                1,
                max)
        node$pvalue <-
          (1 + sum(node$permutation_test > node$max_gain)) / (1 + length(node$permutation_test))
        if (node$pvalue >= control$permutation_test_pvalue) {
          return(NA)
        }
      }
      
      # stop if no best split can be found (e.g. when the gain is nonpositive for each split)
      if (is.na(node$split_point) || node$max_gain <= gamma) {
        return(NA)
      }
      
      stopifnot(node$start < node$split_point &
                  node$split_point < node$end)
      
      # Create left child
      child_left <- node$AddChild(
        paste('(', node$start, ' ', node$split_point, ']', sep = ''),
        start = node$start,
        end = node$split_point,
        lambda = node$lambda
      )
      
      # Create right child
      child_right <- node$AddChild(
        paste('(', node$split_point, ' ', node$end, ']', sep = ''),
        start = node$split_point,
        end = node$end,
        lambda = node$lambda
      )
      
      class(child_left) <-
        c("binary_segmentation_tree", class(child_left))
      class(child_right) <-
        c("binary_segmentation_tree", class(child_right))
      
      # If cross_validation_function is supplied, calculate cv_losses of both subsegments and calculate
      # inner cross-validated gain
      if (!is.null(cross_validation_function)) {
        temp_left <-
          cross_validation_function(
            x,
            start = node$start,
            end = node$split_point,
            lambda = node$lambda,
            folds = node$root$folds,
            control = control
          )
        temp_right <-
          cross_validation_function(
            x,
            start = node$split_point,
            end = node$end,
            lambda = node$lambda,
            folds = node$root$folds,
            control = control
          )
        node$cv_improvement <-
          node$cv_loss - temp_left$cv_loss - temp_right$cv_loss
        child_left$cv_loss <- temp_left$cv_loss
        child_left$lambda <- temp_left$lambda_opt
        child_right$cv_loss <- temp_right$cv_loss
        child_right$lambda <- temp_right$lambda_opt
        rm(temp_left)
        rm(temp_right)
        if (node$cv_improvement < 0) {
          return(NA)
        }
      }
      
      binary_segmentation_recursive(child_left)
      binary_segmentation_recursive(child_right)
    }
    
    binary_segmentation_recursive(tree)
    tree$segments <- segments
    
    tree
    }
