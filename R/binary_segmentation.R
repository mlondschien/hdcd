#' binary_segmentation
#' 
#' @param x An n x p data matrix
#' @inheritParams hdcd
binary_segmentation <- function(x, get_best_split, delta, lambda, model_selection_function = NULL, 
    cross_validation_function = NULL, control = hdcd_control()) {
    
    control$n_obs <- n_obs <- nrow(x)
    
    # Create tree structure to save binary segmentation output. Each
    # segment corresponds to a node.
    tree <- data.tree::Node$new(paste("(", 0, " ", n_obs, "]", sep = ""), 
        start = 0, end = n_obs)
    class(tree) <- c("binary_segmentation_tree", class(tree))
    
    # Store estimation parameters in root node.
    tree$lambda <- lambda
    tree$delta <- delta
    tree$n_obs <- n_obs
    
    # If a cross_validation_function is supplied, calculate cv_error of the
    # segment (0, n_obs].
    if (!is.null(cross_validation_function)) {
        
        tree$folds <- sample_folds(tree$n_obs, control$cv_inner_n_folds, 
            randomize = control$cv_inner_randomize_folds)
        temp <- cross_validation_function(x, start = 0, end = tree$n_obs, 
            lambda = tree$lambda, folds = tree$folds)
        
        # Stop with error if cross_validation_function does not provide
        # necessary output
        if (is.null(temp$cv_loss) | is.null(temp$lambda_opt)) {
            stop("cross_validation_function is not of the required form. Make sure cross_validation_function returns a list with
           attributes cv_loss and lambda_opt")
        }
        
        # Store results of cross_validation_function in tree
        tree$cv_loss <- temp$cv_loss
        tree$lambda <- temp$lambda_opt
        rm(temp)
    }
    
    # Recursive function that creates the binary segmentation tree
    binary_segmentation_recursive <- function(node) {
        
        # stop if segment is not long enough to allow for splitting
        segment_length <- node$end - node$start
        if (segment_length/node$root$n_obs < 2 * node$root$delta) {
            return(NA)
        }
        
        # Set split_candidates to be points within (start, end] with distance
        # of at least delta * n_obs to the boundary
        split_candidates <- (node$start + ceiling(node$root$n_obs * node$root$delta)):(node$end - 
            ceiling(node$root$n_obs * node$root$delta))
        
        # obtain best split via get_best_split function
        temp <- get_best_split(x = x, start = node$start, end = node$end, 
            split_candidates = split_candidates, lambda = node$lambda)
        
        # Stop with error if get_best_split is not of the required form
        if (is.null(temp$gain) | is.null(temp$best_split)) {
            stop("get_best_split if not of the required form. Make sure get_best_split returns a list with attributes gain
           and best_split")
        }
        
        
        # Save output from get_best_split in node
        node$gain <- temp$gain
        node$split_point <- temp$best_split
        node$max_gain <- temp$max_gain
        
        # stop if no best split can be found (e.g. when the gain is nonpositive
        # for each split)
        if (is.na(node$split_point)) {
            return(NA)
        } else if (node$max_gain <= 0) {
            return(NA)
        }
        
        # Create left child
        child_left <- node$AddChild(paste("(", node$start, " ", node$split_point, 
            "]", sep = ""), start = node$start, end = node$split_point, 
            lambda = node$lambda)
        
        # Create right child
        child_right <- node$AddChild(paste("(", node$split_point, " ", 
            node$end, "]", sep = ""), start = node$split_point, end = node$end, 
            lambda = node$lambda)
        
        class(child_left) <- c("binary_segmentation_tree", class(child_left))
        class(child_right) <- c("binary_segmentation_tree", class(child_right))
        
        
        # if a model_selection_function is supplied, test for significance of
        # split
        if (!is.null(model_selection_function)) {
            temp <- model_selection_function(x, node$start, node$split_point, 
                node$end)
            
            # Stop with error if model_selection_function is not of required form
            if (is.null(temp$statistic) | is.null(temp$is_significant)) {
                stop("model_selection_function is not of the required form. Make sure model_selection_function returns a list with
             attritbutes statistic and is_significant")
            }
            
            # Save output of model_selection_function in node
            node$model_selection_statistic <- temp$statistic
            node$is_significant <- temp$is_significant
            
            # stop if split is not significant
            if (!node$is_significant) {
                return(NA)
            }
        }
        
        # If cross_validation_function is supplied, calculate cv_losses of both
        # subsegments and calculate inner cross-validated gain
        if (!is.null(cross_validation_function)) {
            temp_left <- cross_validation_function(x, start = node$start, 
                end = node$split_point, lambda = node$lambda, folds = node$root$folds)
            temp_right <- cross_validation_function(x, start = node$split_point, 
                end = node$end, lambda = node$lambda, folds = node$root$folds)
            
            node$cv_improvement <- node$cv_loss - temp_left$cv_loss - temp_right$cv_loss
            node$relative_cv_improvement <- node$cv_improvement/(node$end - 
                node$start)
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
    
    tree
    
}
