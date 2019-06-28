get_glasso_cross_validation_function <- function(control) {
    
    # get parameters
    n_folds <- control$cv_inner_n_folds
    search_lambda_inner <- control$cv_inner_search_lambda
    lambda_inner_step <- control$cv_inner_lambda_step
    lambda_inner <- control$cv_inner_lambda
    
    function(x, start, end, lambda, folds) {
        
        stopifnot(length(folds) == nrow(x))
        
        folds_current <- folds[(start + 1):end]
        
        # function that returns a matrix of dimension (n_folds, length(lambda))
        # with respective losses
        evaluate_fit <- function(lambda) {
            
            # prepare matrix to store values in later
            out <- array(NA, dim = c(length(unique(folds)), length(lambda)))
            
            # iterate over folds and lambdas
            for (l in 1:length(lambda)) {
                for (i in as.integer(unique(folds_current))) {
                  glasso_fit <- get_glasso_fit(x[(start + 1):end, , drop = F][folds_current != 
                    i, , drop = F], lambda[l], control)
                  out[i, l] <- sum(loglikelihood(x[(start + 1):end, glasso_fit$inds, 
                    drop = F][folds_current == i, , drop = F], glasso_fit$mu[glasso_fit$inds], 
                    glasso_fit$wi[glasso_fit$inds, glasso_fit$inds, drop = F]), 
                    na.rm = T)
                }
            }
            
            # return results
            out
        }
        
        
        if (search_lambda_inner) {
            
            # if search_lambda_inner is true, we search for an optimal lambda by
            # adjusting lambda0 in steps of size lambda_inner_step until a valley
            # is reached
            lambda <- c(1/lambda_inner_step, 1, lambda_inner_step) * lambda
            
            loss <- evaluate_fit(lambda)
            loss_sum <- apply(loss, 2, sum, na.rm = T)
            loss_sd <- apply(loss, 2, sd, na.rm = T)
            i <- which.min(loss_sum)
            
            if (i == 2) {
                # lambda already valley
                list(loss_array = loss[, i], lambda_opt = lambda[2], cv_loss = loss_sum[i])
            } else if (i == 1) {
                while (TRUE) {
                  # bad form, change this
                  lambda <- c(lambda[1]/lambda_inner_step, lambda)
                  loss_cur <- evaluate_fit(lambda[1])
                  loss_sum <- c(sum(loss_cur, na.rm = T), loss_sum)
                  loss_sd <- c(sd(loss_cur, na.rm = T), loss_sd)
                  loss <- cbind(loss_cur, loss)
                  if (loss_sum[1] >= loss_sum[2]) {
                    return(list(loss_array = loss[, 2], lambda_opt = lambda[2], 
                      cv_loss = loss_sum[2]))
                  }
                }
            } else {
                while (TRUE) {
                  k <- length(lambda)
                  lambda <- c(lambda, lambda[k] * lambda_inner_step)
                  loss_cur <- evaluate_fit(lambda[k + 1])
                  loss_sum <- c(loss_sum, sum(loss_cur, na.rm = T))
                  loss_sd <- c(loss_sd, sd(loss_cur, na.rm = T))
                  loss <- cbind(loss, loss_cur)
                  if (loss_sum[k + 1] >= loss_sum[k]) {
                    return(list(loss_array = loss[, k], lambda_opt = lambda[k], 
                      cv_loss = loss_sum[k]))
                  }
                }
            }
        } else {
            stopifnot(!is.null(lambda_inner))
            loss <- evaluate_fit(lambda_inner)
            loss_sum <- apply(loss, 2, sum)
            loss_sd <- apply(loss, 2, sd)
            i <- which.min(loss_sum)
            return(list(loss_array = loss[, i], lambda_opt = lambda_inner[i], 
                cv_loss = loss_sum[i]))
        }
    }
}
