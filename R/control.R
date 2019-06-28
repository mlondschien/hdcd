#' Create an object of class hdcd_control to supply parameters to \link{hdcd}
#' 
#' @export
hdcd_control <- function(cv_inner = FALSE, cv_inner_lambda = NULL, cv_inner_randomize_folds = FALSE, 
    cv_inner_n_lambda = 10, cv_inner_n_folds = 10, cv_inner_nrep_em = 5, 
    cv_inner_search_lambda = F, cv_inner_min_grid_ratio = 0.01, cv_inner_lambda_step = (0.1)^(1/10), 
    cv_inner_stop_early = FALSE, glasso_standardize = TRUE, glasso_penalize_diagonal = FALSE, 
    glasso_threshold = 1e-04, glasso_NA_method = "loh_wainwright_bias_correction", 
    nodewise_regression_node = 1, elastic_net_alpha = 1, elastic_net_family = "gaussian", 
    verbose = F, max_depth = Inf, section_search_min_points = 5, section_search_stepsize = 0.5, 
    section_search_k_sigma = 0, section_search_tolerance = 0, segment_loss_min_points = 2, 
    n_obs = NULL) {
    structure(list(cv_inner = cv_inner, cv_inner_randomize_folds = cv_inner_randomize_folds, 
        cv_inner_lambda = cv_inner_lambda, cv_inner_n_lambda = cv_inner_n_lambda, 
        cv_inner_n_folds = cv_inner_n_folds, cv_inner_nrep_em = cv_inner_nrep_em, 
        cv_inner_search_lambda = cv_inner_search_lambda, cv_inner_min_grid_ratio = cv_inner_min_grid_ratio, 
        cv_inner_lambda_step = cv_inner_lambda_step, cv_inner_stop_early = cv_inner_stop_early, 
        glasso_standardize = glasso_standardize, glasso_penalize_diagonal = glasso_penalize_diagonal, 
        glasso_threshold = glasso_threshold, glasso_NA_method = glasso_NA_method, 
        nodewise_regression_node = nodewise_regression_node, elastic_net_alpha = elastic_net_alpha, 
        elastic_net_family = elastic_net_family, verbose = verbose, max_depth = max_depth, 
        section_search_min_points = section_search_min_points, section_search_stepsize = section_search_stepsize, 
        section_search_k_sigma = section_search_k_sigma, section_search_tolerance = section_search_tolerance, 
        segment_loss_min_points = segment_loss_min_points, n_obs = n_obs), 
        class = "hdcd_control")
}
