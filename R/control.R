#' Create an object of class hdcd_control to supply parameters to \link{hdcd}
#' 
#' @param glasso_cv_inner Should change points be selected according to inner cross validation in the glassocd procedure.
#' @param glasso_cv_inner_lambda Grid of values for lambda to be analyzed during inner cross validation. If NULL, values will
#' be searched for iteratively
#' @param glasso_cv_inner_randomize_folds Should inner folds be selected in a randomized fashion (instead of equispaced)
#' @param glasso_cv_inner_n_folds Number of inner folds
#' @param glasso_NA_method Method to be used for covariance estimation in the presence of missing values. One of 
#' \code{loh_wainwright_bias_correction}, \code{pairwise} or \code{average}. See accompanying paper for more details.
#' @param glasso_standardize Should values be be standardized  before they are fed into the glasso?
#' @param glasso_penalize_diagonal Should values on the diagonal of the precision matrix be penalized in the glasso?
#' @param classifier_oob Does the (custom) classifier return oob predictions?
#' @param section_search_stepsize Stepsize parameter for the section_search optimizer
#' @param random_forest_n_tree n_tree parameter for Random Forest algorithm ranger
#' @param random_forest_mtry mtry parameter for Random Forest algorithm ranger
#' @param random_forest_sample_fraction sample_fraction parameter for Random Forest algorithm ranger
#' @param kNN_k Number of nearest neighbors used for classification. Either a positive integer or a function
#' that reuturns a positive integer given the number of total observations.
#' @param wbs_n_segments Number of segments to be drawn for the WBS procedure
#' @param sbs_alpha Decay parameter for the SBS procedure
#' @param   permutation_test should a permutation test be done (approximated) for model selection (only relevant for kNNcd and RFcd)
#' @param permutation_test_pvalue pvalue threshold used in permutation test
#' @param permutation_test_n number of permutations (approximations) in permutation test
#' @export
hdcd_control <- function(glasso_cv_inner = TRUE,
                         glasso_cv_inner_lambda = NULL,
                         glasso_cv_inner_randomize_folds = FALSE,
                         glasso_cv_inner_n_folds = 5,
                         glasso_NA_method = 'loh_wainwright_bias_correction',
                         glasso_standardize = TRUE,
                         glasso_penalize_diagonal = FALSE,
                         #
                         classifier_oob = TRUE,
                         #
                         # section_search_min_points = 5,
                         section_search_stepsize = 0.1,
                         #
                         random_forest_n_tree = 600,
                         random_forest_mtry = NULL,
                         random_forest_sample_fraction = 1,
                         #
                         kNN_k = function(x) ceiling(sqrt(x)),
                         #
                         wbs_n_segments = 100,
                         sbs_alpha = 1/sqrt(2),
                         #
                         permutation_test = TRUE,
                         permutation_test_pvalue = 0.05,
                         permutation_test_n = 400
){
  structure(list(
    glasso_cv_inner = glasso_cv_inner,
    glasso_cv_inner_lambda = glasso_cv_inner_lambda,
    glasso_cv_inner_randomize_folds = glasso_cv_inner_randomize_folds,
    glasso_cv_inner_n_folds = glasso_cv_inner_n_folds,
    glasso_NA_method = glasso_NA_method,
    glasso_standardize = glasso_standardize,
    glasso_penalize_diagonal = glasso_penalize_diagonal,
    #
    classifier_oob = classifier_oob,
    #
    # section_search_min_points = 5,
    section_search_stepsize = section_search_stepsize,
    #
    random_forest_n_tree = random_forest_n_tree,
    random_forest_mtry = random_forest_mtry,
    random_forest_sample_fraction = random_forest_sample_fraction,
    #
    kNN_k = kNN_k,
    #
    wbs_n_segments = wbs_n_segments,
    sbs_alpha = sbs_alpha,
    #
    permutation_test = permutation_test,
    permutation_test_pvalue = permutation_test_pvalue,
    permutation_test_n = permutation_test_n
  ), class = "hdcd_control")
}
