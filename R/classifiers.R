random_forest <- function(x_train, y_train, control, x_test = NULL) {
  stopifnot(nrow(x_train) == length(y_train))
  
  # fit a random forest
  fit <-
    ranger::ranger(
      y ~ .,
      data = data.frame(x = x_train, y = as.factor(y_train)),
      num.trees = control$random_forest_n_tree,
      mtry = control$random_forest_mtry,
      sample.fraction = control$random_forest_sample_fraction,
      keep.inbag = FALSE,
      probability = TRUE,
      classification = TRUE,
      seed = 0
    )
  
  if (is.null(x_test)) {
    # in case no additional test data was supplied, use fitted OOB values
    list(train = fit$predictions)
  } else {
    list(
      test = stats::predict(fit, data = data.frame(x = x_test))$predictions,
      train = fit$predictions
    )
  }
}