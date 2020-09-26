test_that("glasso line_search", {
  set.seed(0)
  x <- hdcd::simulate_from_model(hdcd::create_model(100, 5, c(20, 60), hdcd::RandomNetwork))
  tree <- hdcd::hdcd(x, method="glasso", optimizer="line_search")
  change_points = hdcd::get_change_points_from_tree(tree)
  testthat::expect_equal(change_points, c(57))
})

test_that("glasso line_search", {
  set.seed(0)
  x <- hdcd::simulate_from_model(hdcd::create_model(100, 5, c(20, 60), hdcd::RandomNetwork))
  tree <- hdcd::hdcd(x, method="glasso", optimizer="section_search")
  change_points = hdcd::get_change_points_from_tree(tree)
  testthat::expect_equal(change_points, c(57))
})
