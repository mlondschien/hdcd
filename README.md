# High-dimensional change point detection
`hdcd` is a package for multivariate, or even **h**igh-**d**imensional **c**hange point **d**etection methods. It implements change point detection for (high-dimensional) Gaussian Graphical Models (GGMs) with possibly missing values as described in the paper [1]: _Londschien, Kovács and Bühlmann: "Change point detection for graphical models in the presence of missing values"_, see  [here](https://arxiv.org/abs/1907.05409) for a preprint. Additionally, the `hdcd` package also implements some ongoing work on multivariate nonparametric change point detection from [2].


## Installation of the hdcd package
The package can be installed directly from github:
```R
# install.packages(devtools)
library(devtools)
devtools::install_github("mlondschien/hdcd")
``` 

## Simulations and plots for the paper [1]
The folders `simulations/main_simulation` and `simulation/plots` contain R-scripts and READMEs to reproduce Table 1 with the main simulation and figures showing gain curves (Figure 2, 3 and 6). Figure 4 can be reproduced using the R-script in `simulations/histograms` and the `data_histograms` stored in the `data` folder. Other raw data generated by the scripts (e.g. for the main simulation of Table 1) are also included in `data` folder. We do not include the raw data for the groundwater example, as the owners only permitted using it as an illustrating example.

## Example: detecting changes in a GGM
The following code snippet generates a `p = 100` dimensional time series of `n = 500` observations with three change points at locations `alpha = c(120, 240, 310)`, where in-segment distributions arise from Gaussian Graphical Models (GGMs), i.e. multivariate Gaussians with mean zero and precision matrices that are generated as chain networks. Note that this is a truly high-dimensional setting, as the number of variables is higher than the length of the smallest segment, 70.

```R
set.seed(0)
model <- hdcd::create_model(n = 500, p = 100, c(120, 240, 310), hdcd::ChainNetwork)
x <- hdcd::simulate_from_model(model)
```

We then use the main algorithm from [1] to find change points in the GGM structure. Note that we specify `method = "glasso"` to use the graphical lasso ([3]) based loss and to estimate within segment precision matrices using the [glasso](https://cran.r-project.org/web/packages/glasso/index.html) package of [3]. We use Optimistic Binary Segmentation (OBS) here, as the `"section_search"` optimizer (i.e. naive Optimistic Search) was selected for finding the "best" split point candidate in each step and Binary Segmentation style algorithm is the default. This reduces computation time drastically. To perform the computationally more expensive full grid search, one can set the optimizer to `"line_search"`. This would then correspond to the traditional Binary Segmentation based change point detection. We use the default value of `delta = 0.1` for the minimal relative segment length and let `hdcd` find a suitable initial regularization parameter `lambda`. Note that running `hdcd` might take a few seconds due to the high-dimensionality of the simulated dataset.
```R
tree <- hdcd::hdcd(x, method = "glasso", optimizer = "section_search")
tree
#                    levelName split_point max_gain   cv_loss cv_improvement     lambda
# 1  (0 500]                           310 9.358010 14023.112     361.784038 0.06944723
# 2   ¦--(0 310]                       120 5.408330  8643.766     225.420988 0.09821322
# 3   ¦   ¦--(0 120]                    65 2.866239  3138.095    -108.436406 0.09821322
# 4   ¦   ¦   ¦--(0 65]                 NA       NA  1719.972             NA 0.19642644
# 5   ¦   ¦   °--(65 120]               NA       NA  1526.559             NA 0.19642644
# 6   ¦   °--(120 310]                 240 4.592640  5280.250     190.629298 0.13889447
# 7   ¦       ¦--(120 240]             189 2.251111  3160.787      60.582240 0.13889447
# 8   ¦       ¦   ¦--(120 189]          NA       NA  1763.735             NA 0.19642644
# 9   ¦       ¦   °--(189 240]          NA       NA  1336.470             NA 0.27778893
# 10  ¦       °--(240 310]              NA       NA  1928.833             NA 0.19642644
# 11  °--(310 500]                     397 3.162574  5017.562      -5.699871 0.06944723
# 12      ¦--(310 397]                  NA       NA  2311.062             NA 0.19642644
# 13      °--(397 500]                  NA       NA  2712.200             NA 0.13889447
```
The returned object is a `binary_segmentation_tree`, which inherits from `data.tree::Node`. Each node in the tree corresponds to one segment, with attributes `start`, `split_point`, `end`, `gain` and `max_gain`. In our setting using the glasso method each node also has an attribute `cv_improvement`, which is calculated as the cross validated increase in likelihood when splitting the segment `(start, end]` at `split_point`. The attribute `lambda` is the optimal regularization parameter (as determined by the cross validation procedure) and it used for the given segment for splitting. `hdcd` stops splitting when `cv_improvement <= 0`. Note that `split_point`s with positive `cv_improvement` are given as 310, 120, 240, 397, which are exactly the true underlying change points up to one false positive (at observation 397). These can be extracted from the saved `tree` with the method `hdcd::get_change_points_from_tree`. 

```R
hdcd::get_change_points_from_tree(tree)
# [1] 120 189 240 310
```

## Example: detecting changes in a GGM with missing values
`hdcd` can also handle missing values according to the methodology described in [1]. We delete 30% of entries of the matrix `x` completely at random and run our procedure again. We take the default imputation method, which is the Loh-Wainwrigth bias correction approach.

```R
x_deleted    <- hdcd::delete_values(x, 0.3, "mcar")
tree_deleted <- hdcd::hdcd(x_deleted, method = "glasso", optimizer = "section_search")
tree_deleted
> tree_deleted
#                    levelName split_point  max_gain  cv_loss cv_improvement     lambda
# 1  (0 500]                           310 5.8346049 9958.733       80.24636 0.08638148
# 2   ¦--(0 310]                       120 2.6951306 6171.951       41.34443 0.12216187
# 3   ¦   ¦--(0 120]                    64 1.4492266 2332.666       64.51703 0.12216187
# 4   ¦   ¦   ¦--(0 64]                 NA        NA 1163.530             NA 0.24432373
# 5   ¦   ¦   °--(64 120]               NA        NA 1104.620             NA 0.24432373
# 6   ¦   °--(120 310]                 240 1.8061935 3797.940       67.46934 0.17276297
# 7   ¦       ¦--(120 240]             174 0.9893486 2331.568       23.00253 0.17276297
# 8   ¦       ¦   ¦--(120 174]          NA        NA 1047.984             NA 0.24432373
# 9   ¦       ¦   °--(174 240]          NA        NA 1260.582             NA 0.24432373
# 10  ¦       °--(240 310]              NA        NA 1398.903             NA 0.17276297
# 11  °--(310 500]                     400 1.5274985 3706.536      -74.70367 0.12216187
# 12      ¦--(310 400]                  NA        NA 1772.515             NA 0.17276297
# 13      °--(400 500]                  NA        NA 2008.725             NA 0.17276297
 hdcd::get_change_points_from_tree(tree)
# [1] 120 189 240 310
```

## References
[1] M. Londschien, S. Kovács and P. Bühlmann (2019), "Change point detection for graphical models in the presence of missing values", [arXiv:1907.05409](https://arxiv.org/abs/1907.05409).

[2] S. Kovács, M. Londschien and P. Bühlmann (2020), "Random Forests and other nonparametric classifiers for multivariate change point detection", working paper.

[3] J. Friedman, T. Hastie and R. Tibshirani (2008), "Sparse inverse covariance estimation with the graphical lasso", Biostatistics, 9, 432–441.
