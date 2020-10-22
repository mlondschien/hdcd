# High-dimensional change point detection
`hdcd` is a package for change point detection in (possibly high-dimensional) Gaussian Graphical Models (GGMs) with missing values. The package accompanies the paper _Change point detection for graphical models in presenceof missing values_ [1]. See  [here](https://arxiv.org/abs/1907.05409) for a preprint.

The package can be installed directly from github:
```R
# install.packages(devtools)
devtools::install_github("mlondschien/hdcd")
```

## Example
The following code snipped generates a `p = 100` dimensional time series of `n = 500` observations with three changepoints at locations `alpha = c(120, 240, 310)`, where in-segment distributions are Gaussian Graphical Models (GGM) of mean zero and where the precision matrices are generated as chain networks. Note that this is a truly high-dimensional setting, as the number of variables is higher than the length of the smallest segment, 70.

```R
set.seed(0)
model <- hdcd::create_model(n=500, p=100, c(120, 240, 310), hdcd::ChainNetwork)
x <- hdcd::simulate_from_model(model)
```

We then use the main algorithm of this package to find change points in the GGM structure. Note that we specify `method="glasso"` to use the [glasso](https://cran.r-project.org/web/packages/glasso/index.html) package [2] to estimate within segment precision matrices and that we use the `"section_search"` optimizer (also known as optimistic binary segmentation, OBS). This reduces computation time drastically. We use the default value of `delta = 0.1` for minimal relative segment length and let `hdcd` find a suitable initial regularization parameter `lamdba`. Note that running `hdcd` might take a few seconds due to the high-dimensionality of the dataset.
```R
tree <- hdcd::hdcd(x, method="glasso", optimizer="section_search")
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
The return object is a `binary_segmentation_tree`, which inherits from `data.tree::Node`. Each node in the tree corresponds to one segment, with notable attributes `start`, `split_point`, `end`, `gain` and `max_gain`. In our setting using the glasso method each node also has an attribute `cv_improvement`, which is calculated as the cross validated increase in likelihood when splitting the segment `(start, end]` at `split_point`. The attribute `lambda` is the regularization parameter optimal as determined by the cross validation procedure and used for the given segment for splitting. The `hdcd` algorithm stops splitting when `cv_improvement <= 0`. Note that `split_point`s with positive `cv_improvement` are given as 120, 189, 240, 310, which are exactly the true underlying change points up to one false positive. These can be extracted from the `tree` with the method `hdcd::get_change_points_from_tree`. 

```R
hdcd::get_change_points_from_tree(tree)
# [1] 120 189 240 310
```

## Missing values
The `hdcd` algorithm can also handle missing values. We delete 30% of entries of the matrix `x` completely at random and run our algorithm again. We use the Loh-Wainwrigth bias corrected approach to impute covariance matrices, which is the default value.

```R
x_deleted = hdcd::delete_values(x, 0.3, "mcar")
tree_deleted = hdcd(x_deleted, method="glasso", optimizer="section_search")
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
## Simulations and plots
The folders `simulations/main_simulation` and `simulation/plots` contain R-scripts and ReadMes to reproduce tables and plots of [1]. Raw data generated by the scripts is included in `data`. We don't include the raw data for the groundwater example, as the owners only permitted using it as an illustrating example.

## References
[1] M. Londschien, S. Kovacs and P. Buhlmann (2019), "Change point detection for graphical models in presenceof missing values", [arXiv:1907.05409](https://arxiv.org/abs/1907.05409)

[2] J. Friedman, T. Hastie and R. Tibshirani (2008), “Sparse inverse covariance estimation with thegraphical lasso", Biostatistics, 9, 432–441
