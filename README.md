# High-dimensional change point detection

`hdcd` is a package for change point detection in (possibly high-dimensional) Gaussian Graphical Models (GGMs) with missing values. The paper is meant to accompany the paper _Change point detection for graphical models in presenceof missing values_ [1]. See  [here](https://arxiv.org/abs/1907.05409) for a preprint.

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
tree <- hdcd(x, method="glasso", optimizer="section_search")
tree
#                    levelName split_point max_gain   cv_loss cv_improvement     lambda
# 1  (0 500]                           310 9.348791 14067.249      379.70836 0.06884284
# 2   ¦--(0 310]                       120 6.478121  8655.451      231.45848 0.06884284
# 3   ¦   ¦--(0 120]                    65 2.864449  3144.510     -127.41210 0.09735847
# 4   ¦   ¦   ¦--(0 65]                 NA       NA  1734.025             NA 0.19471694
# 5   ¦   ¦   °--(65 120]               NA       NA  1537.897             NA 0.27537134
# 6   ¦   °--(120 310]                 240 5.387394  5279.482      183.40675 0.09735847
# 7   ¦       ¦--(120 240]             180 2.281114  3156.330     -110.98839 0.13768567
# 8   ¦       ¦   ¦--(120 180]          NA       NA  1651.751             NA 0.19471694
# 9   ¦       ¦   °--(180 240]          NA       NA  1615.568             NA 0.27537134
# 10  ¦       °--(240 310]              NA       NA  1939.745             NA 0.19471694
# 11  °--(310 500]                     397 2.680704  5032.090       27.17972 0.09735847
# 12      ¦--(310 397]                  NA       NA  2312.857             NA 0.19471694
# 13      °--(397 500]                 450 2.256485  2692.053     -106.99206 0.13768567
# 14          ¦--(397 450]              NA       NA  1403.713             NA 0.27537134
# 15          °--(450 500]              NA       NA  1395.333             NA 0.27537134
```
The return object is a `binary_segmentation_tree`, which inherits from `data.tree::Node`. Each node in the tree corresponds to one segment, with notable attributes `start`, `split_point`, `end`, `gain` and `max_gain`. In our setting each node also has an attribute `cv_improvement`, which is calculated as the cross validated increase in likelihood when splitting the segment `(start, end]` at `split_point`. The attribute `lambda` is the regularization parameter optimal as determined by the cross validation procedure and used for the given segment for splitting. The `hdcd` algorithm stops splitting when `cv_improvement <= 0`. Note that `split_point`s with positive `cv_improvement` are given as 310, 120, 240, 397, which are exactly the true underlying change points up to one false positive. These can be extracted from the `tree` with the method `hdcd::get_change_points_from_tree`. 

```R
hdcd::get_change_points_from_tree(tree)
# [1] 120 240 310 397
```

## Missing values

The `hdcd` algorithm can also handle missing values. We delete 30% of entries of the matrix `x` completely at random and run our algorithm again. We use the Loh-Wainwrigth bias corrected approach to impute covariance matrices, which is the default value.

```R
x_deleted = hdcd::delete_values(x, 0.3, "mcar")
tree_deleted = hdcd(x_deleted, method="glasso", optimizer="section_search")
tree_deleted
#                   levelName split_point max_gain  cv_loss cv_improvement    lambda
# 1  (0 500]                           310 5.965655 9930.557     108.787376 0.1076932
# 2   ¦--(0 310]                       120 3.032775 6186.243      61.549010 0.1076932
# 3   ¦   ¦--(0 120]                    64 1.232258 2364.591     120.089053 0.1523011
# 4   ¦   ¦   ¦--(0 64]                 NA       NA 1165.593             NA 0.2153863
# 5   ¦   ¦   °--(64 120]               NA       NA 1078.908             NA 0.2153863
# 6   ¦   °--(120 310]                 240 2.088026 3760.104      29.236988 0.1523011
# 7   ¦       ¦--(120 240]             180 1.185872 2320.360     -42.192047 0.1523011
# 8   ¦       ¦   ¦--(120 180]          NA       NA 1206.920             NA 0.2153863
# 9   ¦       ¦   °--(180 240]          NA       NA 1155.632             NA 0.1523011
# 10  ¦       °--(240 310]              NA       NA 1410.507             NA 0.2153863
# 11  °--(310 500]                     392 1.747568 3635.527      -3.764511 0.1076932
# 12      ¦--(310 392]                  NA       NA 1590.383             NA 0.2153863
# 13      °--(392 500]                  NA       NA 2048.909             NA 0.1523011
hdcd::get_change_points_from_tree(tree)
# [1]  64 120 240 310
```