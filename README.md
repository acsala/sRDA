sRDA: the R package
========================

[![CRAN Version](https://www.r-pkg.org/badges/version/sRDA)](https://CRAN.R-project.org/package=sRDA) [![Downloads](https://cranlogs.r-pkg.org/badges/sRDA)](https://CRAN.R-project.org/package=sRDA) [![Total Downloads](https://cranlogs.r-pkg.org/badges/grand-total/sRDA?color=orange)](https://CRAN.R-project.org/package=sRDA)

R package providing a multivariate statistical tool for high dimensional (biomedical) data. sRDA is a a directional multivariate analysis to express the maximum variance in the predicted dataset by a linear combination of variables of the predictive dataset. Implemented in a partial least squares framework, for more details see [sRDA on CRAN](https://CRAN.R-project.org/package=sRDA) and [Csala et al. (2017)](https://doi.org/10.1093/bioinformatics/btx374).

Installation
------------

``` r
# from CRAN
install.packages("sRDA")
```

Usage
-----

Run sRDA with predefined penalization parameters.
``` r
# generate data with few highly correlated variahbles
dataXY <- generate_data(nr_LVs = 2,
                           n = 250,
                           nr_correlated_Xs = c(5,20),
                           nr_uncorrelated_Xs = 250,
                           mean_reg_weights_assoc_X =
                             c(0.9,0.5),
                           sd_reg_weights_assoc_X =
                             c(0.05, 0.05),
                           Xnoise_min = -0.3,
                           Xnoise_max = 0.3,
                           nr_correlated_Ys = c(10,15),
                           nr_uncorrelated_Ys = 350,
                           mean_reg_weights_assoc_Y =
                             c(0.9,0.6),
                           sd_reg_weights_assoc_Y =
                             c(0.05, 0.05),
                           Ynoise_min = -0.3,
                           Ynoise_max = 0.3)



# seperate predictor and predicted sets
X <- dataXY$X
Y <- dataXY$Y

# run sRDA
RDA.res <- sRDA(predictor = X, predicted = Y, nonzero = 5,
ridge_penalty = 1, penalization = "ust")

# check first 10 weights of X
RDA.res$ALPHA[1:10]
```

Run sRDA with cross-validation to determine best penalization parameters:

``` r
RDA.res <- sRDA(predictor = X, predicted = Y, nonzero = c(5,10,15),
ridge_penalty = c(0.1,1), penalization = "enet", cross_validate = TRUE,
parallel_CV = TRUE)

# check first 10 weights of X
RDA.res$ALPHA[1:10]

# check the Ridge parameter and the number of nonzeros included in the model
RDA.res$ridge_penalty
RDA.res$nr_nonzeros

# check how much time the cross validation did take
RDA.res$CV_results$stime
```

Obtain multiple latent variables (components):

``` r
RDA.res <- sRDA(predictor = X, predicted = Y, nonzero = c(5,10,15),
ridge_penalty = c(0.1,1), penalization = "enet", cross_validate = TRUE,
parallel_CV = TRUE, multiple_LV = TRUE, nr_LVs = 2, max_iterations = 5)

# check first 20 weights of X in first two component
RDA.res$ALPHA[[1]][1:20]
RDA.res$ALPHA[[2]][1:20]

# components are orthogonal to each other
t(RDA.res$XI[[1]]) %*% RDA.res$XI[[2]]
```

License
-------

This package as well as the source repositories are licensed under MIT.
