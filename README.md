sRDA: the R package
========================

[![CRAN Version](https://www.r-pkg.org/badges/version/sRDA)](https://CRAN.R-project.org/package=sRDA) [![Downloads](https://cranlogs.r-pkg.org/badges/sRDA)](https://CRAN.R-project.org/package=sRDA) [![Total Downloads](https://cranlogs.r-pkg.org/badges/grand-total/sRDA?color=orange)](https://CRAN.R-project.org/package=sRDA)

Multivariate statistical technique for high dimensional (biomedical) data analysis. sRDA is a a directional multivariate method to express the maximum variance in the predicted dataset by a linear combination of variables of the predictive dataset. Implemented in a partial least squares framework, for more details see [sRDA on CRAN](https://CRAN.R-project.org/package=sRDA) and [Csala et al. (2017)](https://doi.org/10.1093/bioinformatics/btx374).

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

Run sRDA with cross-validation to determine optimal penalization parameters:

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

Session Info
-------
``` r
sessionInfo()
R version 3.5.0 (2018-04-23)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Debian GNU/Linux buster/sid

Matrix products: default
BLAS: /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.8.0
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.8.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
[1] sRDA_1.0.0        mvtnorm_1.1-3     elasticnet_1.3    lars_1.3         
[5] doParallel_1.0.17 iterators_1.0.14  foreach_1.5.2     Matrix_1.2-14    

loaded via a namespace (and not attached):
[1] compiler_3.5.0   tools_3.5.0      codetools_0.2-15 grid_3.5.0      
[5] lattice_0.20-35
```

License
-------

This package as well as the source repositories are licensed under MIT.
