#        *         *         *         *         *         *         *         #
# Metadata                                                                  ####
#
# Title:    sparese RDA
# meta:     sparese RDA
#
# code by:  Attila Csala
# place:    Amsterdam
# date:     2016-11
#
# setwd("C:/amc.intra/users/A/acsala/home/sRDA/")
# devtools::load_all()
# devtools::use_package_doc()
# devtools::document()
# devtools::check()
# devtools::use_vignette("sRDA")
# devtools::build_vignettes()
# devtools::build()
# devtools::install()
# install.packages("C:/amc.intra/users/A/acsala/home/sRDA_0.1.tar.gz")
#        *         *         *         *         *         *         *         #

#' @title
#' Sparse Canonical Correlation analysis
#'
#' @description
#' Sparse Canonical Correlation analysis for high dimensional (biomedical) data.
#' The function takes two datasets and returns a linear combination of
#' maximally correlated canonical correlate pairs.  Elastic net penalization
#' (with its variants, UST, Ridge and Lasso penalization) is implemented for
#' sparsity and smoothnesswith a  built in cross validation procedure to obtain
#' the optimal penalization parameters. It is possible to obtain multiple
#' canonical variate pairs that are orthogonal to each other.
#'
#'
#' @param predictor
#' The n*p matrix of the predictor data set
#' @param predicted
#' The n*q matrix of the predicted data set
#' @param penalization
#' The penalization method applied during the analysis (none, enet or ust)
#' @param ridge_penalty
#' The ridge penalty parameter of the predictor set's latent variable used
#' for enet or ust (an integer if cross_validate = FALE, a list otherwise)
#' @param nonzero
#' The number of non-zero weights of the predictor set's latent variable
#' (an integer if cross_validate = FALE, a list otherwise)
#' @param max_iterations
#' The maximum number of iterations of the algorithm
#' @param tolerance
#' Convergence criteria
#' @param cross_validate
#' K-fold cross validation to find best optimal penalty parameters
#' (TRUE or FALSE)
#' @param parallel_CV
#' Run the cross validation parallel (TRUE or FALSE)
#' @param nr_subsets
#' Number of subsets for k-fold cross validation
#' @param multiple_LV
#' Obtain multiple latent variable pairs (TRUE or FALSE)
#' @param nr_LVs
#' Number of latent variables to be obtained
#'
#' @return An object of class \code{"sRDA"}.
#' @return \item{XI}{
#' Predictor set's latent variable scores}
#' @return \item{ETA}{
#' Predictive set's latent variable scores}
#' @return \item{ALPHA}{
#' Weights of the predictor set's latent variable}
#' @return \item{BETA}{
#' Weights of the predicted set's latent variable}
#' @return \item{nr_iterations}{
#' Number of iterations ran before convergence (or max number of iterations)}
#' @return \item{SOLVE_XIXI}{
#' Inverse of the predictor set's latent variable variance matrix}
#' @return \item{iterations_crts}{
#' The convergence criterion value (a small positive tolerance)}
#' @return \item{sum_absolute_betas}{
#' Sum of the absolute values of beta weights}
#' @return \item{ridge_penalty}{
#' The ridge penalty parameter used for the model}
#' @return \item{nr_nonzeros}{
#' The number of nonzero alpha weights in the model}
#' @return \item{nr_latent_variables}{
#' The number of latient variable pairs in the model}
#' @return \item{CV_results}{
#' The detailed results of cross validations (if cross_validate is TRUE)}
#'
#' @author Attila Csala
#' @keywords Redundancy RDA
#' @import elasticnet
#' @import mvtnorm
#' @import parallel
#' @import doParallel
#' @import Matrix
#' @export
#' @examples
#'
#' # generate data with few highly correlated variahbles
#' dataXY <- generate_data(nr_LVs = 2,
#'                            n = 250,
#'                            nr_correlated_Xs = c(5,20),
#'                            nr_uncorrelated_Xs = 250,
#'                            mean_reg_weights_assoc_X =
#'                              c(0.9,0.5),
#'                            sd_reg_weights_assoc_X =
#'                              c(0.05, 0.05),
#'                            Xnoise_min = -0.3,
#'                            Xnoise_max = 0.3,
#'                            nr_correlated_Ys = c(10,15),
#'                            nr_uncorrelated_Ys = 350,
#'                            mean_reg_weights_assoc_Y =
#'                              c(0.9,0.6),
#'                            sd_reg_weights_assoc_Y =
#'                              c(0.05, 0.05),
#'                            Ynoise_min = -0.3,
#'                            Ynoise_max = 0.3)
#'
#' # seperate predictor and predicted sets
#' X <- dataXY$X
#' Y <- dataXY$Y
#'
#' # run sRDA
#' CCA.res <- sCCA(predictor = X, predicted = Y, nonzero = 5,
#' ridge_penalty = 1, penalization = "ust")
#'
#'
#' # check first 10 weights of X
#' CCA.res$ALPHA[1:10]
#'
#'\dontrun{
#' # run sRDA with cross-validation to determine best penalization parameters
#' CCA.res <- sCCA(predictor = X, predicted = Y, nonzero = c(5,10,15),
#' ridge_penalty = c(0.1,1), penalization = "enet", cross_validate = TRUE,
#' parallel_CV = TRUE)
#'
#' # check first 10 weights of X
#' CCA.res$ALPHA[1:10]
#' CCA.res$ridge_penalty
#' CCA.res$nr_nonzeros
#'
#' # obtain multiple latent variables
#' CCA.res <- sCCA(predictor = X, predicted = Y, nonzero = c(5,10,15),
#' ridge_penalty = c(0.1,1), penalization = "enet", cross_validate = TRUE,
#' parallel_CV = TRUE, multiple_LV = TRUE, nr_LVs = 2, max_iterations = 5)
#'
#' # check first 10 weights of X in first two component
#' CCA.res$ALPHA[[1]][1:10]
#' CCA.res$ALPHA[[2]][1:10]
#'
#' # latent variables are orthogonal to each other
#' t(CCA.res$XI[[1]]) %*% CCA.res$XI[[2]]
#'
#'}


#******************************************************************************#
#                               Peanalized CCA                              ####
#******************************************************************************#
#+-------+                                                            +-------+
#|       |                                                            |       |
#|       |                                                            |       |
#|       |                                                            |       |
#|       |   ALPHA       +------+          +------+       BETA        |       |
#|       +--------------->      |          |      <-------------------+       |
#|       |               |      |          |      |                   |       |
#|       |               |      |          |      |                   |       |
#|  X    +--------------->  XI  |          | ETA  <-------------------+   Y   |
#|       |               |      <---------->      |                   |       |
#|       |               |      |          |      |                   |       |
#|       +--------------->      |          |      <-------------------+       |
#|       |               +------+          +------+                   |       |
#|       |                                                            |       |
#|       |                                                            |       |
#|       |                                                            |       |
#+-------+                                                            +-------+




sCCA <- function(predictor,
                 predicted,
                 penalization = "enet",
                 ridge_penalty = 1,
                 nonzero = 1,
                 #stop criteruim in a form of maximum iterations
                 max_iterations = 100,
                 #stop criterium regarding to CRT
                 tolerance = 1*10^-20,
                 cross_validate = FALSE,
                 parallel_CV = TRUE,
                 nr_subsets = 10,
                 multiple_LV = FALSE,
                 nr_LVs = 1
){

  #****************************************************************************#
  #0. ancillary functions
  #*********************
  #check for dependencies
  if (!requireNamespace("elasticnet", quietly = TRUE)) {
    stop("The elasticnet package is needed for this function to work.
         Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("stats", quietly = TRUE)) {
    stop("The stats package is needed for this function to work.
         Please install it.",
         call. = FALSE)
  }

  #if ust penalization, change ridge penalties to Inf
  if (penalization == "ust"){
    ridge_penalty <- rep(Inf, length(ridge_penalty))
  }

  #check for dependencies for parallel computing
  if (parallel_CV == TRUE){

    #check for dependencies
    if (!requireNamespace("parallel", quietly = TRUE)) {
      stop("The parallel package is needed for this function to work.
           Please install it.",
           call. = FALSE)
    }
    if (!requireNamespace("foreach", quietly = TRUE)) {
      stop("The foreach package is needed for this function to work.
           Please install it.",
           call. = FALSE)
    }

    #check for cores
    nr_cores      <-    detectCores()
    if(nr_cores){
      cat("Cross valdiation function on multiple cores is called,
          number of detected cores on this work station:", nr_cores, "\n")
    }else{
      stop("R wasn't able to identify multiple cores on this work station.
           Please make sure your work station is compatible with
           parallel computing.",
           call. = FALSE)
    }

  }

  #Check input values for the function
  if( dim(predictor)[1] != dim(predicted)[1]){
    stop("The sample size for the predictor and predicted data sets is
         different. Please use data sets with equal sample size.",
         call. = FALSE)
  }

  #Check for cross validation *************************************************#

  if (multiple_LV){

    result <- get_multiple_LVs(X = predictor,
                               Y = predicted,
                               penalization = penalization,
                               lambda = ridge_penalty,
                               nonzero = nonzero,
                               nr_latent = nr_LVs,
                               stop_criterium = tolerance,
                               max_iterations = max_iterations,
                               cross_validate = cross_validate)

    return(result)


  }

  CV_results <- "Cross validation function is not called"

  if (cross_validate){

    if (penalization == "none") {
      stop("You do not need cross validation for a non penalized model
           (pelase change the penalization or cross_validate variable)",
           call. = FALSE)
    }

    CV_results <-
      get_cross_validated_penalty_parameters(predictor = predictor,
                                             predicted = predicted,
                                             penalization = penalization,
                                             ridge_penalty = ridge_penalty,
                                             nonzero = nonzero,
                                             nr_subsets = nr_subsets,
                                             max_iterations = max_iterations,
                                             tolerance = tolerance,
                                             parallel_CV = parallel_CV
      )

    ridge_penalty   <- CV_results$best_ridge
    nonzero         <- CV_results$best_nonzero

  } else {
    # if no cross validation and more then 1 ridge/non zero penalty, give error
    if (length(ridge_penalty)>1 || length(nonzero)>1) {
      stop("Multiple ridge_penalty and/or nonzero parameter where only one is
           needed (pelase set cross_validate = TRUE if you would like to run
           cross validation)",
           call. = FALSE)
    }

  }


  #Algorithm for peanalization*************************************************#

  #****************************************************************************#
  #1. Preperation of the data
  #*************************#

  Y.mat <- as.matrix(predicted)
  #scale (SD = 1) and center (mean = 0) Y
  Yc <- scale(Y.mat)
  Yc <- Yc[,!colSums(!is.finite(Yc))]

  X.mat <- as.matrix(predictor)
  #scale (SD = 1) and center (mean = 0) X
  Xcr <- scale(X.mat)
  Xcr <- Xcr[,!colSums(!is.finite(Xcr))]

  #variable length of X and Y
  p <- ncol(Xcr)
  q <- ncol(Yc)

  #   For initalization, let:
  #
  #         ETA^(0) = y_1 + y_2 + .... + y_q = Y*BETA_hat^(0) ( = Y_i)
  #           were BETA^(0) = [1,1....,1]
  # an y-weight column vector of q elements
  BETA = matrix( rep(1, q), nrow = q, ncol = 1, byrow = TRUE)
  # an x-weight column vector of p elements
  ALPHA = matrix( rep(1, p), nrow = p, ncol = 1, byrow = TRUE)


  #****************************************************************************#
  #2. Iterative loop until convergence
  #*************************#
  CRTs          <- c()
  sum_abs_Betas <- c()
  Nr_iterations = 0         #iteration counter

  WeContinnue = TRUE        #T/F value for stop criterium based on CRT values
  #in last two iteration
  CRT = 1                   #convergance measure between alpha and beta

  while(CRT > tolerance && WeContinnue && Nr_iterations < max_iterations) {

    ETA = Yc %*% BETA
    ETA = scale(ETA)


    XI = Xcr %*% ALPHA
    XI = scale(XI)


    #**************************************************************************#
    #Chose between generic RDA, SRDA with ENET or SRDA wiht UST****************#

    # For the value of ahat^(0) (ALPHA) and hence XIhat^(0),
    # regress ETAhat^(0) on X jointly to get ALPHA,
    # ALPH_0 = (X'X)^-1 (X' ETA)

    ALPH_0 <- switch(penalization,
                     #lm without penalization
                     #ALPH_0 = as.matrix(lm(ETA ~ 0+Xcr)$coefficients)
                     "none" = {
                       solve(t(Xcr) %*% Xcr) %*% t(Xcr) %*% ETA
                     },

                     #CALCULATE WITH ELASTIC NET
                     #ALPH_0 = as.matrix(calculateVectorEnet(Xcr,ETA,1,30))
                     "enet" = {
                       as.matrix(get_enet(Xcr,ETA,ridge_penalty,nonzero))
                     },

                     #CALCULATE WITH UST
                     #ALPH_0 = as.matrix(calculateVectorEnet(Xcr,ETA,1,30))
                     "ust" = {
                       as.matrix(get_ust(Xcr,ETA,nonzero))
                     },
                     {
                       stop("Please choose a valid penalization method
                            (i.e. 'ust', 'enet' or 'none')",
                            call. = FALSE)

                     }

                       )

    #***************************************************************************
    #         compute the value of XI^(1):
    #         XIhat^(1) = SUM_running to p where t=1 ( ahat_t^(0) *x_t )

    XI = Xcr %*% ALPH_0

    #         and normalize XIhat^(1) such that
    #           t(XIhat^(1))*XIhat^(1) = t(ahat^(0)) t(X)*X*ahat^(0) = 1
    #           t(XI)%*%XI = t(ALPH_0) %*% t(Xcr) %*% Xcr %*% ALPH_0 = 1
    #           that is its variance is 1
    XI = scale(XI)


    #         For the value BETAhat^(1) and hence ETAhat^(1),
    #           regress y1,y2 ... yq separately on XIhat^(1),
    #
    #           y1 = BETA_1 * XIhat^(1) + Epsilon_1
    #           .
    #           yq = BETA_q * XIhat^(1) + Epsilon_q

    BETA_0 <- switch(penalization,
                     #lm without penalization
                     "none" = {
                       solve(t(Yc) %*% Yc) %*% t(Yc) %*% XI
                     },

                     #CALCULATE WITH ELASTIC NET
                     "enet" = {
                       as.matrix(get_enet(Yc,XI,ridge_penalty,nonzero))
                     },

                     #CALCULATE WITH UST
                     "ust" = {
                       as.matrix(get_ust(Yc,XI,nonzero))
                     },
                     {
                       stop("Please choose a valid penalization method
                            (i.e. 'ust', 'enet' or 'none')",
                            call. = FALSE)

                     }

                       )

    #BETA
    #         compute the vaule of ETAhat^(1),
    #
    #           ETAhat^(1) = SUM_running to q where k=1 ( BETAhat_k^(1) * y_k)
    ETA = Yc %*% BETA_0

    ETA = scale(ETA)

    #check if ETA needed to be scaled here <- nope
    #ETA = scale(ETA)


    #Calculate convergence of Alphas and Betas
    CRT = sum((ALPHA - ALPH_0)^2, (BETA - BETA_0)^2);


    ALPHA=            ALPH_0
    BETA =            BETA_0

    #Check if last two iterations CR converges*********************************#
    Nr_iterations                   =   Nr_iterations + 1
    CRTs[[Nr_iterations]]           =   CRT
    sum_abs_Betas[[Nr_iterations]]  =   sum(abs(BETA))

    if (Nr_iterations>1){

      stop_condition <- abs(CRTs[[Nr_iterations]] - CRTs[[Nr_iterations-1]])
      stop_criterium <- 1 * 10^-6

      if (stop_condition < stop_criterium){
        WeContinnue <- FALSE
      }

    }#END Check if last two iterations CR converges*************************#

  }# End of main loop

  #REDN = REDUNDANCY INDEX: redundancy of y variables
  #  REDN = (sum(BETA^2))/4
  Red_index = (sum(BETA^2))/4

  #************************

  #  t(BETA) %*% t(XYMA) %*% solve(XXMA) %*% XYMA %*% BETA

  #  PATH


  #***************************

  #use this later for second latent variable
  SOLVE_XIXI = solve(t(XI)%*%XI) %*% t(XI)


  result <- list(XI = XI,
                 ETA = ETA,
                 ALPHA = ALPHA,
                 BETA= BETA,
                 nr_iterations = Nr_iterations,
                 inverse_of_XIXI = SOLVE_XIXI,
                 iterations_crts = CRTs,
                 sum_absolute_betas = sum_abs_Betas,
                 ridge_penalty = ridge_penalty,
                 nr_nonzeros = nonzero,
                 nr_latent_variables = nr_LVs,
                 CV_results = CV_results
  )

  class(result) = "sRDA"
  return(result)

                     }
