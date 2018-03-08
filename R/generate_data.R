#        *         *         *         *         *         *         *         #
# Metadata                                                                  ####
#
# Title:    sparese RDA
# meta:     sparese RDA
#
# code by:  Koos Zwinderman, Attila Csala
# place:    Amsterdam
# date:     2016-11
#        *         *         *         *         *         *         *         #

#' @title Generate data sets for sparse multivariate analysis
#'
#' @description
#' Generate two data sets with highly correlated and noise variables modeled
#' in a multiple latent variable structure. The latent variables are orthogonal
#' to each other thus capture a different portion of association between the
#' involved data sets. Thu function generates data that can be used to verify
#' sRDA's ability of finding the highly correlated variables accross multiple
#' latent variables.
#'
#' @param nr_LVs The number of latent variables between the predicitve
#' and predicted data sets. The latent variables model the association between
#' data sets.
#' @param n The number of observations (rows) in the data sets.
#' @param nr_correlated_Xs Number of variables of the
#' predictive data set that are associated with the latent variables.
#' @param nr_uncorrelated_Xs Number of variables of the predictive data
#' set that is not associated with the latent variables.
#' @param mean_reg_weights_assoc_X Mean of the
#' regression weights of the predictive varaibles that are associated with the
#' latent variables.
#' @param sd_reg_weights_assoc_X Standard deviation
#' of the regression weights of the predictive varaibles that are associated
#' with the latent variables.
#' @param Xnoise_min The lower bound of the unifrom distribution that is used to
#' sample the values for the regression weights of the predictive varaibles that
#' are not associated with the latent variables.
#' @param Xnoise_max The upper bound of the unifrom distribution that is used to
#' sample the values for the regression weights of the predictive varaibles that
#' are not associated with the latent variables.
#' @param nr_correlated_Ys Number of variables of the
#' predictive data set that are associated with the latent variables.
#' @param nr_uncorrelated_Ys Number of variables of the predicted data
#' set that is not associated with the latent variables.
#' @param mean_reg_weights_assoc_Y Mean of the
#' regression weights of the predicted varaibles that are associated with the
#' latent variables.
#' @param sd_reg_weights_assoc_Y Standard deviation
#' of the regression weights of the predicted varaibles that are associated
#' with the latent variables.
#' @param Ynoise_min The lower bound of the unifrom distribution that is used to
#' sample the values for the regression weights of the predicted varaibles that
#' are not associated with the latent variables.
#' @param Ynoise_max The upper bound of the unifrom distribution that is used to
#' sample the values for the regression weights of the prediced varaibles that
#' are not associated with the latent variables.
#' @keywords generate data
#'
#' @export
#'
#' @examples
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
#' dim(X);dim(Y)
#'
#' @import mvtnorm

#******************************************************************************#
#                        Generate Data a'la Koos                       #########
#******************************************************************************#

generate_data <- function(nr_LVs = 1,
                          n = 50,
                          nr_correlated_Xs = c(5),
                          nr_uncorrelated_Xs = 250,
                          mean_reg_weights_assoc_X = c(0.7),
                          sd_reg_weights_assoc_X = c(0.05),
                          Xnoise_min = -0.3,
                          Xnoise_max = 0.3,
                          nr_correlated_Ys = c(5),
                          nr_uncorrelated_Ys = 350,
                          mean_reg_weights_assoc_Y =  c(0.7),
                          sd_reg_weights_assoc_Y =    c(0.05),
                          Ynoise_min = -0.3,
                          Ynoise_max = 0.3){

  number_of_ksi                       <- nr_LVs
  number_of_patients                  <- n

  number_of_Xs_associated_with_ksis   <- nr_correlated_Xs
  number_of_not_associated_Xs         <- nr_uncorrelated_Xs

  mean_of_the_regression_weights_of_the_associated_Xs <-mean_reg_weights_assoc_X
  sd_of_the_regression_weights_of_the_associated_Xs   <-sd_reg_weights_assoc_X

  number_of_Ys_associated_with_ksis   <- nr_correlated_Ys
  number_of_not_associated_Ys         <- nr_uncorrelated_Ys

  mean_of_the_regression_weights_of_the_associated_Ys <-mean_reg_weights_assoc_Y
  sd_of_the_regression_weights_of_the_associated_Ys   <-sd_reg_weights_assoc_Y


  #*Generate latent vector ski from lultivariATE NORMAL distribution with
  #    mean 0, and covariance matrix = 1
  ksi=rmvnorm(number_of_patients,mean=rep(0,number_of_ksi),
              sigma=diag(number_of_ksi))

  # make X data
  #*make empty matrix with nrow and ncol
  x=matrix(NA,
           nrow=number_of_patients,
           ncol=(sum(number_of_Xs_associated_with_ksis)
                 +number_of_not_associated_Xs))
  columncount=0

  #i loops trhough number of ksi's
  for (i in 1:length(number_of_Xs_associated_with_ksis)) {

    #j loops throguh number of associated X's
    #8 since its starts from 1:, it will make the
    #  first j elements associated with ksi
    for (j in 1:number_of_Xs_associated_with_ksis[i]) {

      #calculates a single regression weight for the associated X
      regressionweight =
          rnorm(1,
                mean_of_the_regression_weights_of_the_associated_Xs[i],
                sd_of_the_regression_weights_of_the_associated_Xs[i])

      columncount=columncount+1

      #sd cannot be lower than 0
      x[,columncount] = rnorm(number_of_patients,mean=regressionweight*ksi[,i],
                              sd=sqrt(max(0,(1-regressionweight^2))+0.001))

    }

  }

  #fill in number of not associated columns
  for (j in 1:number_of_not_associated_Xs) {

    columncount=columncount+1
    x[,columncount] = runif(number_of_patients,min=Xnoise_min,max=Xnoise_max)
  }



  # make Y data
  y=matrix(NA,
           nrow=number_of_patients,ncol=(sum(number_of_Ys_associated_with_ksis)
                                         +number_of_not_associated_Ys))

  columncount=0

  for (i in 1:length(number_of_Ys_associated_with_ksis)) {

    for (j in 1:number_of_Ys_associated_with_ksis[i]) {

      regressionweight =
          rnorm(1,
                mean_of_the_regression_weights_of_the_associated_Ys[i],
                sd_of_the_regression_weights_of_the_associated_Ys[i])
      columncount=columncount+1
      y[,columncount] =
          rnorm(number_of_patients,mean=regressionweight*ksi[,i],
                sd=max(0,(1-regressionweight^2))+0.001)

    }
  }


  for (j in 1:number_of_not_associated_Ys) {
    columncount=columncount+1
    y[,columncount] = runif(number_of_patients,min=Ynoise_min,max=Ynoise_max)
  }

  X <- x
  Y <- y
  data_info <- data.frame(number_of_ksi,
                          number_of_patients,
                          number_of_Xs_associated_with_ksis,
                          number_of_not_associated_Xs,
                          mean_of_the_regression_weights_of_the_associated_Xs,
                          sd_of_the_regression_weights_of_the_associated_Xs,
                          Xnoise_min, Xnoise_max,
                          number_of_Ys_associated_with_ksis,
                          number_of_not_associated_Ys,
                          mean_of_the_regression_weights_of_the_associated_Ys,
                          sd_of_the_regression_weights_of_the_associated_Ys,
                          Ynoise_min, Ynoise_max)

  result <- list(
    X,Y,data_info
  )

  names(result) <- c(
    "X","Y","data_info"
  )

  result

}
