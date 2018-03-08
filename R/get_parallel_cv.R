
get_parallel_cv <- function(X,
                            Y,
                            lambdas,
                            non_zeros,
                            label,
                            penalization,
                            max_iterations,
                            tolerance){
#parallel n-fold cross validation
#
#input:     X             n*p matrix          - independent variables
#           Y             n*q matrix          - dependent variables
#           lambda        integer             - factor for Ridge penalty
#           labels        vector with n_subset unique elements with length of
#                         n - for crossvalidation
#
#output:    abs_cor       vector of integers -returns the summed abs correlation
#           stime         int/time          -return algorithms running time

    nr_subsets    <-    length(unique(label))


    #get working clusters initialized, otherwise parallel will use only one cpu*
    nr_cores      <-    detectCores()
    cl            <-    makeCluster(nr_cores)
    registerDoParallel(cl)
    getDoParWorkers()
    l <- lambdas

    #measure time
    stime         <-    system.time({

      #Get the alphas with parallel computing
      x <- foreach(l = lambdas, .combine = c, .export= c("enet",
                                                         "sRDA")) %dopar% {

        sub_abs_cor       <- c()
        sub_results       <- c()
        abs_cors          <- c()
        sub_iterations_m  <- c()
        iterations_m      <- c()

        for (nz in 1:length(non_zeros)){

          for (i in 1:nr_subsets){

            X.train   <- X[label!=i,]
            X.test    <- X[label==i,]


            Y.train   <- Y[label!=i,]
            Y.test    <- Y[label==i,]



            sub_results[[i]] <- sRDA(predictor = X.train,
                                     predicted = Y.train,
                                     ridge_penalty = l,
                                     nonzero = non_zeros[nz],
                                     max_iterations = max_iterations,
                                     penalization = penalization,
                                     tolerance = tolerance)


            XI.test = scale(X.test) %*% sub_results[[i]]$ALPH

            sub_iterations_m[[i]] <- sub_results[[i]]$nr_iterations
            sub_abs_cor[[i]] <- sum(abs(cor(XI.test,Y.test)), na.rm = T)

          }#end of subset for loop

          abs_cors        <- cbind(abs_cors, sub_abs_cor)
          iterations_m    <- cbind(iterations_m, sub_iterations_m)

        }#end of non_zeros loop

        list(abs_cors, iterations_m)

      }#end of lambda for loop
    })[3]#end of measure time

    #stop cluster
    stopCluster(cl)



    #mine out the iterations from vector x ************************************
    iterations_m          <- c()
    needed_x_for_abs_cors <- c()
    j                     <- 0

    for (i in (1:length(lambdas))){

      j                       <- j+1
      needed_x_for_abs_cors   <- cbind(needed_x_for_abs_cors, c(x[[j]]))

      j                       <- j+1
      iterations_m            <- cbind(iterations_m, c(x[[j]]))

    }

    #**************************************************************************

    x <- needed_x_for_abs_cors

    abs_cors = matrix(x, nrow=nr_subsets)

    #Figure out lambdas and non-zeros columns in results
    labels_non_zeros  <- rep(non_zeros, dim(abs_cors)[2]/length(non_zeros))

    labels_lambdas    <- rep(lambdas, each=dim(abs_cors)[2]/length(lambdas))

    all_abs_cors  <- rbind(labels_lambdas, labels_non_zeros, abs_cors)


    mean_abs_cors <- c()

    for (i in 1:length(labels_lambdas)){

      sub_result    <-  (c(labels_lambdas[i],
                           labels_non_zeros[i], mean(abs_cors[,i])))
      mean_abs_cors <- rbind(mean_abs_cors, sub_result)

    }
    rownames(mean_abs_cors)   <- NULL
    colnames(mean_abs_cors)   <- c("Ridge Penalty",
                                   "Number of Nonzeros",
                                   "mean Cor over CVs")


    #add number of iterations to mean_abs_cors
    line_in_parameters <- function(iterations_m,lambdas,non_zeros,nr_subsets){
      temp_iters    <- c()
      iterations_m2 <- c(iterations_m)

      for (l in (1:length(lambdas))){

        for (nz in (1:length(non_zeros))){

          temp_iters      <- rbind(temp_iters, iterations_m2[1:nr_subsets])
          iterations_m2    <- iterations_m2[(nr_subsets+1):length(iterations_m2)]

        }
      }
      temp_iters
    }

    mean_abs_cors <- cbind(mean_abs_cors,
                           line_in_parameters(iterations_m,
                                              lambdas,non_zeros,nr_subsets))


  #Return section**********************
  result        <-    list(abs_cors = abs_cors,
                           mean_abs_cors = mean_abs_cors,
                           stime = stime

  )

  result
}
