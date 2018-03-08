# Cross validation for finding the optimal penalization parameters.

get_non_parallel_cv <- function(X,
                                Y,
                                lambdas,
                                non_zeros,
                                label,
                                penalization,
                                max_iterations,
                                tolerance){

  nr_subsets      <-    length(unique(label))
  abs_cors        <-    c()
  iterations_m    <-    c()
  kth_fold        <-    0

  #measure time
  stime <- system.time({
    for (l in 1:length(lambdas)){

      sub_abs_cor       <- c()
      sub_results       <- c()
      sub_iterations_m  <- c()

      for (nz in 1:length(non_zeros)){
        kth_fold <- kth_fold + 1

        for (i in 1:nr_subsets){

          X.train   <- X[label!=i,]
          X.test    <- X[label==i,]

          Y.train   <- Y[label!=i,]
          Y.test    <- Y[label==i,]


          sub_results[[i]] <- sRDA(predictor = X.train,
                                   predicted = Y.train,
                                   ridge_penalty = lambdas[l],
                                   nonzero = non_zeros[nz],
                                   penalization = penalization,
                                   max_iterations = max_iterations,
                                   tolerance = tolerance)

          XI.test = scale(X.test) %*% sub_results[[i]]$ALPH

          #devide with dim(Y.train)[2], exclude NA's from sum
          sub_abs_cor[[i]]            <- sum(abs(cor(XI.test,Y.test)),na.rm = T)

          sub_iterations_m[[i]]       <- sub_results[[i]]$nr_iterations


        }#end of subset for loop


        abs_cors        <- cbind(abs_cors, sub_abs_cor)
        colnames(abs_cors)[kth_fold]  <- paste(kth_fold, "-CVed abs_cors", sep="")
        iterations_m    <- cbind(iterations_m, sub_iterations_m)

      }#end of non_zeros loop

    }#end of lambda for loop
  })[3]#end of measure time

  #Figure out lambdas and non-zeros columns in results
  labels_non_zeros  <- rep(non_zeros, dim(abs_cors)[2]/length(non_zeros))

  labels_lambdas    <- rep(lambdas, each=dim(abs_cors)[2]/length(lambdas))

  all_abs_cors  <- rbind(labels_lambdas, labels_non_zeros, abs_cors)

  mean_abs_cors <- c()

  for (i in 1:length(labels_lambdas)){

    sub_result    <-  (c(labels_lambdas[i], labels_non_zeros[i], mean(abs_cors[,i])))
    mean_abs_cors <- rbind(mean_abs_cors, sub_result)

  }
  rownames(mean_abs_cors)   <- NULL
  colnames(mean_abs_cors)   <- c("Ridge Penalty",
                                 "Number of Nonzeros",
                                 "mean Cor over CVs")


  #Return section**********************
  result        <-    list(abs_cors = abs_cors,
                           mean_abs_cors = mean_abs_cors,
                           stime = stime,
                           iterations_m = iterations_m

  )


  result
}
