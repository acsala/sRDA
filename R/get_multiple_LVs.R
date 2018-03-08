get_multiple_LVs   <-  function(X,
                                Y,
                                penalization,
                                lambda,
                                nonzero,
                                nr_latent=1,
                                stop_criterium = 1 * 10^-5,
                                max_iterations,
                                cross_validate) {


  Res_X         <-      X
  alphas        <-      c()
  betas         <-      c()
  xis           <-      c()
  etas          <-      c()
  iterations    <-      c()
  corr_v        <-      c()
  s_cond_v      <-      c()
  red_indexs    <-      c()
  ridge_penaltys<-      c()
  nr_nonzeros   <-      c()
  CV_results    <-      c()

  iterations_crts <-    c()

  sum_of_sq_betas     <- c()
  sum_of_sq_alphas    <- c()

  i             <- 1
  WeCarryOn     <- TRUE
  cat("Multiple latent variables scenario,
      number of latent variables calculated:",nr_latent, "\n")

  while ( !(i > nr_latent) && WeCarryOn ){

    results       <-    sRDA(predictor = Res_X,
                             predicted = Y,
                             penalization = penalization,
                             ridge_penalty = lambda,
                             nonzero = nonzero,
                             tolerance = stop_criterium,
                             # cross validate for every latent variables
                             cross_validate = cross_validate,
                             multiple_LV = FALSE,
                             max_iterations = max_iterations)

    alphas[[i]]           <-  results$ALPHA
    betas[[i]]            <-  results$BETA
    xis[[i]]              <-  results$XI
    etas[[i]]             <-  results$ETA
    iterations[[i]]       <-  results$nr_iterations
    red_indexs[[i]]       <-  results$redundancy_index
    iterations_crts[[i]]  <-  results$iterations_crts
    ridge_penaltys[[i]]   <-  results$ridge_penalty
    nr_nonzeros[[i]]      <-  results$nr_nonzeros

    if(cross_validate){
      CV_results[[i]]     <- results$CV_results
    }

    reg_coeff     <-    results$inverse_of_XIXI %*% as.matrix(Res_X)

    # calculate the residuals
    calcres = function(Xcol)
      Xcol - results$inverse_of_XIXI %*% Xcol %*% t(xis[[i]])

    Res_X = apply(Res_X, 2,calcres)

    sum_of_sq_betas[[i]]    <- sum(betas[[i]]^2)
    sum_of_sq_alphas[[i]]   <- sum(alphas[[i]]^2)

    if (i>1){

      stop_condition <- abs(sum_of_sq_betas[[i]] - sum_of_sq_betas[[i-1]])

      s_cond_v[[i]] <- stop_condition

      if (stop_condition < stop_criterium){
        WeCarryOn <- FALSE
      }

    }

    i <- i +1
  }


  result <- list(
    XI = xis,
    ETA = etas,
    ALPHA = alphas,
    BETA= betas,
    nr_iterations = iterations,
#   inverse_of_XIXI = SOLVE_XIXI,
    iterations_crts = iterations_crts,
    sum_absolute_betas = sum_of_sq_betas,
    ridge_penalty = ridge_penaltys,
    nr_nonzeros = nr_nonzeros,
    nr_latent_variables = nr_latent,
    CV_results = CV_results

  )

  result

}
