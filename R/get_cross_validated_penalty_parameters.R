get_cross_validated_penalty_parameters <- function(predictor,
                                                   predicted,
                                                   penalization,
                                                   ridge_penalty,
                                                   nonzero,
                                                   nr_subsets,
                                                   max_iterations,
                                                   tolerance,
                                                   parallel_CV = parallel_CV
                                                   ){

  shuffled <-  get_split_sets(X = predictor,
                              Y = predicted,
                              nr_subsets = nr_subsets)

  X.sampled     <-   shuffled$X.sampled
  Y.sampled     <-   shuffled$Y.sampled
  label         <-   shuffled$labels

  if (parallel_CV){
    cv_results <- get_parallel_cv(X = X.sampled,
                                  Y = Y.sampled,
                                  lambdas = ridge_penalty,
                                  non_zeros = nonzero,
                                  label = label,
                                  penalization = penalization,
                                  max_iterations = max_iterations,
                                  tolerance = tolerance)

  } else {
    cv_results <- get_non_parallel_cv(X = X.sampled,
                        Y = Y.sampled,
                        lambdas = ridge_penalty,
                        non_zeros = nonzero,
                        label = label,
                        penalization = penalization,
                        max_iterations = max_iterations,
                        tolerance = tolerance)

  }

    a = cv_results$mean_abs_cors[,3]

  best_values     <- cv_results$mean_abs_cors[which.max(a),]

  best_ridge   <- best_values[1]
  best_nonzero   <- best_values[2]


  #**********************
  result <- list(
                 abs_cors = cv_results$abs_cors,
                 mean_abs_cors = cv_results$mean_abs_cors,
                 stime = cv_results$stime,
                 iterations_m = cv_results$iterations_m,
                 best_ridge = best_ridge,
                 best_nonzero = best_nonzero
  )


  class(result) = "CVsRDA"

  result
}

#'@S3method print CVsRDA
print.CVsRDA <- function(x, ...)
{
  cat("Cross validation (CV) for sRDA or sCCA", "\n")
  cat(rep("-", 45), sep="")
  cat("\n   NAME                ", "DESCRIPTION")
  cat("\n1  $abs_cors           ", "sum of absolute correlations per k-th fold CV")
  cat("\n2  $mean_abs_cors      ", "mean absolute correlation per CV")
  cat("\n3  $best_ridge         ", "best ridge parameter selected for the model")
  cat("\n4  $best_nonzero       ", "best number of nonzero alpha weights selected")
  cat("\n5  $stime              ", "time elapsed in seconds")
  cat("\n6  $iterations_m       ", "number of iterations ran before convergence")
  cat("\n")
  cat(rep("-", 45), sep="")
  invisible(x)
}
