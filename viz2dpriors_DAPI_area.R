

starting_covariance_matrices <- list(
  matrix(c(var_area, 0, 0, var_DAPI_median), ncol = 2),
  matrix(c(var_area, 0, 0, var_DAPI_median), ncol = 2),
  matrix(c(var_area, 0, 0, var_DAPI_median), ncol = 2),
  matrix(c(var_area, 0, 0, var_DAPI_median), ncol = 2)
)

scaling_matrices <- list(
  matrix(c(1/4, 0, 0, 3), ncol = 2),  # Scaling matrix for the first covariance matrix
  matrix(c(1/4, 0, 0, 3), ncol = 2),  # Scaling matrix for the second covariance matrix
  matrix(c(1/4, 0, 0, 2), ncol = 2),  # Scaling matrix for the third covariance matrix
  matrix(c(1/4, 0, 0, 2), ncol = 2)   # Scaling matrix for the fourth covariance matrix
)

