library(ggplot2)
library(MASS) # For mvrnorm, if needed

viz2dpriors_area <- function(data_frame, area = c(420, 650, 980, 1050), dapi_mod = c(0.2, 0, 0, -0.15), scaling_matrices = NULL, rotation_values = NULL ) {
  
  
  # the rationale behind this is explained here 
  # https://www.visiondummy.com/2014/04/geometric-interpretation-covariance-matrix/
  
  # area 4 values, dead cell, G0/G1, G0 large and G2 
  # dapi_mod adds or substracts the value from the bulk mean value 
  
  # Assuming 'area' and 'DAPI_median' are column names in your data_frame
  var_area <- var(data_frame$area)
  var_DAPI_median <- var(data_frame$DAPI_median)
  
  # Assuming a negative covariance for illustration
  negative_covariance <- -sqrt(var_area * var_DAPI_median) / 2
  
  # calculate the mean of the dapi channel as this can cary
  DAPI_mu_all <- mean(data_frame$DAPI_median)
  
  # Example initial means
  initial_mu <- matrix(c(area[1], DAPI_mu_all + dapi_mod[1], area[2], DAPI_mu_all + dapi_mod[2], area[3], DAPI_mu_all + dapi_mod[3], area[4], DAPI_mu_all + dapi_mod[4] ), ncol = 2, byrow = TRUE)
  colnames(initial_mu) <- c("area", "DAPI_median")
  initial_mu_df <- as.data.frame(initial_mu)
  
  original_covariance_matrices <- list(
    matrix(c(var_area/4, -20, -20, var_DAPI_median*3), ncol = 2),
    matrix(c(var_area/2, negative_covariance -15, negative_covariance-15, var_DAPI_median*3), ncol = 2),
    matrix(c(var_area/2, negative_covariance-10, negative_covariance-10, var_DAPI_median*2), ncol = 2),
    matrix(c(var_area*2, negative_covariance-10, negative_covariance-10, var_DAPI_median/2), ncol = 2)
  )
  
  starting_covariance_matrices <- list(
    matrix(c(var_area, 0, 0, var_DAPI_median), ncol = 2),
    matrix(c(var_area, 0, 0, var_DAPI_median), ncol = 2),
    matrix(c(var_area, 0, 0, var_DAPI_median), ncol = 2),
    matrix(c(var_area, 0, 0, var_DAPI_median), ncol = 2)
  )
  
  if(is.null(scaling_matrices)) {
  # Define covariance matrices for each initial mean
    scaling_matrices <- list(
      matrix(c(1/4, 0, 0, 3), ncol = 2),  # Scaling matrix for the first covariance matrix
      matrix(c(1/4, 0, 0, 3), ncol = 2),  # Scaling matrix for the second covariance matrix
      matrix(c(1/4, 0, 0, 2), ncol = 2),  # Scaling matrix for the third covariance matrix
      matrix(c(1/4, 0, 0, 2), ncol = 2)   # Scaling matrix for the fourth covariance matrix
    )
    } else {
  #    if (length(scaling_matrices) == length(starting_covariance_matrices) { 
#TODO test if the matrix size and the list size makes sense
      # 
#      }
    
  }
  
  covariance_matrices <- Map(function(cov_matrix, scale_matrix) cov_matrix * scale_matrix, starting_covariance_matrices, scaling_matrices)
  rotate_covariance_matrices <- rotate_covariance_matrices(covariance_matrices, rotation_values)
  covariance_matrices <- rotate_covariance_matrices
  
  # Generate ellipse points for each cluster
  ellipse_points_list <- lapply(1:4, function(i) generate_ellipse_points(initial_mu[i, ], rotate_covariance_matrices[[i]]))
  
  # Generate noise component
  noise_mu <- c(mean(data_frame$area) + 0.2, mean(data_frame$DAPI_median) + 0.0)
  noise_variance_area <-14 * var(data_frame$area)
  noise_variance_DAPI <- 14 *  var(data_frame$DAPI_median)
  noise_sigma <- matrix(c(noise_variance_area, 0, 0, noise_variance_DAPI), ncol = 2)
  ellipse_noise <- generate_ellipse_points(noise_mu, noise_sigma)
  colnames(ellipse_noise) <- c("area", "DAPI_median")
  
  # Combine initial means and sigma lists with noise component
  initial_mu_combined <- rbind(initial_mu, noise_mu)
  initial_sigma_list_combined <- c(covariance_matrices, list(noise_sigma))
  initial_mu_combined_list <- lapply(1:nrow(initial_mu_combined), function(i) initial_mu_combined[i, ])
  # Plotting
#  plot <- ggplot() +
#    geom_point(data = data_frame, aes(x = area, y = DAPI_median), alpha = 0.1, size = 0.1) +
#    geom_point(data = initial_mu_df, aes(x = area, y = DAPI_median), color = "red", size = 3, alpha = 1) +
#    lapply(ellipse_points_list, function(ep) geom_polygon(data = as.data.frame(ep), aes(x = area, y = DAPI_median), fill = NA, color = "blue")) +
#    geom_polygon(data = as.data.frame(ellipse_noise), aes(x = area, y = DAPI_median), fill = NA, color = "blue") +
#    theme_minimal()
 
  # Start the plot
  plot <- ggplot(data_frame, aes(x = area, y = DAPI_median)) +
    geom_point(alpha = 0.1, size = 0.05) +
    geom_point(data = initial_mu_df, color = "red", size = 3, alpha = 1) +
    theme_minimal()
  
  # Add ellipse layers individually
  for (ep in ellipse_points_list) {
    plot <- plot + geom_polygon(data = as.data.frame(ep), aes(x = area, y = DAPI_median), fill = NA, color = "blue")
  }
  
  # Add noise ellipse
  plot <- plot + geom_polygon(data = as.data.frame(ellipse_noise), aes(x = area, y = DAPI_median), fill = NA, color = "blue")
  
  
  #plot <- plot + xlim(50, 2200) + ylim(2.6, 4.5)
  # Print the plot
  print(plot)
  
  # Return initial means and sigma lists
  return(list(initial_mu_list = initial_mu_combined_list, initial_sigma_list = initial_sigma_list_combined))

  } # end of main function 



# rotate_covariance_matrices: Rotates a list of covariance matrices by specified angles.
#
# This function takes a list of covariance matrices and a vector of rotation angles (in degrees),
# and applies a rotation to each covariance matrix according to its corresponding angle. The function
# ensures that each covariance matrix is rotated individually, allowing for a unique rotation
# for each matrix based on the angles provided in the rotation_angle_degrees vector.
#
# Parameters:
#   covariance_matrices: A list of 2x2 covariance matrices. Each matrix in the list represents
#                        the covariance matrix of a cluster or group that is to be rotated.
#
#   rotation_angle_degrees: A numeric vector containing rotation angles in degrees. Each element
#                           in this vector corresponds to a covariance matrix in the list, specifying
#                           the angle by which that particular matrix should be rotated. The length
#                           of this vector must match the length of the covariance_matrices list.
#
# The function checks to ensure the lengths of the covariance_matrices list and the rotation_angle_degrees
# vector match. It then iterates over each covariance matrix and its corresponding rotation angle,
# converting the angle from degrees to radians and constructing a 2D rotation matrix. This rotation
# matrix is then applied to the covariance matrix, effectively rotating it by the specified angle.
# The rotated covariance matrices are returned as a list, with each element being a rotated version
# of the input covariance matrices.
#
# Returns:
#   A list of rotated covariance matrices, where each matrix has been rotated by its corresponding
#   angle specified in rotation_angle_degrees.
#
# Example usage:
#   rotated_matrices <- rotate_covariance_matrices(list_of_cov_matrices, c(30, 45, 60))
#   This example rotates the first covariance matrix by 30 degrees, the second by 45 degrees,
#   and the third by 60 degrees.

rotate_covariance_matrices <- function(covariance_matrices, rotation_angle_degrees) {
  if (length(covariance_matrices) != length(rotation_angle_degrees)) {
    stop("The lengths of covariance_matrices and rotation_angle_degrees must match.")
  }
  
  rotated_covariance_matrices <- mapply(function(cov_matrix, angle_degree) {
    rotation_angle_radians <- angle_degree * (pi / 180)  # Convert degrees to radians
    rotation_matrix <- matrix(c(cos(rotation_angle_radians), -sin(rotation_angle_radians),
                                sin(rotation_angle_radians), cos(rotation_angle_radians)), nrow = 2, byrow = TRUE)
    
    rotated_cov_matrix <- rotation_matrix %*% cov_matrix %*% t(rotation_matrix)
    return(rotated_cov_matrix)
  }, cov_matrix = covariance_matrices, angle_degree = rotation_angle_degrees, SIMPLIFY = FALSE)
  
  return(rotated_covariance_matrices)
}

  
# Function to generate ellipse points
generate_ellipse_points <- function(mean, covariance, n_points = 100) {
  eigen_values <- eigen(covariance)$values
  eigen_vectors <- eigen(covariance)$vectors
  theta <- seq(0, 2*pi, length.out = n_points)
  # Ellipse points in the principal component space
  ellipse_points <- t(eigen_vectors %*% diag(sqrt(eigen_values)) %*% rbind(cos(theta), sin(theta)))
  # Adjust points to the correct location
  t(apply(ellipse_points, 1, function(x) x + mean))
}