require(ggplot2)
require(dplyr)

viz2dpriors <- function(data_frame, col_x = NULL, col_y = NULL, values_x = NULL, values_y = NULL, 
                        mods_x = rep(0, length(values_x)), mods_y = rep(0, length(values_y)), 
                        scaling_matrices = NULL, rotation_values = NULL, 
                        add_noise = TRUE, noise_stretch = 6, noise_mod_x= 0, noise_mod_y = 0) {
  
  # Automatically assign col_x and col_y if not specified
  if(is.null(col_x)) {
    col_x <- names(data_frame)[1]
  }
  
  if(is.null(col_y)) {
    col_y <- names(data_frame)[2]
  }
  # Calculate the mean of the data along each axis
  mean_x <- mean(data_frame[[col_x]], na.rm = TRUE)
  mean_y <- mean(data_frame[[col_y]], na.rm = TRUE)
  
  # Use values_x and values_y if provided, otherwise adjust based on mods_x and mods_y
  if(is.null(values_x)) {
    adjusted_values_x <- mean_x + mods_x
  } else {
    adjusted_values_x <- values_x
  }
  
  if(is.null(values_y)) {
    adjusted_values_y <- mean_y + mods_y
  } else {
    adjusted_values_y <- values_y
  }
  
  # Prepare initial means based on adjusted values
  initial_mu <- cbind(adjusted_values_x, adjusted_values_y)
  colnames(initial_mu) <- c(col_x, col_y)
  initial_mu_df <- as.data.frame(initial_mu)
  
  
  var_x <- var(data_frame[[col_x]])
  var_y <- var(data_frame[[col_y]])
  
  if(is.null(scaling_matrices)) {
    scaling_matrices <- lapply(1:length(adjusted_values_y), function(i) {
      matrix(c(1, 0, 0, 1), ncol = 2) # Default identity scaling if not provided
    })
  }
  
  if(is.null(rotation_values)) {
    rotation_values <- rep(0, length(adjusted_values_y)) # No rotation if not provided
  }
  
  # Check if the length of scaling_matrices and rotation_values match the number of clusters
  if(length(scaling_matrices) != length(adjusted_values_y) || length(rotation_values) != length(adjusted_values_y)) {
    stop("The lengths of scaling_matrices and rotation_values must match the number of clusters.")
  }
  
  # Apply scaling and rotation
  covariance_matrices <- lapply(scaling_matrices, function(scale_matrix) {
    matrix(c(var_x, 0, 0, var_y), ncol = 2) * scale_matrix
  })
  
  rotated_covariance_matrices <- rotate_covariance_matrices(covariance_matrices, rotation_values)
  
  ellipse_points_list <- lapply(1:length(rotated_covariance_matrices), function(i) {
    generate_ellipse_points(initial_mu[i, ], rotated_covariance_matrices[[i]])
  })
  
  # Start the plot
  plot <- ggplot(data_frame, aes_string(x = col_x, y = col_y)) +
    geom_point(alpha = 0.1, size = 0.1) +
    geom_point(data = initial_mu_df, aes_string(x = col_x, y = col_y), color = "red", size = 3, alpha = 1) +
    theme_minimal()
  
  
  # Add ellipse layers individually
  for (ep in ellipse_points_list) {
    plot <- plot + geom_polygon(data = as.data.frame(ep), aes_string(x = col_x, y = col_y), fill = NA, color = "blue")
  }
  
  
  if (add_noise) {
    # Generate noise component
    # TODO if no noise modifier provided that make sure the noise surrounds the mu values 
    
    noise_mu <- c(mean(data_frame[,col_x]) + noise_mod_x  , mean(data_frame[,col_y]) + noise_mod_y )
    noise_variance_x <-noise_stretch * var(data_frame[,col_x])
    noise_variance_y <- noise_stretch*  var(data_frame[,col_y])
    noise_sigma <- matrix(c(noise_variance_x, 0, 0, noise_variance_y), ncol = 2)
    ellipse_noise <- generate_ellipse_points(noise_mu, noise_sigma)
    colnames(ellipse_noise) <- c(col_x, col_y)
    
    plot <- plot + geom_polygon(data = as.data.frame(ellipse_noise), aes_string(x = col_x, y = col_y), fill = NA, color = "blue")
    
    initial_mu <- rbind(initial_mu, noise_mu)
    initial_sigma_list <- c(covariance_matrices, list(noise_sigma))
    
  }
  
  print(plot)
  
  # Combine initial means and sigma lists with noise component
  
  
  initial_mu_list <- lapply(1:nrow(initial_mu), function(i) initial_mu[i, ])
  
  
  # estimate the lambda values based on the number of points within the elipse 
  # estimate the lambda values based on the number of points within the ellipse
  points_within_ellipse <- function(data, mean, covariance) {
    mahalanobis_distance <- apply(data, 1, function(x) {
      (t(x - mean) %*% solve(covariance) %*% (x - mean))
    })
    threshold <- qchisq(0.95, df = 2) # 95% confidence level
    sum(mahalanobis_distance <= threshold)
  }
  
  # Apply the function to each cluster and calculate the number of points within each ellipse
  counts_within_each_ellipse <- sapply(1:length(adjusted_values_x), function(i) {
    data_subset <- data_frame[, c(col_x, col_y)]
    points_within_ellipse(as.matrix(data_subset), initial_mu[i, ], rotated_covariance_matrices[[i]])
  })
  
  # Return initial means and sigma lists
  return(list(initial_mu_list = initial_mu_list, initial_sigma_list = initial_sigma_list))
  
}



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