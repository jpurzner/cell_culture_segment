
require(dplyr)

getcluster_mixtools <- function(data_frame, model, noise_cluster = 5){


  
# Determine the maximum posterior probability for each observation
max_posterior <- apply(model$posterior, 1, max)


cluster_assignments <- apply(model$posterior, 1, which.max)



temp <- data.frame(cluster =  cluster_assignments, max_post = max_posterior, cell_id = row.names(model$x))


temp <- temp %>%
  dplyr::mutate(cluster_mod = ifelse(max_post <= 0.4, noise_cluster, cluster) )


data_frame <- left_join(data_frame, temp, by = "cell_id")
data_frame$dapi_area_cluster <- data_frame$cluster_mod


plot <- ggplot(data_frame, aes(x = area, y = DAPI_median, color = as.factor(dapi_area_cluster))) +
  geom_point(alpha = 0.1, size = 0.1) + # Adjust alpha for point transparency if needed
  scale_color_viridis_d() + # A different color palette for distinction 
  xlim(0, 2000) +
  ylim(3.0, 4.2) +
  theme_minimal() +
  labs(color = "Cluster Assignment") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme(legend.position = c(0.99, 0.99), # Position the legend inside the plot
        legend.justification = c(1, 1), # Anchor point for the legend position
        legend.direction = "vertical", # Align legend items vertically
        legend.box.just = "right", # Justify the legend box
        legend.background = element_rect(fill = "transparent", colour = NA)) 



print(plot)

return(data_frame)

}