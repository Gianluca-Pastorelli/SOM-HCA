## Load relevant packages.
library(dplyr) # For dataframe manipulation.
library(kohonen) # For self-organising maps (SOM).
library(aweSOM) # For calculating SOM map quality.
library(maptree) # For hierarchical cluster analysis (HCA).

## Read the data.
# Data must be stored in a single data matrix, saved as a CSV file.
# The matrix should consist of multiple observations (rows) representing measured samples or points,
# and a set of variables (columns) associated with X-ray energies.
# Except for the column headings, all cells must only contain numeric values.
# The net intensities for each variable may or may not be normalised within their respective ranges.
My_ds <- as.matrix(read.csv('C:\\File_path\\...\\File_name.csv')) # Specify the path and file name.

## Define the function optimalSOM() to estimate the SOM training grid size based on the data.
optimalSOM <- function(method, increments, iterations) {
  # Calculate the maximum dimension of a square grid size using one of the approaches indicated by the argument "method".
  if(method=="A") { 
    max_dim <- round(sqrt(5*sqrt(nrow(My_ds)))) # Max dimension is determined with the heuristic formula by Vesanto et al.
  }
  else if(method=="B") {
    max_dim <- round((nrow(My_ds))^(1/2.5)) # Max dimension is determined using a slightly different approach.
  }
  else {max_dim <- method # Max dimension is a value manually selected.
  }
  # Execute a growing SOM iterative process to calculate map quality and determine the optimal grid size.
  error_df <-data.frame(NA, NA, NA, NA, NA)
  names(error_df) <-c("Dimension","Quantisation_error","Topographic_error","Kaski-Lagus_error","Explained_variance")
  seq <- seq(2,max_dim,increments)
  progbar <- txtProgressBar(min=0,max=max_dim,style=3)
  set.seed(281122)
  for (i in seq) {
    My_Grid <- somgrid(xdim = i, ydim = i, topo = "hexagonal", toroidal = T)
    My_Model <- som(X = My_ds,
                    grid = My_Grid, 
                    rlen=iterations,
                    alpha=c(0.05,0.01),
                    keep.data = TRUE)
    quant.error <- as.numeric((aweSOM::somQuality(My_Model, My_ds))$err.quant)
    topo.error <- as.numeric((aweSOM::somQuality(My_Model, My_ds))$err.topo)
    kaski.error <- as.numeric((aweSOM::somQuality(My_Model, My_ds))$err.kaski)
    varratio.error <- as.numeric((aweSOM::somQuality(My_Model, My_ds))$err.varratio)
    error_df <- error_df %>% add_row(Dimension = i,
                                     `Quantisation_error` = quant.error,
                                     `Topographic_error` = topo.error,
                                     `Kaski-Lagus_error` = kaski.error,
                                     `Explained_variance` = varratio.error)
    setTxtProgressBar(progbar,value=i)
  }
  close(progbar)
  error_df <- error_df %>%
    na.omit() %>%
    mutate(across(c(Quantisation_error, Topographic_error, `Kaski-Lagus_error`, Explained_variance), scale)) %>%
    mutate(QplusT_error = Quantisation_error + Topographic_error,
           QplusTplusK_error = Quantisation_error + Topographic_error + `Kaski-Lagus_error`,
           all_error = Explained_variance - QplusTplusK_error)
  QTerror_df <- error_df %>% filter(QplusT_error == min(QplusT_error))
  Kerror_df <- error_df %>% filter(`Kaski-Lagus_error` == min(`Kaski-Lagus_error`))
  QTKerror_df <- error_df %>% filter(QplusTplusK_error == min(QplusTplusK_error))
  Verror_df <- error_df %>% filter(Explained_variance == max(Explained_variance))
  All_errors_df <- error_df %>% filter(all_error == max(all_error))
  My_dim_QT <- as.numeric(QTerror_df$Dimension)
  My_dim_K <- as.numeric(Kerror_df$Dimension)
  My_dim_QTK <- as.numeric(QTKerror_df$Dimension)
  My_dim_V <- as.numeric(Verror_df$Dimension)
  My_dim_all <- as.numeric(All_errors_df$Dimension)
  return(data.frame("[Quality measure]"=c("Min nQe+nTe","Min nKLe","Min nQe+nTe+nKLe","Max n%ev", "Max QI"),
                    "[Value]"=c(QTerror_df$QplusT_error, Kerror_df$`Kaski-Lagus_error`, QTKerror_df$QplusTplusK_error,
                            Verror_df$Explained_variance, All_errors_df$all_error),
                    "[Associated grid dimension]"=c(My_dim_QT, My_dim_K, My_dim_QTK, My_dim_V, My_dim_all),
                    check.names = FALSE))
  }

## Call the optimalSOM(method, increments, iterations) function
## and assess the resulting ideal grid dimension based on various quality measures.
# Set the "method" argument to either "A" (Vesanto) or "B", or to a specific numeric value.
# Adjust the "increments" argument to a desired value, such as increasing by 2 or 5 rows/columns at each step.
# Set the "iterations" argument to a relatively low value, e.g., less than 500, to reduce computation time
# (an error message could suggest that reducing the number of iterations might be necessary).
optimalSOM("A",2,300)

## Define the function finalSOM() to re-train the SOM model using the selected optimal grid size.
finalSOM <- function(dimension, iterations) {
  My_dim <- dimension
  My_Grid <- somgrid(xdim = My_dim, ydim = My_dim, topo = "hexagonal", toroidal = T)
  # Create the model
  set.seed(291122)
  My_Model <- som(X = My_ds,
                  grid = My_Grid, 
                  rlen=iterations, 
                  alpha=c(0.05,0.01),
                  keep.data = TRUE)
  return(My_Model)
}

## Call the finalSOM(dimension, iterations) function and store the resulting model.
# Set the "dimension" argument to one of the previously calculated optimal dimension values,
# or choose another value manually if desired.
# Set the "iterations" argument to a large value, e.g., 500 or higher, for improved training
# (an error message could suggest that reducing the number of iterations might be necessary).
My_Model <- finalSOM(6,700)

## Generate different types of plots.
# Training progress for SOM.
plot(My_Model, type="changes")
# Node count plot (for map quality).
plot(My_Model, type="count", main="Node Counts")
# U-matrix visualization (similarity between each node and its neighbors).
plot(My_Model, type="dist.neighbours", main = "SOM neighbour distances")
# Weight vector view (patterns in the distribution of samples and variables).
plot(My_Model, type="codes")
# Kohonen heatmap creation (distribution of single variables across the map);
# if there are improper column names, it could potentially result in an error.
for (i in 1:ncol(My_ds)){
  plot(My_Model, type = "property", property = getCodes(My_Model)[,i], main=colnames(getCodes(My_Model))[i])
}

## Perform clustering on the SOM nodes to group similar units.
# Use HCA to cluster the codebook vectors.
set.seed(231122)
distance <- dist(getCodes(My_Model))
clustering <- hclust(distance)
# The KGS penalty function allows the optimal number of clusters to be estimated automatically by selecting the lowest penalty.
optimal_k <- kgs(clustering, distance, maxclus = 20)
clusters <- as.integer(names(optimal_k[which(optimal_k == min(optimal_k))]))
cat(clusters,"clusters were determined","\n") # Output the calculated number of clusters
som_cluster <- cutree(clustering, clusters) # Instead of "clusters", a specific number of clusters can be chosen and used
# Plot the result.
pretty_palette <- c("#1f77b4", '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2') # Define the colour palette
plot(My_Model, type="mapping", bgcol = pretty_palette[som_cluster], main = "Clusters") 
add.cluster.boundaries(My_Model, som_cluster)

## Assign the generated clusters to the original observations in the data.
# Get the vector with the cluster value for each original observation.
cluster_assignment <- som_cluster[My_Model$unit.classif]
# Add the assigned clusters as a column to a copy of the original data.
My_ds2 <- as.data.frame(My_ds)
My_ds2$cluster <- cluster_assignment
# Export the copy of the original data with the assigned clusters.
write.csv(My_ds2, file='C:\\File_path\\...\\New_file_name.csv', row.names=FALSE) # Specify the path and file name.
