# Self-Organising Maps coupled with Hierarchical Cluster Analysis

## Overview
This R script utilizes the kohonen, aweSOM and maptree packages to perform a SOM-HCA analysis on a given dataset.
The script includes functions to estimate the optimal SOM grid size based on various quality measures and subsequently generates a SOM model with the selected dimensions.
It also performs hierarchical clustering on the SOM nodes to group similar units.
Make sure to adapt the script according to your specific dataset and analysis.

## Instructions
### 1. Load Relevant Packages
The script requires the dplyr, kohonen, aweSOM, and maptree packages. Make sure to install these packages before running the script.

### 2. Read the Data:
Provide the path and file name of your dataset in CSV format.

### 3. Estimate Optimal SOM Grid Size:
Use the optimalSOM function to estimate the optimal SOM grid size based on specified criteria.

### 4. Train the Final SOM Model:
Use the finalSOM function to train the SOM model with the selected optimal grid size.

### 5. Generate Plots:
Create various plots to visualize SOM training progress and results.

### 6. Perform Clustering on SOM Nodes:
Utilize hierarchical clustering and KGS penalty function to determine the optimal number of clusters.

### 7. Assign Clusters to Original Observations:
Add the assigned clusters as a column to a copy of the original data and export the results.
