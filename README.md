# Stability of Topological Clustering in Reproducing Kernel Hilbert Spaces

## Structure:
* COVID_19_Datasets folder consists of (1) daily COVID-related information for each state from January 10, 2021 to January 23, 2021 (Daily counts are updated every day at 4pm EST); (2) weekly COVID-related information for each state (i.e., January 10, 2021 - January 16, 2021 and January 17, 2021 - January 23, 2021).

* Dataset_generation.R includes the procedure for COVID-19 dataset generation from covid19us pacakge.

* TopoKmeans_model.R includes the proposed TopoKmeans function and the codes for reproducing results in Table 1.

* Baselines.R includes the following baselines: (1) TopoCBN, (2) K-means, (3) K-means over vectorized persistence diagram (PD), and (4) Agglomerative Hierarchical clustering.

## Usage
```R
# load a COVID-19 data, e.g., January 11, 2021, from COVID_US_Datasets folder
covid_us_data <- read.csv("COVID_US_Datasets/covid19us_data_01_11_2021.csv",row.names = 1)
# run TopoKmeans model on January 11, 2021 COVID-19 data by setting different hyperparameters
res <- TopoKmeans(covid_us_data,nKNN=38, nClust=4, power=15, sigma =20, preserveOrdering=FALSE, null_dim = TRUE, first_dim = FALSE) 
# compute silhouette information according to a given clustering in 4 clusters
topokmeans_SIL <- cluster::silhouette(res$results$Cls, res$dist_object)
# compute the average silhouette coefficient
avg_topokmeans_SIL <- mean(topokmeans_SIL[, 3])
```

## Requirements
R version 3.6.3; Packages: covid19us, cluster, kernlab, TDA, FCPS, and FNN.
