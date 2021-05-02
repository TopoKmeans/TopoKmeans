require(kernlab)
require(TDA)
require(FCPS)
require(FNN)

TopoKmeans = function(data, nKNN, nClust=2, power = 5, sigma=0.05,dist_matrix=FALSE, preserveOrdering=FALSE, null_dim = FALSE, first_dim= FALSE){
  # main body
  print('Forming PDs (1/3)')
  maxscale = ceiling(max(apply(data, 2, max))) + 1
  
  N <- nrow(data) # number of data points
  if (dist_matrix){
    ind <- t(apply(data,1,order))[,1:nKNN]
    maxscale <- max(t(apply(data,1,sort))[,nKNN])
  } else{
    if(preserveOrdering){
      ind<-list()
      for (i in 1:N){
        ind <- append(ind, list(max(1,i-floor(nKNN/2)):min(N,i+floor(nKNN/2))))
      }
      maxscale<-max(data)
    }
    else{
      ind <- cbind(1:N,knn.index(data,k=nKNN-1)) # indices of k nearest neighbors
      maxscale<-max(knn.dist(data,k=nKNN-1)[,nKNN-1]) # maxscale
    }
  }
  # 1. VR-Complex Filtration #
  dimM = dim0 = dim1 = numeric()

  for (i in 1:N){
    
    if (dist_matrix){
      ripsFltr<-ripsFiltration(data[ind[i,],ind[i,]],maxdimension = 0,maxscale = maxscale,dist = 'arbitrary') 
    } else{
      if(preserveOrdering){
        ripsFltr<-ripsFiltration(matrix(data[ind[[i]],]),maxdimension = 0,maxscale = maxscale,dist = 'euclidean') 
      }
      else{
        ripsFltr<-ripsFiltration(data[ind[i,],],maxdimension = 1,maxscale = maxscale,dist = 'euclidean') 
      }
    }
    # compute PD of Rips filtration
    PD<-filtrationDiag(filtration = ripsFltr, maxdimension = 1)$diagram
    PD[PD[,3]==Inf,3]=maxscale
    dimM[i] = nrow(PD) #67
    
    # pot
    # Pdata0 = PD[PD[,1]==0,2:3]
    # Pdata1 = PD[PD[,1]==1,2:3]
    # pot 2
    dim0[i] = length(which(PD[,1]==0))
    dim1[i] = length(which(PD[,1]==1))    
  }
  
  perst.val_0 = matrix(0, ncol = N, nrow = max(dim0)) 
  perst.val_1 = matrix(0, ncol = N, nrow = max(dim1))
  perst.val = matrix(0, ncol = N, nrow = max(dimM))
  
  for (i in 1:N){
    if (dist_matrix){
      ripsFltr<-ripsFiltration(data[ind[i,],ind[i,]],maxdimension = 0,maxscale = maxscale,dist = 'arbitrary') 
    } else{
      if(preserveOrdering){
        ripsFltr<-ripsFiltration(matrix(data[ind[[i]],]),maxdimension = 0,maxscale = maxscale,dist = 'euclidean') 
      }
      else{
        ripsFltr<-ripsFiltration(data[ind[i,],],maxdimension = 1,maxscale = maxscale,dist = 'euclidean') 
      }
    }
    
    # compute PD of Rips filtration
    PD<-filtrationDiag(filtration = ripsFltr, maxdimension = 1)$diagram
    # replace death=Inf with death=maxscale
    PD[PD[,3]==Inf,3]=maxscale
    aad_dimM = nrow(PD)
    
    # pot 2
    aad_dim_1 = length(which(PD[,1]==0))
    aad_dim_2 = length(which(PD[,1]==1))
    # pot 3 
    if  (aad_dim_1>0){
      if  (aad_dim_1>1){
        perst.val_0[1:aad_dim_1, i] = PD[PD[,1]==0,2:3][,2] - PD[PD[,1]==0,2:3][,1]
      }
      else{
        perst.val_0[1:aad_dim_1, i] = PD[PD[,1]==0,2:3][2] - PD[PD[,1]==0,2:3][1]
      }
    }
    if  (aad_dim_2>0){
      if (aad_dim_2>1){
        perst.val_1[1:aad_dim_2, i] = PD[PD[,1]==1,2:3][,2] - PD[PD[,1]==1,2:3][,1]
      }
      else{
        perst.val_1[1:aad_dim_2, i] = PD[PD[,1]==1,2:3][2] - PD[PD[,1]==1,2:3][1]
      }
    }
    perst.val[1:aad_dimM, i] = PD[,3] - PD[,2]
    
  }
  
  ############################
  # Forming distance matrix
  ############################
  print('Forming distance matrix (2/3)')
  # compute the distance matrix between kernalized PDs
  N<-dim(perst.val)[1]
  n<-dim(perst.val)[2]
  rbf <- rbfdot(sigma = sigma) # rbf = laplacedot(sigma = sigma) ; #rbfdot(sigma = sigma)
  dist<-matrix(0, nr=n, nc=n)
  
  if(null_dim){
    persistence = perst.val_0
  }else if (first_dim){
    persistence = perst.val_1
  }else{
    persistence = perst.val
  }
  print(dim(persistence))
  kM<-rep(0,n)
  for (i in 1:n){
    kM[i]<-sum(kernelMatrix(rbf, persistence[,i]))
  }
  
  for (i in 1:n){
    print(i)
    for (j in i:n){
      if (i!=j){
        dist[i,j]<- (kM[i] + kM[j] -2*sum(kernelMatrix(rbf, x=persistence[,i], y=persistence[,j])))^{power}
        dist[j,i]<-dist[i,j]
      }
    }
  }
  
  
  ############################
  # Performing clustering
  ############################
  print('Performing clustering (3/3)')
  res<-kmeansClustering(dist, ClusterNo=nClust,PlotIt=TRUE,Verbose = T) # Verbose = F
  
  return(list(results=res, persistence=persistence, dist_object = dist))
}


# Peformances of TopoKmeans in Table 1
# 01/10/2021
# based on Euclidean distance #
covid_us_data <- read.csv("COVID_US_Datasets/covid19us_data_01_10_2021.csv",row.names = 1)
res <- TopoKmeans(covid_us_data,nKNN=38, nClust=4, power=15, sigma =20,
                  preserveOrdering=FALSE, null_dim = TRUE, first_dim = FALSE) 
topokmeans_SIL <- cluster::silhouette(res$results$Cls, dist(covid_us_data, method = "euclidean"))
avg_topokmeans_SIL <- mean(topokmeans_SIL[, 3])
avg_topokmeans_SIL

# based on Topological distance #
covid_us_data <- read.csv("COVID_US_Datasets/covid19us_data_01_10_2021.csv",row.names = 1)
res <- TopoKmeans(covid_us_data,nKNN=38, nClust=4, power=15, sigma =20,
                  preserveOrdering=FALSE, null_dim = TRUE, first_dim = FALSE) 
topokmeans_SIL <- cluster::silhouette(res$results$Cls, res$dist_object)
avg_topokmeans_SIL <- mean(topokmeans_SIL[, 3])
avg_topokmeans_SIL


# ************************************* #
# 01/11/2021
# based on Euclidean distance #
covid_us_data <- read.csv("COVID_US_Datasets/covid19us_data_01_11_2021.csv",row.names = 1)
res <- TopoKmeans(covid_us_data,nKNN=48, nClust=4, power=15, sigma =20, 
                  preserveOrdering=FALSE, null_dim = TRUE, first_dim = FALSE) 
topokmeans_SIL <- cluster::silhouette(res$results$Cls, dist(covid_us_data, method = "euclidean"))
avg_topokmeans_SIL <- mean(topokmeans_SIL[, 3])
avg_topokmeans_SIL

# based on Topological distance #
covid_us_data <- read.csv("COVID_US_Datasets/covid19us_data_01_11_2021.csv",row.names = 1)
res <- TopoKmeans(covid_us_data,nKNN=48, nClust=4, power=15, sigma =20, 
                  preserveOrdering=FALSE, null_dim = TRUE, first_dim = FALSE) 
topokmeans_SIL <- cluster::silhouette(res$results$Cls, res$dist_object)
avg_topokmeans_SIL <- mean(topokmeans_SIL[, 3])
avg_topokmeans_SIL


# ************************************* #
# 01/12/2021
# based on Euclidean distance #
covid_us_data <- read.csv("COVID_US_Datasets/covid19us_data_01_12_2021.csv",row.names = 1)
res <- TopoKmeans(covid_us_data,nKNN=38, nClust=4, power=15, sigma =10, 
                  preserveOrdering=FALSE, null_dim = TRUE, first_dim = FALSE) 
topokmeans_SIL <- cluster::silhouette(res$results$Cls, dist(covid_us_data, method = "euclidean"))
avg_topokmeans_SIL <- mean(topokmeans_SIL[, 3])
avg_topokmeans_SIL

# based on Topological distance #
covid_us_data <- read.csv("COVID_US_Datasets/covid19us_data_01_12_2021.csv",row.names = 1)
res <- TopoKmeans(covid_us_data,nKNN=38, nClust=4, power=15, sigma =10, 
                  preserveOrdering=FALSE, null_dim = TRUE, first_dim = FALSE) 
topokmeans_SIL <- cluster::silhouette(res$results$Cls, res$dist_object)
avg_topokmeans_SIL <- mean(topokmeans_SIL[, 3])
avg_topokmeans_SIL


# ************************************* #
# 01/13/2021
# based on Euclidean distance #
covid_us_data <- read.csv("COVID_US_Datasets/covid19us_data_01_13_2021.csv",row.names = 1)
res <- TopoKmeans(covid_us_data,nKNN=42, nClust=4, power=15, sigma =15, 
                  preserveOrdering=FALSE, null_dim = TRUE, first_dim = FALSE) 
topokmeans_SIL <- cluster::silhouette(res$results$Cls, dist(covid_us_data, method = "euclidean"))
avg_topokmeans_SIL <- mean(topokmeans_SIL[, 3])
avg_topokmeans_SIL

# based on Topological distance #
covid_us_data <- read.csv("COVID_US_Datasets/covid19us_data_01_13_2021.csv",row.names = 1)
res <- TopoKmeans(covid_us_data,nKNN=42, nClust=4, power=15, sigma =15, 
                  preserveOrdering=FALSE, null_dim = TRUE, first_dim = FALSE) 
topokmeans_SIL <- cluster::silhouette(res$results$Cls, res$dist_object)
avg_topokmeans_SIL <- mean(topokmeans_SIL[, 3])
avg_topokmeans_SIL


# ************************************* #
# 01/14/2021
# based on Euclidean distance #
covid_us_data <- read.csv("COVID_US_Datasets/covid19us_data_01_14_2021.csv",row.names = 1)
res <- TopoKmeans(covid_us_data,nKNN=42, nClust=4, power=15, sigma =15, 
                  preserveOrdering=FALSE, null_dim = TRUE, first_dim = FALSE) 
topokmeans_SIL <- cluster::silhouette(res$results$Cls, dist(covid_us_data, method = "euclidean"))
avg_topokmeans_SIL <- mean(topokmeans_SIL[, 3])
avg_topokmeans_SIL

# based on Topological distance #
covid_us_data <- read.csv("COVID_US_Datasets/covid19us_data_01_14_2021.csv",row.names = 1)
res <- TopoKmeans(covid_us_data,nKNN=42, nClust=4, power=15, sigma =15, 
                  preserveOrdering=FALSE, null_dim = TRUE, first_dim = FALSE) 
topokmeans_SIL <- cluster::silhouette(res$results$Cls, res$dist_object)
avg_topokmeans_SIL <- mean(topokmeans_SIL[, 3])
avg_topokmeans_SIL


# ************************************* #
# 01/15/2021
# based on Euclidean distance #
covid_us_data <- read.csv("COVID_US_Datasets/covid19us_data_01_15_2021.csv",row.names = 1)
res <- TopoKmeans(covid_us_data,nKNN=38, nClust=4, power=20, sigma =20, 
                  preserveOrdering=FALSE, null_dim = TRUE, first_dim = FALSE) 
topokmeans_SIL <- cluster::silhouette(res$results$Cls, dist(covid_us_data, method = "euclidean"))
avg_topokmeans_SIL <- mean(topokmeans_SIL[, 3])
avg_topokmeans_SIL

# based on Topological distance #
covid_us_data <- read.csv("COVID_US_Datasets/covid19us_data_01_15_2021.csv",row.names = 1)
res <- TopoKmeans(covid_us_data,nKNN=38, nClust=4, power=20, sigma =20, 
                  preserveOrdering=FALSE, null_dim = TRUE, first_dim = FALSE) 
topokmeans_SIL <- cluster::silhouette(res$results$Cls, res$dist_object)
avg_topokmeans_SIL <- mean(topokmeans_SIL[, 3])
avg_topokmeans_SIL


# ************************************* #
# 01/16/2021
# based on Euclidean distance #
covid_us_data <- read.csv("COVID_US_Datasets/covid19us_data_01_16_2021.csv",row.names = 1)
res <- TopoKmeans(covid_us_data,nKNN=38, nClust=4, power=20, sigma =20, 
                  preserveOrdering=FALSE, null_dim = TRUE, first_dim = FALSE) 
topokmeans_SIL <- cluster::silhouette(res$results$Cls, dist(covid_us_data, method = "euclidean"))
avg_topokmeans_SIL <- mean(topokmeans_SIL[, 3])
avg_topokmeans_SIL

# based on Topological distance #
covid_us_data <- read.csv("COVID_US_Datasets/covid19us_data_01_16_2021.csv",row.names = 1)
res <- TopoKmeans(covid_us_data,nKNN=38, nClust=4, power=20, sigma =20, 
                  preserveOrdering=FALSE, null_dim = TRUE, first_dim = FALSE) 
topokmeans_SIL <- cluster::silhouette(res$results$Cls, res$dist_object)
avg_topokmeans_SIL <- mean(topokmeans_SIL[, 3])
avg_topokmeans_SIL


# ************************************* #
# 01/17/2021
# based on Euclidean distance #
covid_us_data <- read.csv("COVID_US_Datasets/covid19us_data_01_17_2021.csv",row.names = 1)
res <- TopoKmeans(covid_us_data,nKNN=38, nClust=4, power=20, sigma =20, 
                  preserveOrdering=FALSE, null_dim = TRUE, first_dim = FALSE) 
topokmeans_SIL <- cluster::silhouette(res$results$Cls, dist(covid_us_data, method = "euclidean"))
avg_topokmeans_SIL <- mean(topokmeans_SIL[, 3])
avg_topokmeans_SIL

# based on Topological distance #
covid_us_data <- read.csv("COVID_US_Datasets/covid19us_data_01_17_2021.csv",row.names = 1)
res <- TopoKmeans(covid_us_data,nKNN=38, nClust=4, power=20, sigma =20, 
                  preserveOrdering=FALSE, null_dim = TRUE, first_dim = FALSE) 
topokmeans_SIL <- cluster::silhouette(res$results$Cls, res$dist_object)
avg_topokmeans_SIL <- mean(topokmeans_SIL[, 3])
avg_topokmeans_SIL


# ************************************* #
# 01/18/2021
# based on Euclidean distance #
covid_us_data <- read.csv("COVID_US_Datasets/covid19us_data_01_18_2021.csv",row.names = 1)
res <- TopoKmeans(covid_us_data,nKNN=39, nClust=4, power=1, sigma =0.5, 
                  preserveOrdering=FALSE, null_dim = TRUE, first_dim = FALSE) 
topokmeans_SIL <- cluster::silhouette(res$results$Cls, dist(covid_us_data, method = "euclidean"))
avg_topokmeans_SIL <- mean(topokmeans_SIL[, 3])
avg_topokmeans_SIL

# based on Topological distance #
covid_us_data <- read.csv("COVID_US_Datasets/covid19us_data_01_18_2021.csv",row.names = 1)
res <- TopoKmeans(covid_us_data,nKNN=39, nClust=4, power=1, sigma =0.5, 
                  preserveOrdering=FALSE, null_dim = TRUE, first_dim = FALSE) 
topokmeans_SIL <- cluster::silhouette(res$results$Cls, res$dist_object)
avg_topokmeans_SIL <- mean(topokmeans_SIL[, 3])
avg_topokmeans_SIL


# ************************************* #
# 01/19/2021
# based on Euclidean distance #
covid_us_data <- read.csv("COVID_US_Datasets/covid19us_data_01_19_2021.csv",row.names = 1)
res <- TopoKmeans(covid_us_data,nKNN=39, nClust=4, power=20, sigma =20, 
                  preserveOrdering=FALSE, null_dim = TRUE, first_dim = FALSE) 
topokmeans_SIL <- cluster::silhouette(res$results$Cls, dist(covid_us_data, method = "euclidean"))
avg_topokmeans_SIL <- mean(topokmeans_SIL[, 3])
avg_topokmeans_SIL

# based on Topological distance #
covid_us_data <- read.csv("COVID_US_Datasets/covid19us_data_01_19_2021.csv",row.names = 1)
res <- TopoKmeans(covid_us_data,nKNN=39, nClust=4, power=20, sigma =20, 
                  preserveOrdering=FALSE, null_dim = TRUE, first_dim = FALSE) 
topokmeans_SIL <- cluster::silhouette(res$results$Cls, res$dist_object)
avg_topokmeans_SIL <- mean(topokmeans_SIL[, 3])
avg_topokmeans_SIL


# ************************************* #
# 01/20/2021
# based on Euclidean distance #
covid_us_data <- read.csv("COVID_US_Datasets/covid19us_data_01_20_2021.csv",row.names = 1)
res <- TopoKmeans(covid_us_data,nKNN=48, nClust=4, power=20, sigma =20, 
                  preserveOrdering=FALSE, null_dim = TRUE, first_dim = FALSE) 
topokmeans_SIL <- cluster::silhouette(res$results$Cls, dist(covid_us_data, method = "euclidean"))
avg_topokmeans_SIL <- mean(topokmeans_SIL[, 3])
avg_topokmeans_SIL

# based on Topological distance #
covid_us_data <- read.csv("COVID_US_Datasets/covid19us_data_01_20_2021.csv",row.names = 1)
res <- TopoKmeans(covid_us_data,nKNN=48, nClust=4, power=20, sigma =20, 
                  preserveOrdering=FALSE, null_dim = TRUE, first_dim = FALSE) 
topokmeans_SIL <- cluster::silhouette(res$results$Cls, res$dist_object)
avg_topokmeans_SIL <- mean(topokmeans_SIL[, 3])
avg_topokmeans_SIL


# ************************************* #
# 01/21/2021
# based on Euclidean distance #
covid_us_data <- read.csv("COVID_US_Datasets/covid19us_data_01_21_2021.csv",row.names = 1)
res <- TopoKmeans(covid_us_data,nKNN=43, nClust=4, power=20, sigma =20, 
                  preserveOrdering=FALSE, null_dim = TRUE, first_dim = FALSE) 
topokmeans_SIL <- cluster::silhouette(res$results$Cls, dist(covid_us_data, method = "euclidean"))
avg_topokmeans_SIL <- mean(topokmeans_SIL[, 3])
avg_topokmeans_SIL

# based on Topological distance #
covid_us_data <- read.csv("COVID_US_Datasets/covid19us_data_01_21_2021.csv",row.names = 1)
res <- TopoKmeans(covid_us_data,nKNN=43, nClust=4, power=20, sigma =20, 
                  preserveOrdering=FALSE, null_dim = TRUE, first_dim = FALSE) 
topokmeans_SIL <- cluster::silhouette(res$results$Cls, res$dist_object)
avg_topokmeans_SIL <- mean(topokmeans_SIL[, 3])
avg_topokmeans_SIL


# ************************************* #
# 01/22/2021
# based on Euclidean distance #
covid_us_data <- read.csv("COVID_US_Datasets/covid19us_data_01_22_2021.csv",row.names = 1)
res <- TopoKmeans(covid_us_data,nKNN=41, nClust=4, power=20, sigma =25, 
                  preserveOrdering=FALSE, null_dim = TRUE, first_dim = FALSE) 
topokmeans_SIL <- cluster::silhouette(res$results$Cls, dist(covid_us_data, method = "euclidean"))
avg_topokmeans_SIL <- mean(topokmeans_SIL[, 3])
avg_topokmeans_SIL

# based on Topological distance #
covid_us_data <- read.csv("COVID_US_Datasets/covid19us_data_01_22_2021.csv",row.names = 1)
res <- TopoKmeans(covid_us_data,nKNN=41, nClust=4, power=20, sigma =25, 
                  preserveOrdering=FALSE, null_dim = TRUE, first_dim = FALSE) 
topokmeans_SIL <- cluster::silhouette(res$results$Cls, res$dist_object)
avg_topokmeans_SIL <- mean(topokmeans_SIL[, 3])
avg_topokmeans_SIL


# ************************************* #
# 01/23/2021
# based on Euclidean distance #
covid_us_data <- read.csv("COVID_US_Datasets/covid19us_data_01_23_2021.csv",row.names = 1)
res <- TopoKmeans(covid_us_data,nKNN=38, nClust=4, power=20, sigma =0.5, 
                  preserveOrdering=FALSE, null_dim = TRUE, first_dim = FALSE) 
topokmeans_SIL <- cluster::silhouette(res$results$Cls, dist(covid_us_data, method = "euclidean"))
avg_topokmeans_SIL <- mean(topokmeans_SIL[, 3])
avg_topokmeans_SIL

# based on Topological distance #
covid_us_data <- read.csv("COVID_US_Datasets/covid19us_data_01_23_2021.csv",row.names = 1)
res <- TopoKmeans(covid_us_data,nKNN=38, nClust=4, power=20, sigma =0.5, 
                  preserveOrdering=FALSE, null_dim = TRUE, first_dim = FALSE) 
topokmeans_SIL <- cluster::silhouette(res$results$Cls, res$dist_object)
avg_topokmeans_SIL <- mean(topokmeans_SIL[, 3])
avg_topokmeans_SIL


# ************************************* #
# 01/10/2021 - 01/16/2021, i.e., weekly dataset
# based on Euclidean distance #
covid_us_data <- read.csv("COVID_US_Datasets/avg_final_data_01_10_01_16.csv",row.names = 1)
res <- TopoKmeans(covid_us_data,nKNN=41, nClust=4, power=20, sigma =20, 
                  preserveOrdering=FALSE, null_dim = TRUE, first_dim = FALSE) 
topokmeans_SIL <- cluster::silhouette(res$results$Cls, dist(covid_us_data, method = "euclidean"))
avg_topokmeans_SIL <- mean(topokmeans_SIL[, 3])
avg_topokmeans_SIL

# based on Topological distance #
covid_us_data <- read.csv("COVID_US_Datasets/avg_final_data_01_10_01_16.csv",row.names = 1)
res <- TopoKmeans(covid_us_data,nKNN=41, nClust=4, power=20, sigma =20, 
                  preserveOrdering=FALSE, null_dim = TRUE, first_dim = FALSE) 
topokmeans_SIL <- cluster::silhouette(res$results$Cls, res$dist_object)
avg_topokmeans_SIL <- mean(topokmeans_SIL[, 3])
avg_topokmeans_SIL


# ************************************* #
# 01/17/2021 - 01/23/2021, i.e., weekly dataset
# based on Euclidean distance #
covid_us_data <- read.csv("COVID_US_Datasets/avg_final_data_01_17_01_23.csv",row.names = 1)
res <- TopoKmeans(covid_us_data,nKNN=43, nClust=4, power=20, sigma =20, 
                  preserveOrdering=FALSE, null_dim = TRUE, first_dim = FALSE) 
topokmeans_SIL <- cluster::silhouette(res$results$Cls, dist(covid_us_data, method = "euclidean"))
avg_topokmeans_SIL <- mean(topokmeans_SIL[, 3])
avg_topokmeans_SIL

# based on Topological distance #
covid_us_data <- read.csv("COVID_US_Datasets/avg_final_data_01_10_01_16.csv",row.names = 1)
res <- TopoKmeans(covid_us_data,nKNN=43, nClust=4, power=20, sigma =20, 
                  preserveOrdering=FALSE, null_dim = TRUE, first_dim = FALSE) 
topokmeans_SIL <- cluster::silhouette(res$results$Cls, res$dist_object)
avg_topokmeans_SIL <- mean(topokmeans_SIL[, 3])
avg_topokmeans_SIL
