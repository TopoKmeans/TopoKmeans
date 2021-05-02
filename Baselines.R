# Topological Clustering using Betti Numbers (TopoCBN)
TopoCBN = function(data,nKNN,filt_len=25,dist_matrix=FALSE){
  
  # extract betti sequence from PD
  extract_betti=function(PD){
    b<-numeric(length = filt_len)
    for (k in 1:filt_len)
      b[k]<-sum((scale_seq[k]>=PD[,1])&(scale_seq[k]<PD[,2]))
    b
  }
  # main body
  
  N <- nrow(data) # number of data points
  if (dist_matrix){
    ind <- t(apply(data,1,order))[,1:nKNN]
    maxscale <- max(t(apply(data,1,sort))[,nKNN])
  } else{
    ind <- cbind(1:N,knn.index(data,k=nKNN-1)) # indices of k nearest neighbors
    maxscale<-max(knn.dist(data,k=nKNN-1)[,nKNN-1]) # maxscale
  }
  
  scale_seq<-seq(0,maxscale,length.out = filt_len) # increasing sequence of scale values
  
  # Compute betti sequences
  print("Extracting Betti sequences (1/4)")
  betti_0<-betti_1<-matrix(nrow = N,ncol = filt_len) ##### using filt_len
  for (i in 1:N){
    
    # construct Rips filtration built on neighborhood around point i
    if (dist_matrix){
      ripsFltr<-ripsFiltration(data[ind[i,],ind[i,]],maxdimension = 0,maxscale = maxscale,dist = 'arbitrary') 
    } else{
      ripsFltr<-ripsFiltration(data[ind[i,],],maxdimension = 0,maxscale = maxscale,dist = 'euclidean') 
    }
    # compute PD of Rips filtration
    PD<-filtrationDiag(filtration = ripsFltr, maxdimension = 1)$diagram
    # replace death=Inf with death=maxscale
    PD[PD[,3]==Inf,3]=maxscale
    
    # extract betti-0 and betti-1 sequences from PD
    betti_0[i,]=extract_betti(PD[PD[,1]==0,2:3,drop=F])
    betti_1[i,]=extract_betti(PD[PD[,1]==1,2:3,drop=F])
    
  }
  
  ##################################################
  # Computing relative change in betti sequences
  ##################################################
  print("Computing relative change in Betti sequences (2/4)")
  delta_betti_0<-delta_betti_1<-matrix(nrow = N,ncol = nKNN)
  for (i in 1:N){
    
    bi_bj<-as.matrix(dist(betti_0[ind[i,],]))[1,]
    denom<-norm(as.matrix(betti_0[i,]),type = 'f')
    delta_betti_0[i,]=bi_bj/denom
    
    bi_bj<-as.matrix(dist(betti_1[ind[i,],]))[1,]
    denom<-norm(as.matrix(betti_1[i,]),type = 'f')
    delta_betti_1[i,]=bi_bj/denom 
  }
  
  ##################################################
  # computing cutoff thresholds
  ##################################################
  
  bp_0<-boxplot(as.vector(delta_betti_0),plot = F)
  bp_1<-boxplot(as.vector(delta_betti_1),plot = F)
  
  cutoff_0<-bp_0$stats[5,] # cutoff threshold for betti-0. Can also be set manually
  cutoff_1<-bp_1$stats[5,] # cutoff threshold for betti-1. Can also be set manually
  
  ############################
  # Forming adjacency matrix
  ############################
  print('Forming adjacency matrix (3/4)')
  A<-matrix(0,ncol = N,nrow = N) 
  for (i in 1:N)
  {
    index_0<-ind[i,which(delta_betti_0[i,]<=cutoff_0)]
    index_1<-ind[i,which(delta_betti_1[i,]<=cutoff_1)]
    
    A[i,intersect(index_0,index_1)]=1
    
  }
  ############################
  # Performing clustering
  ############################
  print('Performing clustering (4/4)')
  g<-graph_from_adjacency_matrix(A,mode='directed') # form graph g from adj. matrix A
  clstrs<-clusters(g,mode = 'strong') # define clusters as strongly connected compents of graph g
  
  return(list(assignments=clstrs$membership,nClust=clstrs$no,cSize=clstrs$csize))
}

# TopoCBN
covid_us_data <- read.csv("COVID_US_Datasets/covid19us_data_01_11_2021.csv",row.names = 1)
dMatrix_0 <- as.matrix(dist(covid_us_data, method = "euclidean"))
# for TopoCBN, user needs to tune nKNN to obtain the expected number of clusters
result <- TopoCBN(dMatrix_0,nKNN=47,dist_matrix = TRUE)
CBN_SIL<- cluster::silhouette(result$assignments, dist(avg_final_data, method = "euclidean"))
avg_CBN_SIL<- mean(CBN_SIL[, 3])
avg_CBN_SIL
result$assignments

# Kmeans over vectorized PDs
# Vectorized PDs generation
VectorizedPDs = function(data, nKNN, dist_matrix=FALSE, preserveOrdering=FALSE, null_dim = FALSE, first_dim= FALSE){
  
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
    dimM[i] = nrow(PD) 
    
    # pot 
    # Pdata0 = PD[PD[,1]==0,2:3]
    # Pdata1 = PD[PD[,1]==1,2:3]
    # pot 2#
    dim0[i] = length(which(PD[,1]==0))
    dim1[i] = length(which(PD[,1]==1))
    #------#
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
  if(null_dim){ # null_dim is True, i.e., only consider 0-dimensional features
    persistence = perst.val_0
  }else if (first_dim){
    persistence = perst.val_1
  }else{
    persistence = perst.val
  }
  return(persistence)
}

# Kmeans over vectorized PDs 
# 01/10/2021
covid_us_data <- read.csv("COVID_US_Datasets/covid19us_data_01_10_2021.csv",row.names = 1)
vectorized_PDs = t(VectorizedPDs(covid_us_data,nKNN=42, preserveOrdering=FALSE, null_dim = TRUE, first_dim = FALSE))
kmeans_PDs_res = kmeansClustering(vectorized_PDs, ClusterNo = 4, PlotIt = FALSE, Verbose = F,Type = 'Steinley')
dMatrix_kmeans_PDs = as.matrix(dist(covid_us_data,method = "euclidean"))
regular_kmeans_PDs_SIL <- cluster::silhouette(kmeans_PDs_res$Cls, dMatrix_kmeans_PDs)
avg_kmeans_PDs_SIL <- mean(regular_kmeans_PDs_SIL[, 3])
avg_kmeans_PDs_SIL


# Kmeans 
covid_us_data <- read.csv("COVID_US_Datasets/covid19us_data_01_11_2021.csv",row.names = 1)
dMatrix_1 <- as.matrix(dist(covid_us_data,method = "euclidean"))
kmeans_res <- kmeansClustering(dMatrix_1, ClusterNo=4,PlotIt=FALSE,Verbose = F)
regular_SIL <- cluster::silhouette(kmeans_res$Cls, dMatrix_1)
avg_kmeans_SIL <- mean(regular_SIL[, 3])
avg_kmeans_SIL


# Agglomerative Hierarchical Clustering 
covid_us_data <- read.csv("COVID_US_Datasets/covid19us_data_01_11_2021.csv",row.names = 1)
hc <- hclust(dist(covid_us_data, method = "euclidean"), "complete")
hc_clusters <- cutree(hc, k = 4) 
hc_SIL <- cluster::silhouette(hc_clusters, dist(covid_us_data, method = "euclidean"))
avg_hc_SIL <- mean(hc_SIL[, 3])
avg_hc_SIL
