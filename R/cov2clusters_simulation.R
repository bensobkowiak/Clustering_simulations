
#### SARS-CoV-2 genomic clusters from phylogenetic trees 
# This script contains the functions required to cluster sequences from simulated outbreaks
# using cov2clusters in Sobkowiak et.al. 2022 https://doi.org/10.1186/s12864-022-08936-4

# Logit probability function
TransProbs<-function(dataInput,probThreshold,dateThreshold,beta){
  PatDist<-dataInput$PatDist #Patristic distance matrix
  dates<-as.numeric(dataInput$dates[,2]) #Date vector
  
  transmat<-foreach(i=seq(ncol(PatDist)-1), .combine = rbind) %dopar% {
    transmatrix<-matrix(ncol = 3)
    for (j in seq(i+1, ncol(PatDist))){
      covariates<-c(1,0,0) #Empty covariates
      covariates[2]<-PatDist[i,j]
      covariates[3]<-abs(as.numeric(dates[i])-as.numeric(dates[j]))
      beta<-beta[1:3]
      
      # Run logit if date threshold met
      if (!is.na(dateThreshold) && covariates[3]>(dateThreshold)){
        trans<-0
      } else {
        trans = 1/(1 + exp(- (sum(beta*covariates))))
      }
      
      # Append to transmatrix if > min probThreshold
      if (trans>=min(probThreshold)){
        transres<-c(colnames(PatDist)[i],colnames(PatDist)[j],trans)
        transmatrix<-rbind(transmatrix,transres)
      }
    }
    transmatrix
  }
  colnames(transmat)<-c("host1","host2","Prob")
  transmat<-na.omit(transmat)
}

# Number clusters
numberClusters<-function(acceptTrans){
  names<-unique(c(as.character(acceptTrans[,1]),as.character(acceptTrans[,2]))) # names of sequences
  cluster_results<-matrix(NA,ncol = 2, nrow = length(names)) # Empty results matrix
  cluster_results[,1]<-names
  if (nrow(cluster_results)>0){
    nums<-1:nrow(cluster_results) # set up list of putative cluster numbers
    for (i in 1:nrow(cluster_results)){
      if (i%in%nums){
        # Find the sequences that cluster with sequence i
        clustnums<-which(cluster_results[,1] %in% unique(unlist(c(acceptTrans[which(acceptTrans[,1]==cluster_results[i,1] | 
                                                                                      acceptTrans[,2]==cluster_results[i,1]),1:2])))) 
        # Find the sequences that cluster with any sequences found to cluster with sequence i
        assocnums<-which(cluster_results[,1] %in% unique(unlist(c(acceptTrans[which(acceptTrans[,1] %in% cluster_results[clustnums,1] | 
                                                                                      acceptTrans[,2] %in% cluster_results[clustnums,1]),1:2]))))
        allnums<-c(clustnums,assocnums)
        prevClust<-unique(cluster_results[allnums,2]) # find all cluster numbers of sequences clustering in this loop
        prevClust<-as.numeric(prevClust[!is.na(prevClust)])
        if (length(prevClust)>0){ # If sequences have been previously clustered, re-assign all sequences to same cluster with the minimum cluster number
          cluster_results[unique(c(which(cluster_results[,2]%in%prevClust),allnums)),2]<-min(prevClust)
        } else{ # Otherwise give cluster new number
          cluster_results[allnums,2]<-min(nums)
        }
        nums<-nums[!nums%in%clustnums] # Remove cluster number used from putative cluster numbers
      }
    }
  }
  return(cluster_results)
}


## Main clustering function

cov2clusters_sim<-function(distmat=distmat,dates=dates,
                           beta=c(3,-0.66,-0.075),
                           dateThreshold=40,
                           probThreshold=0.9,
                           no.Cores=8){
  # Load packages
  require(reshape2)
  require(stringi)
  require(foreach)
  require(doMC)
  registerDoMC(no.Cores)
  options(stringsAsFactors = F)
  
  #Data Setup
  dataInput<-list()
  
  # inputs
  dataInput$PatDist<-as.matrix(distmat)
  dataInput$dates<-dates
  
  # Transmission matrix
  transmat<-TransProbs(dataInput,probThreshold,dateThreshold,beta)
  acceptTrans<-transmat[transmat[,3]>=probThreshold,,drop=FALSE]
  
  # Number clusters
  cluster_results<-numberClusters(acceptTrans)
  cluster_results<-data.frame(cluster_results)
  
  # Output file 
  cluster_results<-cluster_results[order(cluster_results[,2]),]
  row.names(cluster_results)<-NULL
  colnames(cluster_results)<-c("SampleID","Cluster.No")
  noclust<-cbind(as.character(dates[which(!as.character(dates[,1]) %in% cluster_results[,1]),1]),"-1")
  
  if (ncol(noclust)==2){
    cluster_results<-rbind(as.matrix(cluster_results),noclust)
  }
  return(cluster_results)
}

cluster_by_SNPs_sim<-function(distmat, threshold){ 
  options(stringsAsFactors = F)
  distMatrix<-distmat # Matrix of pairwise distances
  distMatrix[lower.tri(distMatrix,diag = T)] = NA 
  distmat<-reshape2::melt(as.matrix(distMatrix), varnames = c('row', 'col'), na.rm = TRUE) # All pairwise distances
  row.names(distmat)<-NULL
  threshold<-as.numeric(threshold) # Clustering treshold
  close_relate<-distmat[which(as.numeric(distmat[,3])<=threshold),] # Find pairwise distances lower than clustering threshold
  names<-unique(c(as.character(close_relate[,1]),as.character(close_relate[,2]))) # All sequence names
  cluster_results<-as.data.frame(matrix(0,ncol = 2, nrow = length(names))) # Empty results matrix
  if (length(names)>0){
  cluster_results[,1]<-as.character(names)
  cluster_results[,2]<-1:length(names) # Cluster numbers
  for (i in 1:nrow(cluster_results)){ # loop through all sequences 
    newclose<-close_relate[which(close_relate$row==cluster_results[i,1] | close_relate$col==cluster_results[i,1]),] # pairwise distances including sequence i
    cluster_results[which(cluster_results[,1] %in% 
                            unique(c(as.character(newclose[,1]),as.character(newclose[,2])))),2]<-cluster_results[i,2] # Find all sequences that cluster with sequence i and assign them the same cluster
    close_relate<-close_relate[-which(close_relate$row==cluster_results[i,1] | close_relate$col==cluster_results[i,1]),] #remove these pairwise distances from the list
  }
  row.names(cluster_results)<-NULL
  colnames(cluster_results)<-c("SampleID",paste0("SNPs",threshold))
  cluster_results<-cluster_results[order(cluster_results[,2]),]
  }
  cluster_results[which(is.na(cluster_results[,2])),2]<-"-1" # Assign unclustered sequences as "-1"
  return(cluster_results)
}
