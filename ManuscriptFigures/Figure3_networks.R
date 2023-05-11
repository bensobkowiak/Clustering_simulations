require(seedy)
require(seqinr)
library(pROC)
require(tidyr)
require(igraph)
require(foreach)
require(doMC)
require(ggplot2)
require(ggpubr)
require(RColorBrewer)
registerDoMC(6)
options(stringsAsFactors = F)
# Run simulations

source("~/R/cov2clusters_simulation.R")
source("~/R/simulate_outbreak.R")

ref<-read.fasta("../MN908947.3.fasta",forceDNAtolower = F) # download Wuhan1 reference strain
ref<-as.character(unlist(ref))
ref[which(ref=="A")]<-1
ref[which(ref=="C")]<-2
ref[which(ref=="G")]<-3
ref[which(ref=="T")]<-4
ref<-as.integer(ref)

sampling.rate.out<-c(0.1,0.25,0.5,1)
inf.rate.out<-c(0.25,0.6)
inf.rate.in<-0.2
removal.rate.all<-0.2
mutation.rate<-0.000003
nums<-0
results<-list()
seq.error<-0.01
threshold<-0.8
plots<-0

Clust_res<-data.frame(TPR=rep(0,20),FPR=rep(0,20))

for (m in 1:2){
  if (m==1){
    init.sus<-500
  } else if (m==2){
    init.sus<-250
  } 
  repeat { 
    W<-simulate_outbreak(init.sus = init.sus,init.sus.in=0,
                         inf.rate=inf.rate.out[m], shape = gamma,
                         samp.freq = 10, mut.rate.out = mutation.rate, 
                         intr.rate = 0,rem.rate=removal.rate.all,
                         min.cases.in=0,min.cases = 100,
                         ref.strain = ref,g.len = 29903)

    size_out<-length(W$epidata[,5])
    sampled.size<-size_out*0.1
    if(sampled.size>=20 & sampled.size<=30) break
  }
  
  for (d in 1:nrow(W$epidata)){
    if (d!=1){
      genome.num<-W$epidata[d,5]
      genome.nuc<-W$nuc[genome.num][[1]]
      change<-which(rbinom(length(genome.nuc),1,seq.error)==1) ## remove SNPs based on seq.error rate
      if (length(change)>0){
        if (length(which(genome.num==W$epidata[,5]))==1){
          W$nuc[genome.num][[1]]<-genome.nuc[-c(change)]
          W$libr[genome.num][[1]]<-W$libr[genome.num][[1]][-c(change)]
        } else { # make new genome if shared by another individual
          newnum<-max(W$epidata[,5])+1
          W$nuc[newnum][[1]]<-genome.nuc[-c(change)]
          W$libr[newnum][[1]]<-W$libr[genome.num][[1]][-c(change)]
          W$epidata[d,5]<-newnum
          W$librstrains<-c(W$librstrains,newnum)
        }
      }
    }
  }
  distmat_all<-gd(W$epidata[,5], W$libr, W$nuc, W$librstrains)
  colnames(distmat_all)<-W$epidata[,1]
  row.names(distmat_all)<-W$epidata[,1]
  for (k in 1:4){
    nums<-nums+1
    sampled.out<-which(rbinom(size_out, 1, sampling.rate.out[k])==1)
    distmat_out <- gd(W$epidata[sampled.out,5], W$libr, W$nuc, W$librstrains)
    colnames(distmat_out)<-W$epidata[sampled.out,1]
    row.names(distmat_out)<-W$epidata[sampled.out,1]
    
    # Make sampling dates 
    dates.out<-data.frame(names=W$epidata[sampled.out,1],dates=W$sampledata[sampled.out,3])
    distmat <- distmat_out
    from <- NULL
    to <- NULL
    gen.dist <- NULL
    time.dist <- NULL
    d1 <- NULL
    index <- 0
    
    for (i in 1:(nrow(distmat)-1)) {
      for (j in (i+1):nrow(distmat)) {
        index <- index + 1
        from[index] <- row.names(distmat)[i]
        to[index] <- colnames(distmat)[j]
        gen.dist[index] <- distmat[i,j]
        time.dist[index] <- abs(dates.out$dates[which(dates.out$names==from[index])]-
                                  dates.out$dates[which(dates.out$names==to[index])])
        inf<-c(which(W$epidata[,1]==from[index]),which(W$epidata[,4]==from[index]))
        rec<-c(which(W$epidata[,1]==to[index]),which(W$epidata[,4]==to[index]))
        d1[index] <- ifelse(test = any(inf%in%rec) , yes = 1, no=0)
      }
    }
    df<- data.frame(from, to, gen.dist, time.dist, d1)
    coef_vec_d1 <- c(3,-0.66 ,-0.075)
    z_d1 <- coef_vec_d1[1] + coef_vec_d1[2]*df$gen.dist + coef_vec_d1[3]*df$time.dist
    df$pred.of.d1 <- 1/(1 + exp(-z_d1))
    
      plots<-plots+1
      res_clust<-cov2clusters_sim(distmat_out,dates.out,probThreshold = threshold[n])
      TP<-rep(0,nrow(df))
      FP<-rep(0,nrow(df))
      for (p in 1:nrow(df)){
        if (df$d1[p]==1 & df$pred.of.d1[p]>=threshold[n]){
          TP[p]<-1
        } else if (df$d1[p]!=1 & df$pred.of.d1[p]>=threshold[n]){
          FP[p]<-1
        }
      }
      TPR<-sum(TP)/sum(df$d1)
      FPR<-sum(FP)/(nrow(df)-sum(df$d1))
      clusters<-data.frame(table(res_clust[,2]))
      if (any(clusters$Var1=="-1")){
        clusters<-clusters[-which(clusters$Var1=="-1"),]
      }
      clusters<-clusters[order(clusters$Freq,decreasing = T),]
      ## setup igraph plots
      colors<-c(rep("grey",length(W$epidata[,1])))
      colors[which(W$epidata[,1] %in% colnames(distmat_out))]<-"lightblue"
      if (nrow(clusters)>0){
        if (nrow(clusters)==1){
          colors[which(W$epidata[,1]%in%res_clust[which(res_clust[,2]==
                                                          as.character(clusters$Var1[1])),1])]<-"red"
        } else if (nrow(clusters)>1){
          colors[which(W$epidata[,1]%in%res_clust[which(res_clust[,2]==
                                                          as.character(clusters$Var1[1])),1])]<-"red"
          color.rand<-brewer.pal(nrow(clusters), "Set3")
          for (col in 2:nrow(clusters)){
            colors[which(W$epidata[,1]%in%res_clust[which(res_clust[,2]==
                                                            as.character(clusters$Var1[col])),1])]<-color.rand[col]
          }
        }
      }
      sizes<-c(rep(2,length(W$epidata[,1])))
      sizes[which(W$epidata[,1] %in% colnames(distmat_out))]<-10
      Clust_res$TPR[plots]<-TPR
      Clust_res$FPR[plots]<-FPR
      results[[plots]]<- graph_from_data_frame(
        W$epidata[-1, c(4,1)],
        directed = T,
        vertices = data.frame(
          name = W$epidata[,1],
          inf.time = W$epidata[,2],
          g.id = W$epidata[,5], # genome IDs
          dist = distmat_all[W$epidata[,5], 1], # distances from genome 1
          color = colors, size=sizes, alpha=0.5,
          label.color = "black")) %>%
        set_graph_attr("layout", layout.reingold.tilford) %>%
        set_graph_attr("par", W$parameters)
  }
}
setwd("~/Documents/SFU/Covid_simulation/Figure3/Test")
saveRDS(results, file="outbreak_plots_figure3.RData")
write.csv(Clust_res,"Figure3A_TPR_FPR.csv",row.names = F)

plots<-readRDS("outbreak_plots_figure3.RData")

nums<-0
for (m in 1:2){
  for (k in 1:4){
    for (n in 1:1){
      nums<-nums+1
      tiff(filename=paste0("R=",inf.rate.out[m]*1/removal.rate.all,"Samp=",
                           sampling.rate.out[k],"Threshold=",threshold[n],
                           ".tiff"), res = 300, units = "cm",width = 20,height = 15)
      plot.igraph(plots[[nums]], 
                  vertex.label = NA, 
                  vertex.label.cex=1,
                  edge.arrow.size=0.1, edge.label.cex=2)

      title(paste0(
        "R  = ", inf.rate.out[m]*1/removal.rate.all ,
        "\nSampling = ", sampling.rate.out[k]*100,"%"),cex.main=1.5)
      dev.off()
    }
  }
}
