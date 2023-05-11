require(seedy)
require(seqinr)
require(pROC)
require(tidyr)
require(igraph)
require(foreach)
require(doMC)
require(ggplot2)
require(ggpubr)
registerDoMC(4)
options(stringsAsFactors = F)
# Run simulations

source("~/R/cov2clusters_simulation.R")
source("~/R/simulate_outbreak.R")

ref<-read.fasta("MN908947.3.fasta",forceDNAtolower = F) # download Wuhan1 reference 
ref<-as.character(unlist(ref))
ref[which(ref=="A")]<-1
ref[which(ref=="C")]<-2
ref[which(ref=="G")]<-3
ref[which(ref=="T")]<-4
ref<-as.integer(ref)

sampling.rate.out<-c(0.1,0.25,0.5,1)
inf.rate.out<-c(0.25,0.6)
removal.rate.all<-0.2
mutation.rate<-0.000003
introduction.rate<-0
nums<-0
results<-list()
seq.error<-0.01
threshold=c(0.5,0.6,0.7,0.8,0.9)

clustering_stats<-matrix(0,ncol = 9)
for (m in 1:2){
  if (m==1){
    init.sus<-500
  } else if (m==2){
    init.sus<-250
  } 
  cluster_reps_stats<-foreach(i=1:100, .combine = rbind) %dopar% {
    cluststat<-matrix(0,ncol = 9,nrow=20)
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
    nums<-0
    for (k in 1:4){
      
      sampled.out<-which(rbinom(size_out, 1, sampling.rate.out[k])==1)
      distmat_out <- gd(W$epidata[sampled.out,5], W$libr, W$nuc, W$librstrains)
      colnames(distmat_out)<-W$epidata[sampled.out,1]
      row.names(distmat_out)<-W$epidata[sampled.out,1]
      
      # Make sampling dates 
      dates.out<-data.frame(names=W$epidata[sampled.out,1],dates=W$sampledata[sampled.out,3])
      #run cov2cluster at different thresholds
      for (p in 1:5){
        nums<-nums+1
        res_clust<-cov2clusters_sim(distmat_out,dates.out,probThreshold = threshold[p])
        cluster_sizes<-data.frame(table(res_clust[,2]))
        cluster_sizes<-cluster_sizes[order(cluster_sizes$Freq,decreasing = T),]
        cluststat[nums,1]<-inf.rate.out[m]*1/removal.rate.all
        cluststat[nums,2]<-sampling.rate.out[k]
        cluststat[nums,3]<-threshold[p]
        if ("-1" %in% cluster_sizes$Var1){
          cluststat[nums,4]<-cluster_sizes$Freq[which(cluster_sizes$Var1=="-1")]/length(sampled.out)
          cluster_sizes<-cluster_sizes[-which(cluster_sizes$Var1=="-1"),]
        }
        if (nrow(cluster_sizes)>0){
          cluststat[nums,5]<-cluster_sizes$Freq[1]/length(sampled.out)
          cluststat[nums,6]<-nrow(cluster_sizes)
          cluststat[nums,7]<-cluster_sizes$Freq[1]
          cluststat[nums,8]<-cluster_sizes$Freq[nrow(cluster_sizes)]
          cluststat[nums,9]<-sd(cluster_sizes$Freq)
        }
      }
    }
    cluststat
  }
  clustering_stats<-rbind(clustering_stats,cluster_reps_stats)
  cat(paste0("\nfinished inf.rate = ",inf.rate.out[m]))
  
}
clustering_stats<-clustering_stats[-1,]
colnames(clustering_stats)<-c("R","sampling.prop","Clustering_threshold","prop.unclustered","prop.largest.cluster","number.clusters","max.clust.size","min.clust.size","sd.clust.size")
write.csv(clustering_stats,"cluster_stats_figure3.csv",row.names = F)


## load results file
clustering_stats_reps<-read.csv("cluster_stats_figure3.csv")
# summarize stats function
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}


### plot 0.8 threshold clustering results

df_prop.non.clustered <- data_summary(clustering_stats_reps, varname="prop.unclustered", 
                                      groupnames=c("R", "sampling.prop","Clustering_threshold"))
df_prop.non.clustered<-df_prop.non.clustered[which(df_prop.non.clustered$Clustering_threshold=="0.8"),]
lowerror1<-df_prop.non.clustered$prop.unclustered-df_prop.non.clustered$sd
higherror1<-df_prop.non.clustered$prop.unclustered+df_prop.non.clustered$sd
lowerror1[which(lowerror1<0)]<-0
df_prop.non.clustered$sampling.prop=as.factor(df_prop.non.clustered$sampling.prop)
df_prop.non.clustered$R=as.factor(df_prop.non.clustered$R)
df_prop.non.clustered$Clustering_threshold<-as.factor(df_prop.non.clustered$Clustering_threshold)
p1<-ggplot(df_prop.non.clustered, aes(x=sampling.prop, y=prop.unclustered, fill=R)) + 
  geom_bar(stat="identity",
           position=position_dodge()) +
  scale_fill_manual(values = c("lightblue","royalblue"),name = "R") +
  geom_errorbar(aes(ymin=lowerror1, ymax=higherror1), width=.2,
                position=position_dodge(.9)) + theme_bw() +
  ylab("Proportion of sequences unclustered") +ylim(c(0,0.3)) +
  xlab("Sampling proportion") 

# prop in largest cluster
df_prop.cluster <- data_summary(clustering_stats_reps, varname="prop.largest.cluster", 
                                groupnames=c("R", "sampling.prop","Clustering_threshold"))
df_prop.cluster<-df_prop.cluster[which(df_prop.cluster$Clustering_threshold=="0.8"),]
lowerror2<-df_prop.cluster$prop.largest.cluster-df_prop.cluster$sd
higherror2<-df_prop.cluster$prop.largest.cluster+df_prop.cluster$sd
higherror2[which(higherror2>1)]<-1
df_prop.cluster$sampling.prop=as.factor(df_prop.cluster$sampling.prop)
df_prop.cluster$R=as.factor(df_prop.cluster$R)
df_prop.cluster$Clustering_threshold<-as.factor(df_prop.cluster$Clustering_threshold)
p2<-ggplot(df_prop.cluster, aes(x=sampling.prop, y=prop.largest.cluster, fill=R)) + 
  geom_bar(stat="identity", 
           position=position_dodge()) +
  scale_fill_manual(values = c("lightblue","royalblue"),name = "R") +
  geom_errorbar(aes(ymin=lowerror2, ymax=higherror2), width=.2,
                position=position_dodge(.9)) + theme_bw() + 
  ylab("Proportion of sequences in largest cluster") +ylim(c(0,1)) +
  xlab("Sampling proportion") 

# Number of clusters
df_clusters <- data_summary(clustering_stats_reps, varname="number.clusters", 
                            groupnames=c("R", "sampling.prop","Clustering_threshold"))
df_clusters<-df_clusters[which(df_clusters$Clustering_threshold=="0.8"),]
df_clusters$sampling.prop=as.factor(df_clusters$sampling.prop)
df_clusters$R=as.factor(df_clusters$R)
p3<-ggplot(df_clusters, aes(x=sampling.prop, y=number.clusters, fill=R)) + 
  geom_bar(stat="identity", 
           position=position_dodge()) +
  scale_fill_manual(values = c("lightblue","royalblue"),name = "R") +
  geom_errorbar(aes(ymin=number.clusters-sd, ymax=number.clusters+sd), width=.2,
                position=position_dodge(.9)) + theme_bw() + 
  scale_y_continuous(breaks=seq(1, 15, 1))+
  ylab("Number of clusters") +
  xlab("Sampling proportion") 

tiff(filename="Figure3A_clusterstats.tiff", res = 300, units = "cm",width = 30,height = 15)
p4<-ggarrange(p1,p2,p3,
              ncol=3, nrow=1, common.legend = TRUE, legend="right")

plot(p4)
dev.off()


### plot different threshold results ####


### proportion unclustered
df_prop.non.clustered <- data_summary(clustering_stats_reps, varname="prop.unclustered", 
                                      groupnames=c("R", "sampling.prop","Clustering_threshold"))
##lowR
df_prop.non.clustered_R<-df_prop.non.clustered[which(df_prop.non.clustered$R=="1.25"),]
lowerror1<-df_prop.non.clustered_R$prop.unclustered-df_prop.non.clustered_R$sd
higherror1<-df_prop.non.clustered_R$prop.unclustered+df_prop.non.clustered_R$sd
lowerror1[which(lowerror1<0)]<-0
df_prop.non.clustered_R$sampling.prop=as.factor(df_prop.non.clustered_R$sampling.prop)
df_prop.non.clustered_R$R=as.factor(df_prop.non.clustered_R$R)
df_prop.non.clustered_R$Clustering_threshold<-as.factor(df_prop.non.clustered_R$Clustering_threshold)
p1<-ggplot(df_prop.non.clustered_R, aes(x=sampling.prop, y=prop.unclustered, fill=Clustering_threshold)) + 
  geom_bar(stat="identity",
           position=position_dodge()) +
  scale_fill_manual(values = c("lightgrey","lightblue","steelblue","royalblue","darkblue"),name = "Logit probability \nthreshold") +
  geom_errorbar(aes(ymin=lowerror1, ymax=higherror1), width=.2,
                position=position_dodge(.9)) + theme_bw() +
  ylab("Proportion of sequences unclustered") +ylim(c(0,0.5)) +
  xlab("Sampling proportion") 


##highR
df_prop.non.clustered_R<-df_prop.non.clustered[which(df_prop.non.clustered$R=="3"),]
lowerror2<-df_prop.non.clustered_R$prop.unclustered-df_prop.non.clustered_R$sd
higherror2<-df_prop.non.clustered_R$prop.unclustered+df_prop.non.clustered_R$sd
lowerror2[which(lowerror2<0)]<-0
df_prop.non.clustered_R$sampling.prop=as.factor(df_prop.non.clustered_R$sampling.prop)
df_prop.non.clustered_R$R=as.factor(df_prop.non.clustered_R$R)
df_prop.non.clustered_R$Clustering_threshold<-as.factor(df_prop.non.clustered_R$Clustering_threshold)
p2<-ggplot(df_prop.non.clustered_R, aes(x=sampling.prop, y=prop.unclustered, fill=Clustering_threshold)) + 
  geom_bar(stat="identity",
           position=position_dodge()) +
  scale_fill_manual(values = c("lightgrey","lightblue","steelblue","royalblue","darkblue"),name = "Logit probability \nthreshold") +
  geom_errorbar(aes(ymin=lowerror2, ymax=higherror2), width=.2,
                position=position_dodge(.9)) + theme_bw() +
  ylab("Proportion of sequences unclustered") +ylim(c(0,0.5)) +
  xlab("Sampling proportion") 

## prop largest cluster
#low R
df_prop.cluster <- data_summary(clustering_stats_reps, varname="prop.largest.cluster", 
                                groupnames=c("R", "sampling.prop","Clustering_threshold"))
df_prop.cluster<-df_prop.cluster[which(df_prop.cluster$R=="1.25"),]
lowerror3<-df_prop.cluster$prop.largest.cluster-df_prop.cluster$sd
higherror3<-df_prop.cluster$prop.largest.cluster+df_prop.cluster$sd
higherror3[which(higherror3>1)]<-1
df_prop.cluster$sampling.prop=as.factor(df_prop.cluster$sampling.prop)
df_prop.cluster$R=as.factor(df_prop.cluster$R)
df_prop.cluster$Clustering_threshold<-as.factor(df_prop.cluster$Clustering_threshold)
p3<-ggplot(df_prop.cluster, aes(x=sampling.prop, y=prop.largest.cluster, fill=Clustering_threshold)) + 
  geom_bar(stat="identity", 
           position=position_dodge()) +
  scale_fill_manual(values = c("lightgrey","lightblue","steelblue","royalblue","darkblue"),name = "Logit probability \nthreshold") +
  geom_errorbar(aes(ymin=lowerror3, ymax=higherror3), width=.2,
                position=position_dodge(.9)) + theme_bw() + 
  ylab("Proportion of sequences in largest cluster") +ylim(c(0,1)) +
  xlab("Sampling proportion") 

#high R
df_prop.cluster <- data_summary(clustering_stats_reps, varname="prop.largest.cluster", 
                                groupnames=c("R", "sampling.prop","Clustering_threshold"))
df_prop.cluster<-df_prop.cluster[which(df_prop.cluster$R=="3"),]
lowerror4<-df_prop.cluster$prop.largest.cluster-df_prop.cluster$sd
higherror4<-df_prop.cluster$prop.largest.cluster+df_prop.cluster$sd
higherror4[which(higherror4>1)]<-1
df_prop.cluster$sampling.prop=as.factor(df_prop.cluster$sampling.prop)
df_prop.cluster$R=as.factor(df_prop.cluster$R)
df_prop.cluster$Clustering_threshold<-as.factor(df_prop.cluster$Clustering_threshold)
p4<-ggplot(df_prop.cluster, aes(x=sampling.prop, y=prop.largest.cluster, fill=Clustering_threshold)) + 
  geom_bar(stat="identity", 
           position=position_dodge()) +
  scale_fill_manual(values = c("lightgrey","lightblue","steelblue","royalblue","darkblue"),name = "Logit probability \nthreshold") +
  geom_errorbar(aes(ymin=lowerror4, ymax=higherror4), width=.2,
                position=position_dodge(.9)) + theme_bw() + 
  ylab("Proportion of sequences in largest cluster") +ylim(c(0,1)) +
  xlab("Sampling proportion") 

## Number of clusters
# lowR
df_clusters <- data_summary(clustering_stats_reps, varname="number.clusters", 
                            groupnames=c("R", "sampling.prop","Clustering_threshold"))
df_clusters<-df_clusters[which(df_clusters$R=="1.25"),]
df_clusters$sampling.prop=as.factor(df_clusters$sampling.prop)
df_clusters$R=as.factor(df_clusters$R)
df_clusters$Clustering_threshold<-as.factor(df_clusters$Clustering_threshold)

p5<-ggplot(df_clusters, aes(x=sampling.prop, y=number.clusters, fill=Clustering_threshold)) + 
  geom_bar(stat="identity", 
           position=position_dodge()) +
  scale_fill_manual(values = c("lightgrey","lightblue","steelblue","royalblue","darkblue"),name = "Logit probability \nthreshold") +
  geom_errorbar(aes(ymin=number.clusters-sd, ymax=number.clusters+sd), width=.2,
                position=position_dodge(.9)) + theme_bw() + 
  scale_y_continuous(breaks=seq(1, 15, 1))+ylim(c(0,10)) +
  ylab("Number of clusters") +
  xlab("Sampling proportion") 

#high R
df_clusters <- data_summary(clustering_stats_reps, varname="number.clusters", 
                            groupnames=c("R", "sampling.prop","Clustering_threshold"))
df_clusters<-df_clusters[which(df_clusters$R=="3"),]
df_clusters$sampling.prop=as.factor(df_clusters$sampling.prop)
df_clusters$R=as.factor(df_clusters$R)
df_clusters$Clustering_threshold<-as.factor(df_clusters$Clustering_threshold)
p6<-ggplot(df_clusters, aes(x=sampling.prop, y=number.clusters, fill=Clustering_threshold)) + 
  geom_bar(stat="identity", 
           position=position_dodge()) +
  scale_fill_manual(values = c("lightgrey","lightblue","steelblue","royalblue","darkblue"),name = "Logit probability \nthreshold") +
  geom_errorbar(aes(ymin=number.clusters-sd, ymax=number.clusters+sd), width=.2,
                position=position_dodge(.9)) + theme_bw() + 
  scale_y_continuous(breaks=seq(1, 10, 1))+ ylim(c(0,10)) +
  ylab("Number of clusters") +
  xlab("Sampling proportion") 


tiff(filename="Supplementary_threshold_singleOutbreak.tiff", res = 300, units = "cm",width = 25,height = 25)
p7<-ggarrange(p1 + xlab(NULL),
              p2 + ylab(NULL) + xlab(NULL),
              p3 + xlab(NULL),
              p4 + ylab(NULL) + xlab(NULL),
              p5,
              p6 + ylab(NULL),
              ncol=2, nrow=3, common.legend = TRUE, legend="right")

plot(p7)
dev.off()
