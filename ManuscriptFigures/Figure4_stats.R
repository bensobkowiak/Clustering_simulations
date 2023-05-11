require(seedy)
require(seqinr)
library(pROC)
require(tidyr)
require(igraph)
require(ggplot2)
require(ggpubr)
require(RColorBrewer)
require(doMC)
registerDoMC(cores=4)
options(stringsAsFactors = F)

source("~/R/simulate_outbreak.R")
source("~/R/cov2clusters_simulation.R")

# Run simulations
ref<-read.fasta("../MN908947.3.fasta",forceDNAtolower = F)
ref<-as.character(unlist(ref))
ref[which(ref=="A")]<-1
ref[which(ref=="C")]<-2
ref[which(ref=="G")]<-3
ref[which(ref=="T")]<-4
ref<-as.integer(ref)

# set up simulation conditions
sampling.rate.in<-c(0.25,1)
inf.rate.out<-0.4
inf.rate.in<-c(0.25,0.6)
removal.rate.all<-0.2
mutation.rate.start<-0.000003
introduction.rate<-c(0.5,0.1)
results<-list()
nums<-0
seq.error<-0.01
time.lag<-NULL
min.perc<-20
init.infection<-20
mindiversity<-c(0,15)
maxdiversity<-c(0,25)
init.sus.out<-300
nreps<-100

clustering_stats<-matrix(0,ncol = 20)

for (m in 1:2){
  for (p in 1:2){
    for (k in 1:2){
      
      
      if (m==1){
        init.sus.in<-550
      } else if (m==2){
        init.sus.in<-250
      }
      cluster_reps_stats<-foreach(i=1:nreps, .combine = rbind) %dopar% {
        cluststat<-matrix(0,ncol = 20,nrow=2)
        repeat {
          W<-simulate_outbreak(init.sus = init.sus.out,init.sus.in=init.sus.in,
                               inf.rate=inf.rate.out, inf.rate.in = inf.rate.in[m],
                               shape = gamma,
                               samp.freq = 10, time.lag.outside = time.lag,
                               min.perc.outside.inf = min.perc,
                               mut.rate.out = mutation.rate.start,
                               mut.rate.in = mutation.rate.start,
                               intr.rate = introduction.rate[k],
                               shape.infect = 6,rate.infect = 3,
                               init.inf=init.infection,
                               min.init.dist = mindiversity[p],
                               max.init.dist = maxdiversity[p],
                               rem.rate=removal.rate.all,min.cases.in=100,
                               min.cases = 100,
                               ref.strain = ref,g.len = 29903)
          
          size_out<-length(W$epidata[,5])
          size_in<-length(W$epidata_in[,6])
          sampled.size<-size_in*0.1
          if (sampled.size>=20 & sampled.size<=30 & size_out>=150 & size_out<=250)
            break
        }
        for (d in 1:nrow(W$epidata_in)){
          if (d!=1){
            genome.num<-W$epidata_in[d,6]
            genome.nuc<-W$nuc[genome.num][[1]]
            change<-which(rbinom(length(genome.nuc),1,seq.error)==1) ## remove SNPs based on seq.error rate
            if (length(change)>0){
              if (length(which(genome.num==W$epidata_in[,6]))==1){
                W$nuc[genome.num][[1]]<-genome.nuc[-c(change)]
                W$libr[genome.num][[1]]<-W$libr[genome.num][[1]][-c(change)]
              } else { # make new genome if shared by another individual
                newnum<-max(c(W$epidata[,5],W$epidata_in[,6]))+1
                W$nuc[newnum][[1]]<-genome.nuc[-c(change)]
                W$libr[newnum][[1]]<-W$libr[genome.num][[1]][-c(change)]
                W$epidata_in[d,6]<-newnum
                W$librstrains<-c(W$librstrains,newnum)
              }
            }
          }
        }
        distmat_all<-gd(c(W$epidata[,5],W$epidata_in[,6]), W$libr, W$nuc, W$librstrains)
        colnames(distmat_all)<-c(W$epidata[,1],W$epidata_in[,1])
        row.names(distmat_all)<-c(W$epidata[,1],W$epidata_in[,1])
        distmat_in <- gd(W$epidata_in[,6], W$libr, W$nuc, W$librstrains)
        colnames(distmat_in)<-W$epidata_in[,1]
        row.names(distmat_in)<-W$epidata_in[,1]
        
        # Make sampling dates from epidata in
        dates.in<-data.frame(names=colnames(distmat_in),dates=0,intro=0,outbreak.in=0)
        for (i in 1:nrow(dates.in)){ # pick random sampling rate between infection and removal time
          dates.in$dates[i]<-sample(W$epidata_in[i,2]:W$epidata_in[i,3],1) 
        }
        
        # Assign all sampled and unsampled cases inside to an outbreak
        introduced_cases<-unique(W$epidata_in[which(W$epidata_in[,5]==1),4])
        no.introductions<-length(introduced_cases)
        dates.in$intro[which(W$epidata_in[,4] %in% introduced_cases)]<-1
        for (i in 1:no.introductions){
          intro<-which(W$epidata_in[,4]==introduced_cases[i])
          dates.in$outbreak.in[intro]<-i
          repeat {
            cluster.ids<-W$epidata_in[which(dates.in$outbreak.in==i),1]
            new.seq<-which(W$epidata_in[,4] %in% cluster.ids & dates.in$outbreak.in==0)
            dates.in$outbreak.in[new.seq]<-i
            if(length(new.seq)==0) break
          }
        }
        
        nums<-0
        for (j in 1:2){
          nums<-nums+1
          sampled.in<-which(rbinom(size_in, 1, sampling.rate.in[j])==1)
          distmat_in.sampled<-distmat_in[sampled.in,sampled.in]
          # Run cov2clusters
          dates.in.sampled<-dates.in[sampled.in,]
          results_sim.in<-cov2clusters_sim(distmat = distmat_in.sampled,dates = dates.in.sampled,probThreshold = 0.8) # run logit model
          
          # size of true inside outbreaks
          inside.outbreaks<-data.frame(table(dates.in.sampled$outbreak.in))       
          
          # Sampled outbreaks stats
          cluststat[nums,1]<-inf.rate.in[m]*1/removal.rate.all
          cluststat[nums,2]<-inf.rate.out*1/removal.rate.all
          cluststat[nums,3]<-sampling.rate.in[j]
          cluststat[nums,4]<-introduction.rate[k]
          cluststat[nums,5]<-no.introductions
          cluststat[nums,6]<-median(distmat_in.sampled)
          # True sampled inside outbreaks stats
          cluststat[nums,7]<-length(which(inside.outbreaks$Freq==1))
          cluststat[nums,8]<-length(which(inside.outbreaks$Freq!=1))
          cluststat[nums,9]<-max(inside.outbreaks$Freq)
          cluststat[nums,10]<-min(inside.outbreaks$Freq[which(inside.outbreaks$Freq!=1)])
          cluststat[nums,11]<-mean(inside.outbreaks$Freq[which(inside.outbreaks$Freq!=1)])
          cluststat[nums,12]<-sd(inside.outbreaks$Freq[which(inside.outbreaks$Freq!=1)])
          
          # clustering stats
          cluster_sizes<-data.frame(table(results_sim.in[,2]))
          cluster_sizes<-cluster_sizes[order(cluster_sizes$Freq,decreasing = T),]
          if ("-1" %in% cluster_sizes$Var1){
            cluststat[nums,13]<-cluster_sizes$Freq[which(cluster_sizes$Var1=="-1")]/length(sampled.in)
            cluster_sizes<-cluster_sizes[-which(cluster_sizes$Var1=="-1"),]
          }
          if (nrow(cluster_sizes)>0){
            cluststat[nums,14]<-cluster_sizes$Freq[1]/length(sampled.in)
            cluststat[nums,15]<-nrow(cluster_sizes)
            cluststat[nums,16]<-max(cluster_sizes$Freq)
            cluststat[nums,17]<-min(cluster_sizes$Freq)
            cluststat[nums,18]<-mean(cluster_sizes$Freq)
            cluststat[nums,19]<-sd(cluster_sizes$Freq)
          }
          cluststat[nums,20]<-maxdiversity[p]
        }
        
        cluststat
      }
      clustering_stats<-rbind(clustering_stats,cluster_reps_stats) 
      
    }
  }
}



clustering_stats<-clustering_stats[-1,]
colnames(clustering_stats)<-c("R_in", "R_out","sampling.prop.in","intro.rate",
                              "no.introductions","Median.SNP.distance.inside",
                              "no.singletons","no.inside.outbreaks",
                              "max.inside.outbreak",
                              "min.inside.outbreak","mean.inside.outbreak","sd.inside.outbreak",
                              "prop.unclustered",
                              "prop.largest.cluster","number.clusters","max.clust.size",
                              "min.clust.size","mean.clust.size","sd.clust.size",
                              "Max.initial.diversity")

write.csv(clustering_stats,"Figure4_stats.csv",row.names = F)

clustering_stats_reps<-read.csv("Figure4_stats.csv")
clustering_stats_reps$clusters_outbreaks<-clustering_stats_reps$number.clusters/clustering_stats_reps$no.inside.outbreaks

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



### Figure 4C
clustering_stats_intro<-clustering_stats_reps[which(clustering_stats_reps$intro.rate==0.5),]

df_prop.non.clustered <- data_summary(clustering_stats_intro, varname="prop.unclustered", 
                                      groupnames=c("R_in","Max.initial.diversity",
                                                   "sampling.prop.in"))
df_prop.non.clustered$sampling.prop.in=factor(df_prop.non.clustered$sampling.prop.in,levels=c("0.25","1"))
df_prop.non.clustered$R_in=factor(df_prop.non.clustered$R_in,levels=c("1.25","3"))
df_prop.non.clustered$Max.initial.diversity<-factor(df_prop.non.clustered$Max.initial.diversity,levels=c("0","25"))

lowerror1<-df_prop.non.clustered$prop.unclustered-df_prop.non.clustered$sd
higherror1<-df_prop.non.clustered$prop.unclustered+df_prop.non.clustered$sd
lowerror1[which(lowerror1<0)]<-0

p1<-ggplot(df_prop.non.clustered, aes(x=sampling.prop.in, 
                                      y=prop.unclustered,color=Max.initial.diversity, fill=R_in)) + 
  geom_col(position = position_dodge(0.95),size=1) + 
  theme_bw() + 
  ylab("Proportion of unclustered sequences") +ylim(c(0,0.2)) +
  xlab("Sampling proportion") +
  scale_color_manual(values=c("lightgrey", "black"),name="Max initial diversity") +
  scale_fill_manual(values=c("lightblue", "royalblue"),name="R")  +
  geom_errorbar(aes(ymin=lowerror1, ymax=higherror1), width=.2,
                position=position_dodge(0.95))  + 
  guides(color=guide_legend(override.aes=list(fill=NA,size=3),title = "Max initial \ndiversity"))
plot(p1)

# number clusters

df_prop.non.clustered <- data_summary(clustering_stats_intro, varname="prop.largest.cluster", 
                                      groupnames=c("R_in","Max.initial.diversity",
                                                   "sampling.prop.in"))
df_prop.non.clustered$sampling.prop.in=factor(df_prop.non.clustered$sampling.prop.in,levels=c("0.25","1"))
df_prop.non.clustered$R_in=factor(df_prop.non.clustered$R_in,levels=c("1.25","3"))
df_prop.non.clustered$Max.initial.diversity<-factor(df_prop.non.clustered$Max.initial.diversity,levels=c("0","25"))

lowerror2<-df_prop.non.clustered$prop.largest.cluster-df_prop.non.clustered$sd
higherror2<-df_prop.non.clustered$prop.largest.cluster+df_prop.non.clustered$sd
lowerror2[which(lowerror2<0)]<-0
higherror2[which(higherror2>1)]<-1
p2<-ggplot(df_prop.non.clustered, aes(x=sampling.prop.in, 
                                      y=prop.largest.cluster,color=Max.initial.diversity, fill=R_in)) + 
  geom_col(position = position_dodge(0.95),size=1) + 
  theme_bw() + 
  ylab("Proportion of sequences in largest cluster") +ylim(c(0,1)) +
  xlab("Sampling proportion") +
  scale_color_manual(values=c("lightgrey", "black"),name="Max initial diversity") +
  scale_fill_manual(values=c("lightblue", "royalblue"),name="R")  +
  geom_errorbar(aes(ymin=lowerror2, ymax=higherror2), width=.2,
                position=position_dodge(0.95))  + 
  guides(color=guide_legend(override.aes=list(fill=NA,size=3),title = "Max initial \ndiversity"))
plot(p2)

## number clusters per outbreak

df_prop.non.clustered <- data_summary(clustering_stats_intro, varname="clusters_outbreaks", 
                                      groupnames=c("R_in","Max.initial.diversity",
                                                   "sampling.prop.in"))
df_prop.non.clustered$sampling.prop.in=factor(df_prop.non.clustered$sampling.prop.in,levels=c("0.25","1"))
df_prop.non.clustered$R_in=factor(df_prop.non.clustered$R_in,levels=c("1.25","3"))
df_prop.non.clustered$Max.initial.diversity<-factor(df_prop.non.clustered$Max.initial.diversity,levels=c("0","25"))

lowerror3<-df_prop.non.clustered$clusters_outbreaks-df_prop.non.clustered$sd
higherror3<-df_prop.non.clustered$clusters_outbreaks+df_prop.non.clustered$sd


p3<-ggplot(df_prop.non.clustered, aes(x=sampling.prop.in, 
                                      y=clusters_outbreaks,color=Max.initial.diversity, fill=R_in)) + 
  geom_col(position = position_dodge(0.95),size=1) + 
  theme_bw() + 
  ylab("Number of genomic clusters per introduced outbreak") +ylim(c(0,3)) +
  xlab("Sampling proportion") +
  scale_color_manual(values=c("lightgrey", "black"),name="Max initial diversity") +
  scale_fill_manual(values=c("lightblue", "royalblue"),name="R")  +
  geom_errorbar(aes(ymin=lowerror3, ymax=higherror3), width=.2,
                position=position_dodge(0.95))  + 
  guides(color=guide_legend(override.aes=list(fill=NA,size=3),title = "Max initial \ndiversity"))
plot(p3)


tiff(filename="Figure4C.tiff", res = 300, units = "cm",width = 30,height = 13)
p5<-ggarrange(p1,p2,p3,
              ncol=3, nrow=1, common.legend = TRUE, legend="right")
plot(p5)
dev.off()


### Supplementaty figure 2C
clustering_stats_intro<-clustering_stats_reps[which(clustering_stats_reps$intro.rate==0.1),]

df_prop.non.clustered <- data_summary(clustering_stats_intro, varname="prop.unclustered", 
                                      groupnames=c("R_in","Max.initial.diversity",
                                                   "sampling.prop.in"))
df_prop.non.clustered$sampling.prop.in=factor(df_prop.non.clustered$sampling.prop.in,levels=c("0.25","1"))
df_prop.non.clustered$R_in=factor(df_prop.non.clustered$R_in,levels=c("1.25","3"))
df_prop.non.clustered$Max.initial.diversity<-factor(df_prop.non.clustered$Max.initial.diversity,levels=c("0","25"))

lowerror1<-df_prop.non.clustered$prop.unclustered-df_prop.non.clustered$sd
higherror1<-df_prop.non.clustered$prop.unclustered+df_prop.non.clustered$sd
lowerror1[which(lowerror1<0)]<-0

p1<-ggplot(df_prop.non.clustered, aes(x=sampling.prop.in, 
                                      y=prop.unclustered,color=Max.initial.diversity, fill=R_in)) + 
  geom_col(position = position_dodge(0.95),size=1) + 
  theme_bw() + 
  ylab("Proportion of unclustered sequences") +ylim(c(0,0.2)) +
  xlab("Sampling proportion") +
  scale_color_manual(values=c("lightgrey", "black"),name="Max initial diversity") +
  scale_fill_manual(values=c("lightblue", "royalblue"),name="R")  +
  geom_errorbar(aes(ymin=lowerror1, ymax=higherror1), width=.2,
                position=position_dodge(0.95))  + 
  guides(color=guide_legend(override.aes=list(fill=NA,size=3),title = "Max initial \ndiversity"))
plot(p1)

# number clusters

df_prop.non.clustered <- data_summary(clustering_stats_intro, varname="prop.largest.cluster", 
                                      groupnames=c("R_in","Max.initial.diversity",
                                                   "sampling.prop.in"))
df_prop.non.clustered$sampling.prop.in=factor(df_prop.non.clustered$sampling.prop.in,levels=c("0.25","1"))
df_prop.non.clustered$R_in=factor(df_prop.non.clustered$R_in,levels=c("1.25","3"))
df_prop.non.clustered$Max.initial.diversity<-factor(df_prop.non.clustered$Max.initial.diversity,levels=c("0","25"))

lowerror2<-df_prop.non.clustered$prop.largest.cluster-df_prop.non.clustered$sd
higherror2<-df_prop.non.clustered$prop.largest.cluster+df_prop.non.clustered$sd
lowerror2[which(lowerror2<0)]<-0
higherror2[which(higherror2>1)]<-1
p2<-ggplot(df_prop.non.clustered, aes(x=sampling.prop.in, 
                                      y=prop.largest.cluster,color=Max.initial.diversity, fill=R_in)) + 
  geom_col(position = position_dodge(0.95),size=1) + 
  theme_bw() + 
  ylab("Proportion of sequences in largest cluster") +ylim(c(0,1)) +
  xlab("Sampling proportion") +
  scale_color_manual(values=c("lightgrey", "black"),name="Max initial diversity") +
  scale_fill_manual(values=c("lightblue", "royalblue"),name="R")  +
  geom_errorbar(aes(ymin=lowerror2, ymax=higherror2), width=.2,
                position=position_dodge(0.95))  + 
  guides(color=guide_legend(override.aes=list(fill=NA,size=3),title = "Max initial \ndiversity"))
plot(p2)

## number clusters per outbreak

df_prop.non.clustered <- data_summary(clustering_stats_intro, varname="clusters_outbreaks", 
                                      groupnames=c("R_in","Max.initial.diversity",
                                                   "sampling.prop.in"))
df_prop.non.clustered$sampling.prop.in=factor(df_prop.non.clustered$sampling.prop.in,levels=c("0.25","1"))
df_prop.non.clustered$R_in=factor(df_prop.non.clustered$R_in,levels=c("1.25","3"))
df_prop.non.clustered$Max.initial.diversity<-factor(df_prop.non.clustered$Max.initial.diversity,levels=c("0","25"))

lowerror3<-df_prop.non.clustered$clusters_outbreaks-df_prop.non.clustered$sd
higherror3<-df_prop.non.clustered$clusters_outbreaks+df_prop.non.clustered$sd


p3<-ggplot(df_prop.non.clustered, aes(x=sampling.prop.in, 
                                      y=clusters_outbreaks,color=Max.initial.diversity, fill=R_in)) + 
  geom_col(position = position_dodge(0.95),size=1) + 
  theme_bw() + 
  ylab("Number of genomic clusters per introduced outbreak") +ylim(c(0,3)) +
  xlab("Sampling proportion") +
  scale_color_manual(values=c("lightgrey", "black"),name="Max initial diversity") +
  scale_fill_manual(values=c("lightblue", "royalblue"),name="R")  +
  geom_errorbar(aes(ymin=lowerror3, ymax=higherror3), width=.2,
                position=position_dodge(0.95))  + 
  guides(color=guide_legend(override.aes=list(fill=NA,size=3),title = "Max initial \ndiversity"))
plot(p3)


tiff(filename="Supplementary_figure2C.tiff", res = 300, units = "cm",width = 30,height = 13)
p5<-ggarrange(p1,p2,p3,
              ncol=3, nrow=1, common.legend = TRUE, legend="right")
plot(p5)
dev.off()
