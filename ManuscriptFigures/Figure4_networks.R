require(seedy)
require(seqinr)
library(pROC)
require(tidyr)
require(igraph)
require(RColorBrewer)
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
inf.rate.out<-c(0.4,0.4)
inf.rate.in<-c(0.25,0.6)
removal.rate.all<-0.2
mutation.rate.start<-0.000003
introduction.rate<-c(0.1,0.5)
results<-list()
nums<-0
seq.error<-0.01
time.lag<-NULL
min.perc<-20
init.infection<-20
mindiversity<-c(0,15)
maxdiversity<-c(0,25)

for (m in 1:2){
  for (p in 1:2){
    for (k in 1:2){
      
      if (p==1){
        init.sus.out<-550
      } else if (p==2){
        init.sus.out<-250
      } 
      if (m==1){
        init.sus.in<-550
      } else if (m==2){
        init.sus.in<-250
      } 
      
      repeat {
        W<-simulate_outbreak(init.sus = init.sus.out,init.sus.in=init.sus.in,
                             inf.rate=inf.rate.out[p], inf.rate.in = inf.rate.in[m],
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
                             min.cases = 10,
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
      
      W1<-W
      distmat_all<-gd(c(W$epidata[,5],W$epidata_in[,6]), W$libr, W$nuc, W$librstrains)
      colnames(distmat_all)<-c(W$epidata[,1],W$epidata_in[,1])
      row.names(distmat_all)<-c(W$epidata[,1],W$epidata_in[,1])
      
      ## Real outbreaks
      nums<-nums+1
      W$epidata_in[,1]<-W$epidata_in[,1]*10000
      W$epidata_in[,4]<-W$epidata_in[,4]*10000
      W$epidata_in[which(W$epidata_in[,5]==1),4]<-W$epidata_in[which(W$epidata_in[,5]==1),4]/10000
      distmat_in <- gd(W$epidata_in[,6], W$libr, W$nuc, W$librstrains)
      colnames(distmat_in)<-W$epidata_in[,1]
      row.names(distmat_in)<-W$epidata_in[,1]
      introduced_cases<-unique(W$epidata_in[which(W$epidata_in[,5]==1),4])
      no.introductions<-length(introduced_cases)
      dates.in<-data.frame(names=colnames(distmat_in),dates=0,intro=0,outbreak.in=0)
      for (i in 1:nrow(dates.in)){ # pick random sampling rate between infection and removal time
        dates.in$dates[i]<-sample(W$epidata_in[i,2]:W$epidata_in[i,3],1) 
      }
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
      ## setup igraph plots
      #change names
      allnames<-c(W$epidata[,1],W$epidata_in[,1])
      qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
      col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
      color.rand<- sample(col_vector,30)
      colors<-c(rep("grey",length(allnames)))
      outbreaks<-data.frame(table(dates.in$outbreak.in))
      outbreaks<-outbreaks[order(outbreaks$Freq,decreasing = T),]
      color.rand<-brewer.pal(nrow(outbreaks), "Set3")
      if (nrow(outbreaks)>12){
        color.rand<-rep(color.rand,10)
      }
      if (nrow(outbreaks)>0){
        if (nrow(outbreaks)==1){
          colors[which(allnames%in%dates.in[which(dates.in[,4]==
                                                    as.character(outbreaks$Var1[1])),1])]<-"red"
        } else if (nrow(outbreaks)>1){
          colors[which(allnames%in%dates.in[which(dates.in[,4]==
                                                    as.character(outbreaks$Var1[1])),1])]<-"red"
            #color.rand<-brewer.pal(nrow(outbreaks), "Set3")
          for (col in 2:nrow(outbreaks)){
            colors[which(allnames%in%dates.in[which(dates.in[,4]==
                                                      as.character(outbreaks$Var1[col])),1])]<-color.rand[col]
          }
        }
      }
      
      sizes<-c(rep(1,length(allnames)))
      shapes<-c(rep("circle",length(allnames)))
      sizes[which(allnames %in% W$epidata_in[,1])]<-10
      introductions<-unique(W$epidata_in[which(W$epidata_in[,5]==1),4])
      sizes[which(allnames %in% W$epidata[which(W$epidata[,1] %in% introductions),1])]<-5
      colors[which(allnames %in% W$epidata[which(W$epidata[,1] %in% introductions),1])]<-"yellow"
        shapes[which(allnames %in% W$epidata[which(W$epidata[,1] %in% introductions),1])]<-"square"
        results[[nums]]<- graph_from_data_frame(
          rbind(W$epidata[-c(1:init.infection), c(4,1)],W$epidata_in[,c(4,1)]),
          directed = T,
          vertices = data.frame(
            name = c(W$epidata[,1],W$epidata_in[,1]),
            inf.time =c(W$epidata[,2], W$epidata_in[,2]),
            g.id = c(W$epidata[,5],W$epidata_in[,6]), # genome IDs
            dist = distmat_all[c(W$epidata[,5],W$epidata_in[,6]), 1], # distances from genome 1
            color = colors, size=sizes, alpha=0.5,shape=shapes,
            label.color = "black")) %>%
          set_graph_attr("layout", layout.reingold.tilford) %>%
          set_graph_attr("par", W$parameters)
        
        
        
        ## Clustering at 100%
        colors1<-c(rep("grey",length(allnames)))
        results_sim.in<-cov2clusters_sim(distmat = distmat_in,dates = dates.in,probThreshold = 0.8) # run logit model
        clusters<-data.frame(table(results_sim.in[,2]))
        if (any(clusters$Var1=="-1")){
          clusters<-clusters[-which(clusters$Var1=="-1"),]
        }
        clusters<-clusters[order(clusters$Freq,decreasing = T),]
        if (nrow(clusters)==0){
          colors1[which(allnames %in% results_sim.in[,1])]<-"lightblue"
        } else {
          qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
          col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
          color.rand<- sample(col_vector,30)
          
          for (col in 1:nrow(clusters)){
            clustered_cases<-which(allnames%in%results_sim.in[which(results_sim.in[,2]==
                                                                      as.character(clusters$Var1[col])),1])
            
            colors_cluster<-data.frame(table(colors[clustered_cases]))
            colors_cluster<-colors_cluster[order(colors_cluster$Freq,decreasing = T),]
            used_cols<-which(colors_cluster$Var1 %in% colors1)
            if (length(used_cols)==0){
              colors1[clustered_cases]<-as.character(colors_cluster$Var1[1])
            } else if (length(used_cols)==nrow(colors_cluster)) {
              colors1[clustered_cases]<-color.rand[1]
              color.rand<-color.rand[-1]
            } else {
              colors_cluster<-colors_cluster[-used_cols,]
              colors1[clustered_cases]<-as.character(colors_cluster$Var1[1])
            }
          }
        }
        
        non_clustered<-which(results_sim.in[,2]=="-1")
        if (length(non_clustered)>0){
          for (i in 1:length(non_clustered)){
            outbreak.in.unclust<-dates.in$outbreak.in[which(dates.in$names %in% results_sim.in[non_clustered[i],1])]
            outbreak.length<-length(which(dates.in$outbreak.in==outbreak.in.unclust))
            if (outbreak.length>1){
              colors1[which(allnames %in% results_sim.in[non_clustered[i],1])]<-"lightblue"
            }
          } 
        }
        colors1[which(allnames %in% W$epidata[which(W$epidata[,1] %in% introductions),1])]<-"yellow"
          ##plot
        nums<-nums+1
        results[[nums]]<- graph_from_data_frame(
          #plot_test<-graph_from_data_frame(
          rbind(W$epidata[-c(1:init.infection), c(4,1)],W$epidata_in[,c(4,1)]),
          directed = T,
          vertices = data.frame(
            name = c(W$epidata[,1],W$epidata_in[,1]),
            inf.time =c(W$epidata[,2], W$epidata_in[,2]),
            g.id = c(W$epidata[,5],W$epidata_in[,6]), # genome IDs
            dist = distmat_all[c(W$epidata[,5],W$epidata_in[,6]), 1], # distances from genome 1
            color = colors1, size=sizes, alpha=0.5,shape=shapes,
            label.color = "black")) %>%
          set_graph_attr("layout", layout.reingold.tilford) %>%
          set_graph_attr("par", W$parameters)
        
        
        ## Clustering at 25%
        nums<-nums+1
        sampled.in<-which(rbinom(size_in, 1, 0.25)==1)
        distmat_in <- gd(W$epidata_in[sampled.in,6], W$libr, W$nuc, W$librstrains)
        colnames(distmat_in)<-W$epidata_in[sampled.in,1]
        row.names(distmat_in)<-W$epidata_in[sampled.in,1]
        
        colors1<-c(rep("gray90",length(allnames)))
        dates.in<-dates.in[sampled.in,]
        results_sim.in<-cov2clusters_sim(distmat = distmat_in,dates = dates.in,probThreshold = 0.8) # run logit model
        clusters<-data.frame(table(results_sim.in[,2]))
        if (any(clusters$Var1=="-1")){
          clusters<-clusters[-which(clusters$Var1=="-1"),]
        }
        clusters<-clusters[order(clusters$Freq,decreasing = T),]
        if (nrow(clusters)==0){
          colors1[which(allnames %in% results_sim.in[,1])]<-"lightblue"
        } else {
          qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
          col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
          color.rand<- sample(col_vector,30)
          
          for (col in 1:nrow(clusters)){
            clustered_cases<-which(allnames%in%results_sim.in[which(results_sim.in[,2]==
                                                                      as.character(clusters$Var1[col])),1])
            
            colors_cluster<-data.frame(table(colors[clustered_cases]))
            colors_cluster<-colors_cluster[order(colors_cluster$Freq,decreasing = T),]
            used_cols<-which(colors_cluster$Var1 %in% colors1)
            if (length(used_cols)==0){
              colors1[clustered_cases]<-as.character(colors_cluster$Var1[1])
            } else if (length(used_cols)==nrow(colors_cluster)) {
              colors1[clustered_cases]<-color.rand[1]
              color.rand<-color.rand[-1]
            } else {
              colors_cluster<-colors_cluster[-used_cols,]
              colors1[clustered_cases]<-as.character(colors_cluster$Var1[1])
            }
          }
        }
        
        non_clustered<-which(results_sim.in[,2]=="-1")
        if (length(non_clustered)>0){
          for (i in 1:length(non_clustered)){
            outbreak.in.unclust<-dates.in$outbreak.in[which(dates.in$names %in% results_sim.in[non_clustered[i],1])]
            outbreak.length<-length(which(dates.in$outbreak.in==outbreak.in.unclust))
            if (outbreak.length>1){
              colors1[which(allnames %in% results_sim.in[non_clustered[i],1])]<-"lightblue"
            } else {
              colors1[which(allnames %in% results_sim.in[non_clustered[i],1])]<-colors[which(allnames %in% results_sim.in[non_clustered[i],1])]
            }
          } 
        }
        colors1[which(allnames %in% W$epidata[which(W$epidata[,1] %in% introductions),1])]<-"yellow"
          sizes[which(sizes==10)]<-4
          sizes[which(allnames %in% dates.in$names)]<-10
          ##plot
          results[[nums]]<- graph_from_data_frame(
            #plot_test<-graph_from_data_frame(
            rbind(W$epidata[-c(1:init.infection), c(4,1)],W$epidata_in[,c(4,1)]),
            directed = T,
            vertices = data.frame(
              name = c(W$epidata[,1],W$epidata_in[,1]),
              inf.time =c(W$epidata[,2], W$epidata_in[,2]),
              g.id = c(W$epidata[,5],W$epidata_in[,6]), # genome IDs
              dist = distmat_all[c(W$epidata[,5],W$epidata_in[,6]), 1], # distances from genome 1
              color = colors1, size=sizes, alpha=0.5,shape=shapes,
              label.color = "black")) %>%
            set_graph_attr("layout", layout.reingold.tilford) %>%
            set_graph_attr("par", W$parameters)
          
          
    }
  }
}

setwd("./Test/")

saveRDS(results, file="Figure4Plots.RData")

plots<-readRDS("Figure4Plots.RData")

sampling<-c("Truth","100","25")
nums<-0
for (m in 1:2){
  for (p in 1:2){
    for (k in 1:2){
      for (n in 1:3){
        nums<-nums+1
        tiff(filename=paste0("R=",inf.rate.in[m]*1/removal.rate.all,"Samp=",
                             sampling[n],"max.div=",maxdiversity[p],"intro.rate",introduction.rate[k],
                             ".tiff"), res = 300, units = "cm",width = 20,height = 15)
        plot.igraph(plots[[nums]], 
                    vertex.label = NA, 
                    vertex.label.cex=1,
                    edge.arrow.size=0.1, edge.label.cex=2)
        
        title(paste0(
          "R  = ", inf.rate.in[m]*1/removal.rate.all ,
          "\nSampling = ", sampling[n], "\nIntro.rate = ", introduction.rate[k],
          "max diversity = ",maxdiversity[p] ,"SNP"),cex.main=1.5)
        dev.off()
      }
    }
  }
}
