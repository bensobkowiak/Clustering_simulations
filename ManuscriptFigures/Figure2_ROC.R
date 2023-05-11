require(seedy)
require(seqinr)
library(pROC)
require(ggplot2)
require(ggpubr)
options(stringsAsFactors = F)

setwd("~/R/Figure2")
source("~/R/simulate_outbreak.R")

# Run simulations
ref<-read.fasta("../MN908947.3.fasta",forceDNAtolower = F) # download Wuhan1 reference
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
mutation.rate<-c(0.000001,0.000005)
introduction.rate<-0
nums<-0
seq.error<-0.05
seq.error.rem<-0.05
seq.error.add<-0.01
plots<-list()

for (k in 1:4){
  for (m in 1:2){
    for (n in 1:2){
      nums<-nums+1
      init.sus = 200
      df<- data.frame(from=0, to=0, gen.dist=0, time.dist=0, d1=0)
      repeat {
        W<-simulate_outbreak(init.sus = init.sus,init.sus.in=0,
                             inf.rate=inf.rate.out[m], inf.rate.in = inf.rate.in, shape = gamma,
                             samp.freq = 10, mut.rate.out = mutation.rate[n], 
                             intr.rate = introduction.rate,
                             rem.rate=removal.rate.all,min.cases.in=0,min.cases = 2,
                             ref.strain = ref,g.len = 29903)
        
        size_out<-length(W$epidata[,5])
        
        ## sampling by sampling proportion and set up snp distance
        sampled.out<-which(rbinom(size_out, 1, sampling.rate.out[k])==1)
        if (length(sampled.out)>1){
          
          for (i in 1:length(sampled.out)){ 
            if (sampled.out[i]!=1){  # add random noise
              genome.num<-W$epidata[sampled.out[i],5]
              genome.nuc<-W$nuc[genome.num][[1]]
              change<-which(rbinom(length(genome.nuc),1,seq.error.rem)==1) ## remove SNPs based on seq.error.rem rate
              if (length(change)>0){
                if (length(which(genome.num==W$epidata[sampled.out,5]))==1){
                  W$nuc[genome.num][[1]]<-genome.nuc[-c(change)]
                  W$libr[genome.num][[1]]<-W$libr[genome.num][[1]][-c(change)]
                } else { # make new genome if shared by another individual
                  newnum<-max(W$epidata[,5])+1
                  W$nuc[newnum][[1]]<-genome.nuc[-c(change)]
                  W$libr[newnum][[1]]<-W$libr[genome.num][[1]][-c(change)]
                  W$epidata[sampled.out[i],5]<-newnum
                  W$librstrains<-c(W$librstrains,newnum)
                }
              }
            }
          }
          
          distmat_out <- gd(W$epidata[sampled.out,5], W$libr, W$nuc, W$librstrains)
          colnames(distmat_out)<-W$epidata[sampled.out,1]
          row.names(distmat_out)<-W$epidata[sampled.out,1]
          
          # Make sampling dates 
          dates.out<-data.frame(names=W$epidata[sampled.out,1],dates=W$sampledata[sampled.out,3])
          
          ## ROC curve setup
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
          df_res<- data.frame(from, to, gen.dist, time.dist, d1)
          df<-rbind(df,df_res)
        }
        if (length(unique(c(df$from,df$to)))>200) break
      }
      df<-df[-1,]
      
      # logit model with date
      coef_vec_d1 <- c(3,-0.66 ,-0.075)
      z_d1 <- coef_vec_d1[1] + coef_vec_d1[2]*df$gen.dist + coef_vec_d1[3]*df$time.dist
      df$pred.of.d1 <- 1/(1 + exp(-z_d1))
      
      # logit model without date
      coef_vec_d3 <- c(3,-0.66)
      z_d3 <- coef_vec_d3[1] + coef_vec_d3[2]*df$gen.dist 
      df$pred.of.d3 <- 1/(1 + exp(-z_d3))
      
      # SNP model 2 SNP threshold
      df$pred.of.d4<-0
      df$pred.of.d4[which(df$gen.dist<=2)]<-1
      
      # SNP model 1 SNP threshold
      df$pred.of.d5<-0
      df$pred.of.d5[which(df$gen.dist<=1)]<-1
      
      # SNP model 3 SNP threshold
      df$pred.of.d6<-0
      df$pred.of.d6[which(df$gen.dist<=3)]<-1
      
      # SNP model 4 SNP threshold
      df$pred.of.d7<-0
      df$pred.of.d7[which(df$gen.dist<=4)]<-1
      
      # SNP model 5 SNP threshold
      df$pred.of.d8<-0
      df$pred.of.d8[which(df$gen.dist<=5)]<-1
      
      roc_info_logit_d1 <- roc(df$d1, df$pred.of.d1, legacy.axes=TRUE)
      roc_df_logit_d1 <- data.frame( 
        TPR_d1 = roc_info_logit_d1$sensitivities,
        FPR_d1 = 1- roc_info_logit_d1$specificities,
        thr_d1 = roc_info_logit_d1$thresholds)
      
      roc_info_logit_d3 <- roc(df$d1, df$pred.of.d3, legacy.axes=TRUE)
      roc_df_logit_d3 <- data.frame( 
        TPR_d3 = roc_info_logit_d3$sensitivities, 
        FPR_d3 = 1 - roc_info_logit_d3$specificities,
        thr_d3 = roc_info_logit_d3$thresholds)
      
      roc_info_SNP_d4 <- roc(df$d1, df$pred.of.d4, legacy.axes=TRUE)
      roc_df_SNP_d4 <- data.frame( 
        TPR_d4 = roc_info_SNP_d4$sensitivities, 
        FPR_d4 = 1 - roc_info_SNP_d4$specificities,
        thr_d4 = roc_info_SNP_d4$thresholds)
      
      roc_info_SNP_d5 <- roc(df$d1, df$pred.of.d5, legacy.axes=TRUE)
      roc_df_SNP_d5 <- data.frame( 
        TPR_d5 = roc_info_SNP_d5$sensitivities, 
        FPR_d5 = 1 - roc_info_SNP_d5$specificities,
        thr_d5 = roc_info_SNP_d5$thresholds)
      
      roc_info_SNP_d6 <- roc(df$d1, df$pred.of.d6, legacy.axes=TRUE)
      roc_df_SNP_d6 <- data.frame( 
        TPR_d6 = roc_info_SNP_d6$sensitivities, 
        FPR_d6 = 1 - roc_info_SNP_d6$specificities,
        thr_d6 = roc_info_SNP_d6$thresholds)
      
      roc_info_SNP_d7 <- roc(df$d1, df$pred.of.d7, legacy.axes=TRUE)
      roc_df_SNP_d7 <- data.frame( 
        TPR_d7 = roc_info_SNP_d7$sensitivities, 
        FPR_d7 = 1 - roc_info_SNP_d7$specificities,
        thr_d7 = roc_info_SNP_d7$thresholds)
      
      roc_info_SNP_d8 <- roc(df$d1, df$pred.of.d8, legacy.axes=TRUE)
      roc_df_SNP_d8 <- data.frame( 
        TPR_d8 = roc_info_SNP_d8$sensitivities, 
        FPR_d8 = 1 - roc_info_SNP_d8$specificities,
        thr_d8 = roc_info_SNP_d8$thresholds)
      
      
      # Plot ROC
      plots[[nums]]<-ggplot() +
        geom_line(data = roc_df_logit_d1,
                  mapping = aes(x = FPR_d1, y = TPR_d1, colour = "Logit model with dates",linetype="Logit model with dates")) +
        geom_line(data = roc_df_logit_d3,
                  mapping = aes(x = FPR_d3, y = TPR_d3, colour = "Logit model without dates",linetype="Logit model without dates")) +
        geom_line(data = roc_df_SNP_d6,
                  mapping = aes(x = FPR_d6, y = TPR_d6, colour = "3 SNP threshold",linetype="3 SNP threshold")) +
        geom_line(data = roc_df_SNP_d4,
                  mapping = aes(x = FPR_d4, y = TPR_d4, colour = "2 SNP threshold",linetype="2 SNP threshold")) +
        geom_line(data = roc_df_SNP_d5,
                  mapping = aes(x = FPR_d5, y = TPR_d5, colour = "1 SNP threshold",linetype="1 SNP threshold")) +
        
        labs(title = paste0("R0=",inf.rate.out[m]*1/removal.rate.all,
                            ", Mut.rate=", mutation.rate[n],
                            ", sampling.prop=",sampling.rate.out[k]),
             x = "False positive rate", y = "True positive rate",
             colour = "Clustering method",linetype = "Clustering method") +
        theme(plot.title = element_text(hjust = 0.5,size = 8)) + 
        scale_color_manual(values=c("blue","blue","blue","red","red")) +
        scale_linetype_manual(values=c("solid", "dashed","dotted","solid", "dashed"))
      
    }
  }
}

saveRDS(plots, file="ROC_curves_test.RData")
print("finished")

plots<-readRDS("ROC_curves_test.RData")
tiff(filename="ROC_curves_final_figure2_mutrate.tiff", res = 300, units = "cm",width = 30,height = 20)
p5<-ggarrange(plots[[1]]+ xlab(NULL) +labs(title=NULL) + annotate("text",x = 0.7, y = 0.25, label ="R = 1.25 \nMutation rate low",size = 3) +theme_bw(),
              plots[[2]]+ xlab(NULL) + ylab(NULL)+labs(title=NULL) + annotate("text",x = 0.7, y = 0.25, label ="R = 1.25 \nMutation rate high",size = 3) +theme_bw(),
              plots[[3]]+ xlab(NULL) + ylab(NULL)+labs(title=NULL) + annotate("text",x = 0.7, y = 0.25, label ="R = 3 \nMutation rate low",size = 3) +theme_bw(),
              plots[[4]]+ xlab(NULL) + ylab(NULL)+labs(title=NULL) + annotate("text",x = 0.7, y = 0.25, label ="R = 3 \nMutation rate high",size = 3) +theme_bw(),
              plots[[5]]+ xlab(NULL)+labs(title=NULL) + annotate("text",x = 0.7, y = 0.25, label ="R = 1.25 \nMutation rate low",size = 3) +theme_bw(),
              plots[[6]]+ xlab(NULL) + ylab(NULL)+labs(title=NULL) + annotate("text",x = 0.7, y = 0.25, label ="R = 1.25 \nMutation rate high",size = 3) +theme_bw(),
              plots[[7]]+ xlab(NULL) + ylab(NULL)+labs(title=NULL) + annotate("text",x = 0.7, y = 0.25, label ="R = 3 \nMutation rate low",size = 3) +theme_bw(),
              plots[[8]]+ xlab(NULL) + ylab(NULL)+labs(title=NULL) + annotate("text",x = 0.7, y = 0.25, label ="R = 3 \nMutation rate high",size = 3) +theme_bw(),
              plots[[9]]+ xlab(NULL)+labs(title=NULL) + annotate("text",x = 0.7, y = 0.25, label ="R = 1.25 \nMutation rate low",size = 3) +theme_bw(),
              plots[[10]]+ xlab(NULL) + ylab(NULL)+labs(title=NULL) + annotate("text",x = 0.7, y = 0.25, label ="R = 1.25 \nMutation rate high",size = 3) +theme_bw(),
              plots[[11]]+ xlab(NULL) + ylab(NULL)+labs(title=NULL) + annotate("text",x = 0.7, y = 0.25, label ="R = 3 \nMutation rate low",size = 3) +theme_bw(),
              plots[[12]]+ xlab(NULL) + ylab(NULL)+labs(title=NULL) + annotate("text",x = 0.7, y = 0.25, label ="R = 3 \nMutation rate high",size = 3) +theme_bw(),
              plots[[13]]+labs(title=NULL) + annotate("text",x = 0.7, y = 0.25, label ="R = 1.25 \nMutation rate low",size = 3) +theme_bw(),
              plots[[14]] + ylab(NULL)+labs(title=NULL) + annotate("text",x = 0.7, y = 0.25, label ="R = 1.25 \nMutation rate high",size = 3) +theme_bw(),
              plots[[15]]+ ylab(NULL)+labs(title=NULL) + annotate("text",x = 0.7, y = 0.25, label ="R = 3 \nMutation rate low",size = 3) +theme_bw(),
              plots[[16]]+ ylab(NULL)+labs(title=NULL) + annotate("text",x = 0.7, y = 0.25, label ="R = 3 \nMutation rate high",size = 3) +theme_bw(),
              ncol=4, nrow=4, common.legend = TRUE, legend="bottom",label.x = "test")
plot(p5)
dev.off()


p1<-ggarrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]],plots[[6]],
              plots[[7]],plots[[8]],plots[[9]],plots[[10]],plots[[11]],plots[[12]],
              ncol=3, nrow=4, common.legend = TRUE, legend="bottom")
annotate_figure(p1, top = text_grob("Sampling proportion = 10%", 
                                    color = "black", face = "bold", size = 14))

p2<-ggarrange(plots[[13]],plots[[14]],plots[[15]],plots[[16]],plots[[17]],plots[[18]],
              plots[[19]],plots[[20]],plots[[21]],plots[[22]],plots[[23]],plots[[24]],
              ncol=3, nrow=4, common.legend = TRUE, legend="bottom")
annotate_figure(p2, top = text_grob("Sampling proportion = 25%", 
                                    color = "black", face = "bold", size = 14))

p3<-ggarrange(plots[[25]],plots[[26]],plots[[27]],plots[[28]],plots[[29]],plots[[30]],
              plots[[31]],plots[[32]],plots[[33]],plots[[34]],plots[[35]],plots[[36]],
              ncol=3, nrow=4, common.legend = TRUE, legend="bottom")
annotate_figure(p3, top = text_grob("Sampling proportion = 50%", 
                                    color = "black", face = "bold", size = 14))

p4<-ggarrange(plots[[37]],plots[[38]],plots[[39]],plots[[40]],plots[[41]],plots[[42]],
              plots[[43]],plots[[44]],plots[[45]],plots[[46]],plots[[47]],plots[[48]],
              ncol=3, nrow=4, common.legend = TRUE, legend="bottom")
annotate_figure(p4, top = text_grob("Sampling proportion = 100%", 
                                    color = "black", face = "bold", size = 14))

p<-plots[[2]]

p + 
  scale_linetype_manual(values = c("1 SNP threshold" = "solid", 
                                   "2 SNP threshold" = "solid",
                                   "Logit model weak date" = "dotted",
                                   "Logit model no date" = "dotted"))

for (i in 1:16){
  p<-plots[[i]]
  p$layers<-p$layers[-2]
  plots[[i]]<-p
}
