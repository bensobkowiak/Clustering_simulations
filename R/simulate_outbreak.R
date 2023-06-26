# This file is based on the simulate_outbreak function by here:
# https://github.com/Pouya-Haghmaram/Clustering-outbreaks, which is an altered version 
# of simulateoutbreak.R in seedy https://github.com/cran/seedy/blob/master/R/simulateoutbreak.R

simulate_outbreak <- # Outbreak simulation function 
  function(init.sus=10, init.sus.in=10, inf.rate=1, inf.rate.in=1, intr.rate=0,rem.rate=0.5, 
           mut.rate.out=0.0008,mut.rate.in=0.0008, 
           nmat=NULL, equi.pop=10000, 
           init.inf=1, inoc.size=1, time.lag.outside=NULL, 
           min.init.dist=0,max.init.dist=0,
           min.perc.outside.inf = NULL,shape.infect=6,rate.infect=3,
           samples.per.time=1, samp.schedule="random", 
           samp.freq=500, full=TRUE, min.cases=1, min.cases.in=1, 
           feedback=500, g.len=10000, 
           ref.strain=NULL, ...) {
    
    # Warning messages
    
    if (!full) {
      stop("full must be TRUE.")
    }
    
    if (init.sus%%1!=0 || init.sus<1) {
      stop("init.sus must be a postive integer.")
    }
    
    if (init.sus.in%%1!=0 || init.sus.in<0) {
      stop("init.sus.in must be a non-negative integer.")
    }
    
    if (inoc.size%%1!=0 || inoc.size<1) {
      stop("inoc.size must be a postive integer.")
    }
    
    if (samp.freq%%1!=0 || samp.freq<1) {
      stop("samp.freq must be a postive integer.")
    }
    if (feedback%%1!=0 || feedback<1) {
      stop("feedback must be a postive integer.")
    }
    if (!is.null(time.lag.outside) & !is.null(min.perc.outside.inf)) {
      stop("Can only specify time lag or minimum percent infected outside.")
    }
    if (g.len%%1!=0 || g.len<1) {
      stop("g.len (genome length) must be a postive integer.")
    }
    
    if (min.cases%%1!=0 || min.cases<1 || min.cases>init.sus+init.inf) {
      stop("min.case must be a postive integer not greater than init.sus+init.inf.")
    }
    
    if (min.cases.in%%1!=0 || min.cases.in<0 || min.cases.in>init.sus.in) {
      stop("min.cases.in must be a non-negative integer not greater than init.sus.in.")
    }
    
    if (samples.per.time%%1!=0 || samples.per.time<1) {
      stop("samples.per.time must be a postive integer.")
    }
    
    if (inf.rate<=0) {
      stop("inf.rate (infection rate) must be greater than zero.")
    }
    if (inf.rate.in<=0) {
      stop("inf.rate.in (infection rate inside) must be greater than zero.")
    }
    
    if (rem.rate<=0) {
      stop("rem.rate (removal rate) must be greater than zero.")
    }
    
    if (intr.rate<0 | intr.rate>1) {
      stop("intr.rate (introduction rate) must be between 0 and 1.")
    }
    
    if (!is.null(nmat)) {
      if (!is.matrix(nmat)) {
        stop("nmat must be a matrix.")
      } else {
        if (nrow(nmat)!=init.sus+init.inf || ncol(nmat)!=init.sus+init.inf) {
          stop("nmat must have init.sus+init.inf rows and columns.")
        } else if (sum(nmat<0)>0) {
          stop("All entries of nmat must be >= 0.")
        }
      }
    }
    
    if (equi.pop%%1!=0 || equi.pop<=0) {
      stop("equi.pop (equilibrium population size) must be a postive integer.")
    }
    
    if (!samp.schedule=="random") {
      stop("samp.schedule must be 'random'.")
    }
    
    #########################
    
    trigger <- FALSE
    trigger1 <- FALSE
    
    at <- 0 #attempt
    
    while (!trigger | !trigger1) { # Repeat if < min.cases are infected or 
      # < min.cases.in are infected inside
      newinfect <- 0
      
      at <- at+1
      #cat("Attempt ", at, "\n", sep="")
      
      eff.cur.inf <- NULL
      
      cur.inf <- 1:init.inf 
      cur.inf.in <- NULL 
      
      cur.sus <- init.inf+(1:init.sus) # vector of susceptible person IDs (outside)
      cur.sus.in <- 1:init.sus.in # vector of susceptible person IDs (inside)
      
      time <- 0 # in days
      
      inf.times <- rep(0,init.inf) # vector of infection times (outside)
      inf.times.in <- NULL # vector of infection times (inside)
      
      rec.times <- rgeom(init.inf,rem.rate)+2 # vector of removal times (outside)
      rec.times.in <- NULL # vector of removal times (inside)
      
      tot.inf <- init.inf # number of all infected people (outside)
      tot.inf.in <- 0 # number of all infected people (inside)
      
      inf.ID <- 1:init.inf 
      inf.ID.in <- NULL
      
      if (samp.schedule=="random") {
        sample.times <- NULL
        for (i in 1:init.inf) {
          sample.times <- c(sample.times, sample(1:(rec.times[i]-1),1)) 
        }
      } else {
        sample.times <- rep(samp.freq, init.inf)
      }
      
      inf.source <- rep(0,init.inf) # source of infection for each individual
      inf.source.in <- NULL
      introduction <- NULL # yes=1, no=0
      
      if (is.null(ref.strain)) {
        ref.strain <- sample(1:4, g.len, replace=T) # reference strain
      } else {
        g.len <- length(ref.strain)
      }
      
      totcurstrains <- 1:init.inf # current list of strains
      uniquestrains <- 1:init.inf # Number of unique strain types
      libr <- list() # list of mutation locations for each genotype
      mut.nuc <- list() # nucleotides at mutation locations
      freq.log <- list() # List of strain frequencies for each infective
      strain.log <- list() # Strain IDs for each within host population
      genomeID.in <- NULL # List of genome IDs inside
      
      # Initialize logs
      for (i in 1:init.inf) {
        if (sample.times[i]>rec.times[i]) {
          sample.times[i] <- Inf 
        }
        if (init.inf == 1) {
          libr[[i]] <- NA 
          mut.nuc[[i]] <- NA      
        } else {
          start.dist<-sample(min.init.dist:max.init.dist,1)# random distance from reference up to max distance
          if (start.dist==0){
            libr[[i]] <- NA 
            mut.nuc[[i]] <- NA  
          } else {
          libr[[i]] <- sample(g.len,start.dist) # take start.dist number of loci
          startsubs<-numeric()
          ref_nucs<-ref.strain[libr[[i]]]
          for (subs in 1:start.dist){
            newnuc<-sample((1:4)[-ref_nucs[subs]],1)
            startsubs<-c(startsubs,newnuc)
          }
          mut.nuc[[i]] <- startsubs
          }
        }
        freq.log[[i]] <- 1
        strain.log[[i]] <- i
      }
      
      for (i in (init.inf+1):(init.sus+init.inf)) {
        freq.log[[i]] <- 0
        strain.log[[i]] <- 0
      }
      types <- init.inf # Cumulative number of strain types
      
      #Sample logs
      if (full) {
        obs.freq <- list()
        obs.strain <- list()
        pID <- NULL
      } else {
        sampleWGS <- NULL
        samplepick <- NULL
      }
      sampletimes <- NULL 
      sampleID <- NULL 
      
      # Cycle through bacterial generations until epidemic ceases
      while (length(cur.inf) + length(cur.inf.in) > 0) { 
        time <- time+1
        
        if (time%in%rec.times) { # recovery? (outside)
          recover <- inf.ID[which(rec.times==time)] # who has recovered?
          cur.inf <- cur.inf[-which(cur.inf%in%recover)] # remove infective(s)
          
          for (r in 1:length(recover)) {
            strain.log[[recover[r]]] <- 0
            freq.log[[recover[r]]] <- 0
          }
          
          if (length(cur.inf)==0) { # If no more infectives
            #cat("t=", time, ", S=", length(cur.sus), ", I=", length(cur.inf), 
            #    ", total genotypes=0\n", sep="")
          }
        }
        
        if (time%in%rec.times.in) { # recovery? (inside)
          recover.in <- inf.ID.in[which(rec.times.in==time)] # who has recovered?
          cur.inf.in <- cur.inf.in[-which(cur.inf.in%in%recover.in)] # remove infective(s)
          
          if (length(cur.inf.in)==0) { # If no more infectives
            #cat("\nt=", time, ", I(inside)=", length(cur.inf.in),"\n\n")
          }
        }
        
        if (time%%feedback==0 && length(cur.inf)>0) { # output current status every x generations
          #cat("t=", time, ", S=", length(cur.sus), ", I=", length(cur.inf), 
          #    ", total genotypes=", length(unique(as.numeric(unlist(strain.log)))), ", next rec time=", 
          #    min(rec.times[which(rec.times>time)]), sep="")
          if (length(cur.sus)>0) {
            #  cat("\n")
          } else {
            #  cat(", final removal time=", max(rec.times), "\n", sep="")
          }
        }
        
        # calculate force of infection (outside)
        pinf <- rep(0, length(cur.sus)) 
        
        if(length(cur.sus)!=0){
          if (is.null(nmat)) {
            curinfrate <- inf.rate*length(cur.inf)/(init.inf+init.sus)
            pinf <- rbinom(length(cur.sus),1,curinfrate)
          } else {
            inf.mat <- nmat*inf.rate/(init.inf+init.sus)
            if (length(cur.inf)>1) {
              curinfrate <- apply(inf.mat[cur.sus,cur.inf],1,sum)
            } else {
              curinfrate <- inf.mat[cur.sus,]
            }
            pinf <- rbinom(length(cur.sus),1,curinfrate)
          }
        }
        
        # infection? (outside)
        if (sum(pinf)>0) {
          for (kp in 1:sum(pinf)) {
            tot.inf <- tot.inf+1
            newinfect <- cur.sus[(which(pinf==1)[kp])-(kp-1)] 
            
            if (length(cur.inf)==1) {
              inf.source <- c(inf.source, cur.inf)
            } else if (is.null(nmat)) {
              #inf.source <- c(inf.source, sample(cur.inf,1)) # sample source at random
              # new code - non-random instead based on gamma distribution
              put.inf.time<- round(rgamma(n=1, shape=shape.infect, rate=rate.infect),0) # putative infection time
              removed<-which(!inf.ID %in% cur.inf)
              if (length(removed)>0){
                curr.inf.time<-inf.times[-which(!inf.ID %in% cur.inf)]
                time.since.infect<- time-curr.inf.time
              } else {
                time.since.infect<- time-inf.times
              }
              ## Find infection source(s) closest to estimated infection time
              best.inf.time<-which(abs(time.since.infect[which(time.since.infect!=0)] - put.inf.time)==
                                     min(abs(time.since.infect[which(time.since.infect!=0)] - put.inf.time)))
              if (length(best.inf.time)==1){
                inf.source <- c(inf.source, cur.inf[best.inf.time])
              } else { ## pick random strain if multiple infected on same day
                inf.source <- c(inf.source, sample(cur.inf[best.inf.time],1))
              }
              
            } else {
              inf.source <- c(inf.source, sample(cur.inf, 1, prob=inf.mat[newinfect,cur.inf]))
            }
            
            inf.ID <- c(inf.ID, newinfect)
            cur.inf <- c(cur.inf, newinfect) # add to current infectives
            cur.sus <- cur.sus[-which(cur.sus==newinfect)] # Remove susceptible
            inf.times <- c(inf.times, time) # new infection time
            recover.time<-time + 1 + rgeom(1,rem.rate)
            rec.times <- c(rec.times, recover.time) # Don't recover today!
            
            if (samp.schedule == "individual") {
              sample.times <- c(sample.times, time+samp.freq)
            } else if (samp.schedule == "calendar") {
              sample.times <- c(sample.times, ceiling(time/samp.freq)*samp.freq)
            } else if (samp.schedule == "random") {
              if (rec.times[tot.inf]>time+1) {
                sample.times <- c(sample.times, sample(time:(rec.times[tot.inf]-1),1))
              } else {
                sample.times <- c(sample.times, time)
              }
            }
            if (sample.times[tot.inf]>=rec.times[tot.inf]) {
              sample.times[tot.inf] <- Inf
            }
            # pass on strain
            src <- inf.ID[which(inf.ID==inf.source[tot.inf])] # Source of infection
            if (length(strain.log[[src]])==1) { # if source has clonal infection
              inoc.samp <- rep(strain.log[[src]], inoc.size)
              if (0%in%inoc.samp) {
                stop("Zeroes in inoculum")
              }
            } else {
              inoc.samp <- sample(strain.log[[src]], inoc.size, 
                                  prob=freq.log[[src]], replace=T) # take random sample
              if (0%in%inoc.samp) {
                stop("Zeroes in inoculum")
              }
            }
            
            strain.log[[newinfect]] <- unique(inoc.samp) # distinct types in new infection
            f <- numeric(length(unique(inoc.samp)))
            k <- 1
            for (i in unique(inoc.samp)) {
              f[k] <- sum(inoc.samp==i)
              k <- k+1
            }
            freq.log[[newinfect]] <- f # frequency of types
            
            for (str in 1:length(strain.log[[newinfect]])) {
              for (fre in 1:length(freq.log[[newinfect]])) {
                inf.times[which(inf.ID==inf.source[tot.inf])]
                num_sub <- rbinom(1, g.len, (mut.rate.out)*
                                    (time-inf.times[which(inf.ID==inf.source[tot.inf])])) # scale mutation rate 
                                                                                          # by number of days infection persists
                
                if(num_sub>0){
                  types <- types+1
                  temp.libr <- libr[[which(totcurstrains==strain.log[[newinfect]][str])]]
                  temp.nuc <- mut.nuc[[which(totcurstrains==strain.log[[newinfect]][str])]]
                  for(count1 in 1:num_sub){
                    mut.loc <- sample(g.len, 1)
                    
                    if (mut.loc %in% temp.libr) { # if mutation at existing location
                      kn <- which(temp.libr==mut.loc)
                      temp.nuc[kn] <- sample((1:4)[-temp.nuc[kn]], 1)
                    } else { 
                      temp.nuc <- c(temp.nuc, 
                                    sample((1:4)[-ref.strain[mut.loc]], 1))
                      temp.libr <- c(temp.libr, mut.loc)
                    }
                    if (sum(is.na(temp.libr))>0) {
                      temp.nuc <- temp.nuc[!is.na(temp.nuc)]
                      temp.libr <- temp.libr[!is.na(temp.libr)]
                    }
                  }
                  
                  libr[[length(totcurstrains)+1]] <- temp.libr
                  mut.nuc[[length(totcurstrains)+1]] <- temp.nuc
                  
                  strain.log[[newinfect]] <- c(strain.log[[newinfect]], types)
                  freq.log[[newinfect]] <- c(freq.log[[newinfect]], 1)
                  totcurstrains <- c(totcurstrains, types)
                }
              }
            }
            
            if(length(strain.log[[newinfect]])>1){
              if(length(strain.log[[newinfect]])==2){
                ran_index <- 2
              } else{
                ran_index <- sample(2:length(strain.log[[newinfect]]), 1)
              }
              strain.log[[newinfect]] <- strain.log[[newinfect]][ran_index]
              freq.log[[newinfect]] <- freq.log[[newinfect]][ran_index]
            }
          }
        }
        
        ### Add lag time or minimum proportion of outside pop. infected before inside outbreaks begin
        if ((length(time.lag.outside)==0 || (length(time.lag.outside)!=0 & time>= 
                                             time.lag.outside)) & 
            (length(min.perc.outside.inf)==0 || 
             (length(min.perc.outside.inf)!=0 & length(inf.ID)>= 
              init.sus/100*min.perc.outside.inf))){
          
          pinf.in <- rep(0, length(cur.sus.in))
          
          # calculate force of infection (inside)
          if(length(cur.sus.in)!=0 & init.sus.in!=0){
            curinfrate1 <- inf.rate.in*length(cur.inf.in)/(init.sus.in)
            pinf.in <- rbinom(length(cur.sus.in),1,curinfrate1)
          }
          
          #infection inside?
          if (sum(pinf.in)>0) {
            for (kp in 1:sum(pinf.in)) {
              tot.inf.in <- tot.inf.in+1
              
              newinfect.in <- cur.sus.in[(which(pinf.in==1)[kp])-(kp-1)]
              
              if (length(cur.inf.in)==1) {
                inf.source.in <- c(inf.source.in, cur.inf.in)
                introduction <- c(introduction, 0)
              } 
              if (length(cur.inf.in)>1){
                #inf.source.in <- c(inf.source.in, sample(cur.inf.in,1)) # sample source at random
                # Same as outside infection - infection source based on gamma distribution of infection times
                put.inf.time<- round(rgamma(n=1, shape=shape.infect, rate=rate.infect),0) # putative infection time
                removed<-which(!inf.ID.in %in% cur.inf.in)
                if (length(removed)>0){
                  curr.inf.time<-inf.times.in[-which(!inf.ID.in %in% cur.inf.in)]
                  time.since.infect<- time-curr.inf.time
                } else {
                  time.since.infect<- time-inf.times.in
                }
                best.inf.time<-which(abs(time.since.infect[which(time.since.infect!=0)] - put.inf.time)==
                                       min(abs(time.since.infect[which(time.since.infect!=0)] - put.inf.time)))
                if (length(best.inf.time)==1){
                  inf.source.in <- c(inf.source.in, cur.inf.in[best.inf.time])
                } else {
                  inf.source.in <- c(inf.source.in, sample(cur.inf.in[best.inf.time],1))
                }
                introduction <- c(introduction, 0)
              }
              
              inf.ID.in <- c(inf.ID.in, newinfect.in)
              cur.inf.in <- c(cur.inf.in, newinfect.in) # add to current infectives
              cur.sus.in <- cur.sus.in[-which(cur.sus.in==newinfect.in)] # Remove susceptible
              inf.times.in <- c(inf.times.in, time) # new infection time
              recover.time <- time + 1 + rgeom(1,rem.rate)
              rec.times.in <- c(rec.times.in, recover.time) 
              
              # pass on strain
              src1 <- inf.ID.in[which(inf.ID.in==inf.source.in[tot.inf.in])] # Source of infection
              genomeID.in[tot.inf.in] <- genomeID.in[which(inf.ID.in==inf.source.in[tot.inf.in])]
              # mutate existing strains for each individual
              num_sub <- rbinom(1, g.len, (mut.rate.in)*
                                  (time-inf.times.in[which(inf.ID.in==inf.source.in[tot.inf.in])]))
              
              
              if(num_sub>0){
                types <- types+1
                temp.libr <- libr[[which(totcurstrains==genomeID.in[tot.inf.in])]]
                temp.nuc <- mut.nuc[[which(totcurstrains==genomeID.in[tot.inf.in])]]
                for(count1 in 1:num_sub){
                  mut.loc <- sample(g.len, 1)
                  
                  if (mut.loc %in% temp.libr) { # if mutation at existing location
                    kn <- which(temp.libr==mut.loc)
                    temp.nuc[kn] <- sample((1:4)[-temp.nuc[kn]], 1)
                  } else { 
                    temp.nuc <- c(temp.nuc, 
                                  sample((1:4)[-ref.strain[mut.loc]], 1))
                    temp.libr <- c(temp.libr, mut.loc)
                  }
                  if (sum(is.na(temp.libr))>0) {
                    temp.nuc <- temp.nuc[!is.na(temp.nuc)]
                    temp.libr <- temp.libr[!is.na(temp.libr)]
                  }
                }
                
                libr[[length(totcurstrains)+1]] <- temp.libr
                mut.nuc[[length(totcurstrains)+1]] <- temp.nuc
                
                genomeID.in[tot.inf.in] <- types
                totcurstrains <- c(totcurstrains, types)
              }
            }
          }
          
          #introduction?
          introduction_time <- rbinom(1,1,intr.rate)*time
          
          if(time==introduction_time & length(cur.sus.in)!=0 & length(cur.inf) !=0 & init.sus.in !=0){# time of a new introduction
            
            if(length(cur.inf)==1){
              infector_introduction <- cur.inf
            } else{
              infector_introduction <- sample(cur.inf, 1)
            }
            
            if(length(cur.sus.in)==1){
              infectee_introduction <- cur.sus.in
            } else{
              infectee_introduction <- sample(cur.sus.in, 1)
            }
            
            tot.inf.in <- tot.inf.in+1
            inf.ID.in <- c(inf.ID.in, infectee_introduction)
            cur.inf.in <- c(cur.inf.in, infectee_introduction) # add to current infectives
            cur.sus.in <- cur.sus.in[-which(cur.sus.in==infectee_introduction)] # Remove susceptible
            inf.times.in <- c(inf.times.in, time) # new infection time
            inf.source.in <- c(inf.source.in, infector_introduction)
            introduction <- c(introduction, 1)
            recover.time <- time + 1 + rgeom(1,rem.rate)
            rec.times.in <- c(rec.times.in, recover.time) 
            genomeID.in[tot.inf.in] <- strain.log[[infector_introduction]]
            
            # mutate existing strains for each individual
            num_sub <- rbinom(1, g.len, (mut.rate.in)*
                                (time-inf.times[which(inf.ID==inf.source.in[tot.inf.in])]))
            
            if(num_sub>0){
              types <- types+1
              temp.libr <- libr[[which(totcurstrains==genomeID.in[tot.inf.in])]]
              temp.nuc <- mut.nuc[[which(totcurstrains==genomeID.in[tot.inf.in])]]
              
              for(count1 in 1:num_sub){
                mut.loc <- sample(g.len, 1)
                
                if (mut.loc %in% temp.libr) { # if mutation at existing location
                  kn <- which(temp.libr==mut.loc)
                  temp.nuc[kn] <- sample((1:4)[-temp.nuc[kn]], 1)
                } else { 
                  temp.nuc <- c(temp.nuc, 
                                sample((1:4)[-ref.strain[mut.loc]], 1)
                  )
                  temp.libr <- c(temp.libr, mut.loc)
                }
                if (sum(is.na(temp.libr))>0) {
                  temp.nuc <- temp.nuc[!is.na(temp.nuc)]
                  temp.libr <- temp.libr[!is.na(temp.libr)]
                }
              }
              
              libr[[length(totcurstrains)+1]] <- temp.libr
              mut.nuc[[length(totcurstrains)+1]] <- temp.nuc
              
              genomeID.in[tot.inf.in] <- types
              totcurstrains <- c(totcurstrains, types)
            }
          }
        }
        ### end outbreak
        
        # take samples, make observations (outside)
        if (time%in%sample.times) {
          smpat <- inf.ID[which(sample.times==time)]
          for (i in smpat) {
            if (full) {
              n <- length(obs.freq)+1
              obs.freq[[n]] <- freq.log[[i]]
              obs.strain[[n]] <- strain.log[[i]]
              sampleID <- c(sampleID, n)
              pID <- c(pID, i)
              sampletimes <- c(sampletimes, time)
            } else {
              for (j in 1:samples.per.time) {
                if (length(strain.log[[i]])==1) {
                  pickgrp <- strain.log[[i]]
                } else {
                  pickgrp <- sample(strain.log[[i]], 1, prob=freq.log[[i]])
                }
                sampleWGS <- c(sampleWGS, pickgrp)
                if (0%in%sampleWGS) {
                  stop("Sampled zeroes")
                }
                sampleID <- c(sampleID, i)
                samplepick <- c(samplepick, j)
                sampletimes <- c(sampletimes, time)
              }
            }
            if (samp.schedule != "random" && rec.times[which(inf.ID==i)] > time + samp.freq) {
              sample.times[which(inf.ID==i)] <- time + samp.freq
            }
          }
        }
        
        # clean up libr etc.
        if (!full) {
          deleters <- NULL
          uniquestrains <- 0
          for (j in 1:length(totcurstrains)) {
            tottype <- 0
            for (k in cur.inf) {
              if (totcurstrains[j]%in%strain.log[[k]]) { # if strain is extant
                tottype <- tottype+1
              }
              if (sum(!strain.log[[k]]%in%totcurstrains)>0) {
                
                #cat("\nk=",k," strain.log[[k]]=",strain.log[[k]],"\n")
                stop("Deleted sequencee for observed sample")
              }
            }
            if (tottype>0) { # don't delete if still around
              uniquestrains <- uniquestrains+1
            } else if (tottype==0 && !totcurstrains[j]%in%sampleWGS) { # if not around AND not logged
              deleters <- c(deleters, j) # delete
            }
          }
          if (length(deleters)>0) {
            for (i in sort(deleters, decreasing=T)) {
              libr[[i]] <- NULL
              mut.nuc[[i]] <- NULL
            }
            deletegroup <- totcurstrains[deleters]
            totcurstrains <- totcurstrains[-deleters]
          }
        }
        if (length(cur.sus)==0) {
          eff.cur.inf <- inf.ID[which(sample.times>time)]
        }
        if (length(cur.sus)==0 && sum(sample.times>time)==0) {
          break
        }
      }
      if (tot.inf >= min.cases) {
        trigger <- TRUE
      } else {
        #cat("Insufficient number of infections outside! (min.cases=", min.cases, ")\n", sep="")
      }
      if (tot.inf.in >= min.cases.in) {
        trigger1 <- TRUE
      } else {
        #cat("Insufficient number of infections inside! (min.cases.in=", min.cases.in, ")\n", sep="")
      }
    }
    
    genomeID <- NULL
    for (i in 1:length(pID)) {
      genomeID <- c(genomeID, obs.strain[[which(pID%in%inf.ID[i])]])
    }
    
    if (full) {
      return(invisible(list(epidata=cbind(inf.ID, inf.times, rec.times, inf.source, genomeID), 
                            epidata_in=cbind(inf.ID.in, inf.times.in, rec.times.in, inf.source.in, introduction,genomeID.in),
                            sampledata=cbind(pID, sampleID, sampletimes),
                            
                            obs.freq=obs.freq, obs.strain=obs.strain,
                            
                            libr=libr, nuc=mut.nuc, librstrains=totcurstrains, endtime=time,
                            parameters=c(inf.rate, inf.rate.in, rem.rate, g.len, mut.rate.in,
                                         mut.rate.out, intr.rate)
      )))
    } else {
      return(invisible(list(epidata=cbind(inf.ID, inf.times, rec.times, inf.source), 
                            sampledata=cbind(sampleID, sampletimes, sampleWGS),
                            libr=libr, nuc=mut.nuc, librstrains=totcurstrains, endtime=time,
                            parameters=c(inf.rate, inf.rate.in, rem.rate, g.len, mut.rate.out,mut.rate.in, intr.rate)
      )))
    }
  }

