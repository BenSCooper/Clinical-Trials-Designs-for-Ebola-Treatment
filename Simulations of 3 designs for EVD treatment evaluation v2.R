  
  # This is based on version 3 of John Whitehead's MSA protocol for 
  # possible Ebola treatments. It simulates whole approach and compares outcomes with two other designs: a conventional RCT
  # and a sequential RCT.
  # This will produce figure 3 in the manuscript "Evaluating clinical trial designs for investigational treatments
  # of Ebola virus disease. " as a pdf and a number of other figures for component trials 
  # See manuscript for further details
  # To run "boundary3.csv" needs to be in the working directory.
  # Authors: Ben Cooper, Lisa White, Wirichada Pan-ngum, Nov 2014 
  
  prob.survival.if.treated<-seq(0.1,0.9,0.02)
  # uncomment 2 lines below to get  approximate numbers for table 1
   #  prob.survival.if.treated[29]<-0.667
    # prob.survival.if.treated[40]<-0.889
  patients.per.day<-5 # i.e. number who can be recruited to the trial per day
  
  nsim<-100000
  neff<-length(prob.survival.if.treated)
  out_size.phase2<-matrix(NA,nrow=nsim,ncol=neff)      #Design 3 MSA: number of samples adaptive phase II (1 arm)
  outcome.phase2<-matrix(NA,nrow=nsim,ncol=neff)      # Design 3 MSA:outcome of phase 2 (numbered 1 to 4)
  
  
  bound<-read.csv("boundary3.csv") # used in Design 3 MSA:adaptive phase II (1 arm)
  bound<-bound[-1,]
  bound[is.na(bound)]<- -1  # set NAs to -1 to avoid NA problems
  for (y in 1:neff){
    
    ########################################################################
    # Design 3 MSA: PHASE 2 SINGLE ARM
    ########################################################################
    maxboundary<-140 # maximum number of 14 day reports needed before a decision is made
    
    # create matrix of 14 day reports for each simulation
    survivors.day14<-matrix(rbinom(nsim*maxboundary,1,prob.survival.if.treated[y]), nrow=nsim)
    cum.survivors.day14<-t(apply(survivors.day14,1, "cumsum"))
    #sum(ifelse(1==bound[2,2:5],1,0),na.rm=T)
    a <- NULL
    num <- NULL
    total <- NULL
    btest<- 0
    for (i in 1: nsim){   
      # two ways to code this...method 1 is to evaluate each new person in trial using a while loop
      # (which maybe slow because of for loop)
      # 
      # method 2 much much faster (but gives same results as method 1) despite doing more calculations
      # Takes about 1 min to do 41*10,000 trial simulations on desktop
      
      method<-2
      
      if(method==1){
        n <- 1 # number of 14 day reports received
        stopped <-0 # set to 1 when trial stopped 
        while (n <=maxboundary && stopped==0){
          btest<-which(ifelse(cum.survivors.day14[i,n]==bound[n,2:5],1,0)==1)
          if(length(btest)==2) btest<-btest[1]
          if(sum(btest)>0) { # boundary hit
            boundary.hit<-btest
            n.at.boundary<-n
            stopped<-1
          } else {
            n<-n+1
          }     
        } 
        if(n>maxboundary) {
          # no boundary hit so we draw conclusion b (so set n.at.boundary to maxboundary and boundary.hit to 2) 
          n.at.boundary<-maxboundary
          boundary.hit=2
        }
      } # end of method 1
      
      
      if(method==2){
        # alternative methods of working out which boundary is hit and when
        # potentially faster as it avoids the while loop
        
        b1.hits<- which(cum.survivors.day14[i,]==bound$B1)
        b2.hits<- which(cum.survivors.day14[i,]==bound$B2)
        b3.hits<- which(cum.survivors.day14[i,]==bound$B3)
        b4.hits<- which(cum.survivors.day14[i,]==bound$B4)
        
        b1.min<- min(b1.hits,999) #no of recruits before hits boundary (999 if never hits it)
        b2.min<- min(b2.hits,999) #no of recruits before hits boundary (999 if never hits it)
        b3.min<- min(b3.hits,999) #no of recruits before hits boundary (999 if never hits it)
        b4.min<- min(b4.hits,999) #no of recruits before hits boundary (999 if never hits it)
        
        boundary.hit<-NA
        n.at.boundary<-b1.min
        if(n.at.boundary<999) boundary.hit<-1
        
        if(b2.min<n.at.boundary) {
          n.at.boundary <-b2.min
          boundary.hit<-2
        }
        if(b3.min<n.at.boundary) {
          n.at.boundary <-b3.min
          boundary.hit<-3
        }  
        
        if(b4.min<n.at.boundary) {
          n.at.boundary <-b4.min
          boundary.hit<-4
        }
        
        if(n.at.boundary==999){ #so no boundary hit then set boundary hit to 2 so decision is 2 
          n.at.boundary<-maxboundary
          boundary.hit<-2
        }
      }     # end of method 2
      
      out_size.phase2[i,y]<-n.at.boundary  # number recruited
      # three possible outcomes
      # outcome 1 if boundary 1 hit : outcome is draw conclusion a (use drug)
      # outcome 2 if boundary 2 or 3 hit : outcome is draw conclusion b (evaluate drug in RCT)
      # outcome 3 if boundary 4 hit : outcome is draw conclusion c (abandon drug)  
      if(boundary.hit==1) outcome<-1
      if(boundary.hit %in% 2:3) outcome<-2
      if(boundary.hit ==4) outcome<-3
      outcome.phase2[i,y]<-outcome
      
    }    # end for (i in 1: nsim)          
  } #  end for (y in 1:neff)
  
  # plot probs of different conclusions
  prob.1<-apply(outcome.phase2==1,2,mean) #  prob of adopting treatment
  prob.2<-apply(outcome.phase2==2,2,mean) #  prob of conducting placebo controlled trial
  prob.3<-apply(outcome.phase2==3,2,mean) #  prob of abandoning treatment
  
  #  some plots
  
  plot(prob.survival.if.treated,prob.1, type='l', col="green",xlim=c(0.2,.8),ylab="Prob of reaching conclusion")
  lines(prob.survival.if.treated,prob.2,col="orange")
  lines(prob.survival.if.treated,prob.3,col="red")
  
  med<-apply(out_size.phase2,2,median)
  avg<-apply(out_size.phase2,2,mean)
  avg<-apply(out_size.phase2,2,mean)
  pc90<-apply(out_size.phase2,2,quantile,0.9)
  pc10<-apply(out_size.phase2,2,quantile,0.1)
  
  plot(prob.survival.if.treated,med, type='l', col="blue",xlim=c(0.2,.8),ylim=c(0,120),ylab="Number 14 day results when trial stops")
  lines(prob.survival.if.treated,avg, col="black")
  lines(prob.survival.if.treated,pc90, col="red")
  lines(prob.survival.if.treated,pc10, col="red")
  
  ################ End of Design 3 MSA:phase II single arm simulation ###########################
  
  
  ##################################################
  # Design 3 MSA: Phase III single arm confirmatory trial
  ##################################################
  
  #The phase III confirmatory trial of Treatment A is designed to confirm that it is very effective with 
  #probability 0.900 if p = 0.8 
  #It will recruit a maximum of 132  patients, and will stop to abandon the experimental treatment if it is ever 
  #observed that that S ≤ 5.2425 + 0.7747n 
  
  
  out_size.phase3conf<-matrix(NA,nrow=nsim,ncol=neff)      # number of samples adaptive phase III confirmatory (1 arm)
  outcome.phase3conf<-matrix(NA,nrow=nsim,ncol=neff)      # outcome of phase 3 confirmatory (number 1 is cofirmed, 2 not confirmed)
  for (y in 1:neff){
    maxboundary<-132 # maximum number of 14 day reports needed before a decision is made (old version was 280)
    # create matrix of 14 day reports for each simulation
    survivors.day14.p3conf<-matrix(rbinom(nsim*maxboundary,1,prob.survival.if.treated[y]), nrow=nsim)
    cum.survivors.day14.p3conf<-t(apply(survivors.day14.p3conf,1, "cumsum"))
    
    a <- NULL
    num <- NULL
    total <- NULL
    btest<- 0
    for (i in 1: nsim){ 
      #    b.hits<- which(cum.survivors.day14.p3conf[i,]<= -8.22509+ 0.6779*1:maxboundary)
      b.hits<- which(cum.survivors.day14.p3conf[i,]<= -5.2425+ 0.7747*1:maxboundary)
      if(length(b.hits)>0) { # if boundary hit
        outcome.phase3conf[i,y]<-2
        out_size.phase3conf[i,y]<-b.hits[1]
      } else { # boundary not hit
        outcome.phase3conf[i,y]<-1
        out_size.phase3conf[i,y]<-maxboundary
      }
    } 
  }
  prob.phase3.confirmed<-apply(outcome.phase3conf==1,2,mean)
  
  prob.phase3.notconfirmed<-apply(outcome.phase3conf==2,2,mean)
  
  plot(prob.survival.if.treated,prob.phase3.notconfirmed, type='l', col="blue",xlim=c(0.2,.8),ylab="Prob of failing to confirm",xlab="Prob surviving to day 14")
  
  med<-apply(out_size.phase3conf,2,median)
  avg<-apply(out_size.phase3conf,2,mean)
  pc90<-apply(out_size.phase3conf,2,quantile,0.9)
  pc10<-apply(out_size.phase3conf,2,quantile,0.1)
  plot(prob.survival.if.treated,med, type='l', col="blue",xlim=c(0.2,.8),ylab="Final sample size",xlab="Prob surviving to day 14")
  points(prob.survival.if.treated,avg, type='l', col="black")
  points(prob.survival.if.treated,pc90, type='l', col="red")
  points(prob.survival.if.treated,pc10, type='l', col="red")
  
  ################ End of Design 3 MSA:phase III single arm confirmatory trial ###########################
  
  ##################################################
  # Design 2: Sequential RCT
  # also used as part of Design 3 MSA:Phase III placebo controlled trial 
  ##################################################
  
  
  prob.survival.if.not.treated<- 0.5 # needed for sims but not for analysis as there is a control group
  log.odds.ratio.for.treatment<-log((prob.survival.if.treated/(1-prob.survival.if.treated))/(prob.survival.if.not.treated/(1-prob.survival.if.not.treated)))
  neff<-length(prob.survival.if.treated)
  # interim analysis performed for every 25 new records (up to 20 in total)
  
  interim.analysis.times<-1:500 %% 25==0
  
  neff<-length(prob.survival.if.treated)
  outcome.phase3.rct<-matrix(NA,nrow=nsim,ncol=neff)
  out_size.phase3.rct<-matrix(NA,nrow=nsim,ncol=neff)      
  for (y in 1:neff){
    prob.surv.treated<-prob.survival.if.treated[y]
    trt<-matrix(rbinom(nsim*20*25,1,.5), nrow=nsim) # treatment allocation
    outcome<-matrix(rbinom(nsim*20*25,1, (trt==1)*prob.surv.treated  + (trt==0)*prob.survival.if.not.treated), nrow=nsim)    
    # outcome is 1 if survived, 0 if died  
    # now calc stats
    nA<-t(apply(trt,1, "cumsum")) #cumulative treated
    nC<-t(apply(trt==0,1, "cumsum")) #cumulative controls
    n<-nA+nC
    SA<-t(apply(trt==1 & outcome==1,1, "cumsum")) #cumulative treated and survivied
    SC<-t(apply(trt==0 & outcome==1,1, "cumsum")) #cumulative not treated and survivied
    Failure<-t(apply(outcome==0,1, "cumsum")) #cumulative died
    Success<-t(apply(outcome==1,1, "cumsum")) #cumulative survived
    
    Z=(nC*SA-nA*SC)/n
    SF<-Failure*Success
    V.1<-SF/(n^3) 
    V<-V.1*nA*nC
    
    stop.conclude.treatment.better<- Z>=6.39903 + 0.21049*V
    stop.conclude.treatment.not.better<- Z<= -6.39903 + 0.63148*V
    for (i in 1: nsim){ 
      b.trt.better<- which(stop.conclude.treatment.better[i,interim.analysis.times])
      b.trt.notbetter<- which(stop.conclude.treatment.not.better[i,interim.analysis.times])
      b.better.first.hit <-20
      b.nobetter.first.hit <-20      
      if(length(b.trt.better)>0) { # if boundary hit at interim time
        b.better.first.hit<-b.trt.better[1] #interim analysis number when boundary first hit
      }
      if(length(b.trt.notbetter)>0) { # if boundary hit at interim time
        b.nobetter.first.hit<-b.trt.notbetter[1] #interim analysis number when boundary first hit
      }  
      if(b.better.first.hit<b.nobetter.first.hit){  # outcome is that treatment is  better
        outcome.phase3.rct[i,y]<-1 
        out_size.phase3.rct[i,y]<-b.better.first.hit*25
      } else { # outcome is that treatment is not better
        outcome.phase3.rct[i,y]<-2 
        out_size.phase3.rct[i,y]<-b.nobetter.first.hit*25
      }   
      
    }
  } 
  
  prob.phase3rct.recommend.treatment<-apply(outcome.phase3.rct==1,2,mean)
  # Note that sample sizes below represent the number of patients with 14 day outcomes
  # in fact more will be recruited because we carry on recruiting until we reach a stpping rule for 14 day outcomess
  # So if we recruit 5 patients per day the actual number recruited will be sample size + 5*14
  mean.sample.size.phase3rct<-apply(out_size.phase3.rct,2,mean)
  prob.sample.size.phase3rct.lteq50<-apply(out_size.phase3.rct<=50,2,mean)
  prob.sample.size.phase3rct.lteq100<-apply(out_size.phase3.rct<=100,2,mean)
  prob.sample.size.phase3rct.lteq200<-apply(out_size.phase3.rct<=200,2,mean)
  prob.sample.size.phase3rct.lteq300<-apply(out_size.phase3.rct<=300,2,mean)
  prob.sample.size.phase3rct.lteq400<-apply(out_size.phase3.rct<=400,2,mean)
  prob.sample.size.phase3rct.lteq500<-apply(out_size.phase3.rct<=500,2,mean)
  
  ##################################################
  # End of Design 2: Sequential RCT
  ##################################################
  
  
  plot(log.odds.ratio.for.treatment,prob.phase3rct.recommend.treatment, ylab="Prob recommending treatment A",col="blue" ,type='l')
  plot(prob.survival.if.treated,prob.phase3rct.recommend.treatment, ylab="Prob recommending treatment A",col="blue" ,type='l')
  plot(prob.survival.if.treated,mean.sample.size.phase3rct, ylab="Expected sample size",col="black" ,type='l',xlim=c(0.2,.8),ylim=c(0,300),xaxs="i") 
  plot(prob.survival.if.treated, prob.sample.size.phase3rct.lteq50,type='l',ylim=c(0,1),xlim=c(0.2,0.8),xaxs="i",col="red")
  lines(prob.survival.if.treated, prob.sample.size.phase3rct.lteq100,col="blue")
  lines(prob.survival.if.treated, prob.sample.size.phase3rct.lteq200,col="green")
  lines(prob.survival.if.treated, prob.sample.size.phase3rct.lteq300,col="purple")
  lines(prob.survival.if.treated, prob.sample.size.phase3rct.lteq400,col="brown")
  lines(prob.survival.if.treated, prob.sample.size.phase3rct.lteq500,col="black")
  
  #  plot sample size against log OR 
  plot(log.odds.ratio.for.treatment,mean.sample.size.phase3rct, ylab="Expected sample size",col="black" ,type='l',xlim=c(-1,2),ylim=c(0,300),xaxs="i") 
  plot(log.odds.ratio.for.treatment, prob.sample.size.phase3rct.lteq50,type='l',ylim=c(0,1),xlim=c(0-1,2),xaxs="i",col="red")
  lines(log.odds.ratio.for.treatment, prob.sample.size.phase3rct.lteq100,col="blue")
  lines(log.odds.ratio.for.treatment, prob.sample.size.phase3rct.lteq200,col="green")
  lines(log.odds.ratio.for.treatment, prob.sample.size.phase3rct.lteq300,col="purple")
  lines(log.odds.ratio.for.treatment, prob.sample.size.phase3rct.lteq400,col="brown")
  lines(log.odds.ratio.for.treatment, prob.sample.size.phase3rct.lteq500,col="black")
  
  # expected recruits to stopping phase II trial of treatment A  if decision a, b or c and probs of decisions
  
  T1  <- apply(out_size.phase2,2,mean)
  T1a <- apply(out_size.phase2*(outcome.phase2==1), 2, sum)/apply(outcome.phase2==1, 2, sum) 
  T1b <- apply(out_size.phase2*(outcome.phase2==2), 2, sum)/apply(outcome.phase2==2, 2, sum) 
  T1c <- apply(out_size.phase2*(outcome.phase2==3), 2, sum)/apply(outcome.phase2==3, 2, sum) 
  P1a <- apply(outcome.phase2==1, 2, sum)/length(outcome.phase2[,1])  
  P1b <- apply(outcome.phase2==2, 2, sum)/length(outcome.phase2[,1])  
  P1c <- apply(outcome.phase2==3, 2, sum)/length(outcome.phase2[,1])  
  
  T2 <- apply(out_size.phase3conf,2,mean)
  T2a <- apply(out_size.phase3conf*(outcome.phase3conf==1),2,sum)/apply(outcome.phase3conf==1, 2, sum)
  T2b <- apply(out_size.phase3conf*(outcome.phase3conf==2),2,sum)/apply(outcome.phase3conf==2, 2, sum)
  P2a <- apply(outcome.phase3conf==1, 2, sum)/length(outcome.phase3conf[,1])  
  P2b <- apply(outcome.phase3conf==2, 2, sum)/length(outcome.phase3conf[,1])  
  
  T3  <- apply(out_size.phase3.rct,2,mean)
  T3a <- apply(out_size.phase3.rct*(outcome.phase3.rct==1),2,sum)/apply(outcome.phase3.rct==1, 2, sum)
  T3b <-apply(out_size.phase3.rct*(outcome.phase3.rct==2),2,sum)/apply(outcome.phase3.rct==2, 2, sum)
  P3a <- apply(outcome.phase3.rct==1, 2, sum)/length(outcome.phase3.rct[,1])  
  P3b <- apply(outcome.phase3.rct==2, 2, sum)/length(outcome.phase3.rct[,1])  
  
  T1a[is.nan(T1a)]<-0
  T1b[is.nan(T1b)]<-0
  T1c[is.nan(T1c)]<-0
  T2a[is.nan(T2a)]<-0
  T2b[is.nan(T2b)]<-0
  T3a[is.nan(T3a)]<-0
  T3b[is.nan(T3b)]<-0
  
  
  Prob.A.confirmed<-P1a*P2a + P1a*P2b*P3a + P1b*P3a
  Mean.recruits.to.confirmation.if.confirmed<-P1a*P2a*(T1a+T2a) +P1a*P2b*P3a *(T1a+T2b+T3a) + P1b*P3a*(T1b+T3a) /(P1a*P2a + P1a*P2b*P3a + P1b*P3a)
  Mean.recruits.to.outcome<-T1 + P1a*T2 +  P1a*P2b*T3 + P1b*T3 
  Mean.recruits.to.rollout.or.abandonment<- T1 + P1b*T3  # where rollout assumed if phaseII trial gives result a (ignoring fact that in some cases confirmation following rollout will not be achieved)
  Prob.RCT.performed<-P1a*P2b + P1b
  Expected.deaths.amongst.evaluated.patients<-T1*(1-prob.survival.if.treated) + P1a*T2*(1-prob.survival.if.treated) +  P1a*P2b*T3*(0.5*(1-prob.survival.if.treated) + 0.5*(1-prob.survival.if.not.treated)) + P1b*T3* (0.5*(1-prob.survival.if.treated) + 0.5*(1-prob.survival.if.not.treated))
  # note that above does not included patients who have been recruited but been included in the 14 day evaluation because the trial stopped
  
  plot(prob.survival.if.treated,Prob.A.confirmed, ylab="Prob recommending treatment A over all pathways",col="blue" ,type='l')
  plot(prob.survival.if.treated,Prob.RCT.performed, ylab="Prob RCT performed",col="blue" ,type='l')
  
  # Now for each value of prob.survival.if treated simulate N  trial networks by sampling from stored values and for each
  # record i) number of patients recruited i) number of deaths
  N<-nsim
  neff<-length(prob.survival.if.treated)
  # start with phase 2 one arm trial 
  # sample numper of patients in phase 2
  sim.num.patients.phase2<-apply(out_size.phase2,2,sample,replace=T)
  # sample numper of patients in phase 3 conf arm
  do.phase3.conf<-matrix(rbinom(neff*N, 1, P1a),nrow=N,byrow=T)
  sim.num.patients.phase3conf<-apply(out_size.phase3conf,2,sample,replace=T)*do.phase3.conf
  
  probofRCT<-P1a*P2b+P1b
  do.phase3.RCT<-matrix(rbinom(neff*N, 1, probofRCT),nrow=N,byrow=T)
  sim.num.patients.phase3RCT<-apply(out_size.phase3.rct,2,sample,replace=T)*do.phase3.RCT #  this is when the adaptive RCT is used after phase 2
  sim.num.patients.group.sequentialRCT<-apply(out_size.phase3.rct,2,sample,replace=T) #  this is when the adaptive RCT is used for everything
  
  sim.num.patients.total<-sim.num.patients.phase2 + sim.num.patients.phase3conf + sim.num.patients.phase3RCT
  
  sim.num.deaths.phase2<-matrix(rep(NA,neff*N),nrow=N)
  sim.num.deaths.phase3conf<-matrix(rep(NA,neff*N),nrow=N)
  sim.num.deaths.phase3RCT<-matrix(rep(NA,neff*N),nrow=N)
  sim.num.deaths.group.sequentialRCT<-matrix(rep(NA,neff*N),nrow=N)
  
  
  # calculate the number of deaths in the studies accounting for the fact that patients.per.day*14 are recruited in adaptive designs but not analysed (as the trials stops first)
  for(i in 1: N){
    # phase 2 one arm deaths 
    sim.num.deaths.phase2[i,]<-rbinom(neff,min(sim.num.patients.phase2[i,]+patients.per.day*14,140), (1-prob.survival.if.treated)) # sample number of deaths
    # phase 3 conf deaths
    zeros<-sim.num.patients.phase3conf[i,]==0
    sim.num.deaths.phase3conf[i,]<-rbinom(neff,pmin(sim.num.patients.phase3conf[i,]+patients.per.day*14,132), (1-prob.survival.if.treated)) # sample number of deaths
    sim.num.deaths.phase3conf[i,zeros]<-0 # reset to zero since no phase 3conf in these sims
    # phase 3 RCT deaths
    zeros<-sim.num.patients.phase3RCT[i,]==0
    prob.of.death.inRCT<-(1-prob.survival.if.treated)*0.5 + (1-prob.survival.if.not.treated)*0.5
    sim.num.deaths.phase3RCT[i,]<-rbinom(neff,pmin(sim.num.patients.phase3RCT[i,]+patients.per.day*14,500), prob.of.death.inRCT) # sample number of deaths in phase 3 RCT
    sim.num.deaths.phase3RCT[i,zeros]<-0 # reset to zero since no phase 3RCT in these sims
    sim.num.deaths.group.sequentialRCT[i,]<-rbinom(neff,sim.num.patients.group.sequentialRCT[i,]+patients.per.day*14, prob.of.death.inRCT) # sample number of deaths if adaptive RCT for everything 
    
  }
  
  num.patients.distribution<-apply(sim.num.patients.total,2, quantile,c(0.025,.05,.1,.5,.9,.95,.975))
  num.patients.distribution.RCTalone<-apply(sim.num.patients.group.sequentialRCT,2, quantile,c(0.025,.05,.1,.5,.9,.95,.975))
  
  
  sim.num.deaths.total.MSA<-sim.num.deaths.phase2 + sim.num.deaths.phase3conf + sim.num.deaths.phase3RCT
  num.deaths.mean.MSA<-apply(sim.num.deaths.total.MSA,2,mean)
  
  num.deaths.distribution.MSA<-apply(sim.num.deaths.total.MSA,2, quantile,c(0.025,.05,.1,.5,.9,.95,.975))
  num.deaths.mean.phase3RCT<-apply(sim.num.deaths.phase3RCT,2,mean)
  num.deaths.distribution.phase3RCT<-apply(sim.num.deaths.phase3RCT,2, quantile,c(0.025,.05,.1,.5,.9,.95,.975))
  
  num.deaths.mean.group.sequentialRCT<-apply(sim.num.deaths.group.sequentialRCT,2,mean)
  num.deaths.distribution.group.sequentialRCT<-apply(sim.num.deaths.group.sequentialRCT,2, quantile,c(0.025,.05,.1,.5,.9,.95,.975))
  
  
  
  # 
  # From STATA
  # . power twoproportions .5 .667, test(chi2) n(360)
  # 
  # Estimated power for a two-sample proportions test
  # Pearson's chi-squared test 
  # Ho: p2 = p1  versus  Ha: p2 != p1
  # 
  # Study parameters:
  # 
  #         alpha =    0.0500
  #             N =       360
  #   N per group =       180
  #         delta =    0.1670  (difference)
  #            p1 =    0.5000
  #            p2 =    0.6670
  # 
  # Estimated power:
  # 
  #         power =    0.8983
  
  ##################################################
  # Design 1: Conventional RCT (no interim)
  ##################################################
  
  
  fixedRCT.n<-360
  Mean.recruits.to.outcome.fixedRCT<- fixedRCT.n
  # nsim<-1000
  neff<-length(prob.survival.if.treated)
  prob.recommendingtreatment.fixedrct<-rep(0,neff)
  deaths.fixedrct<-matrix(rep(NA,neff*N),nrow=N)
  adopt<-rep(0,neff)
  for(i in 1:nsim){ 
    if(i %% 100 ==0) print(i)
    number.trt<-rbinom(neff,fixedRCT.n,.5)
    number.ctrl<-fixedRCT.n-number.trt
    number.trt.died<-rbinom(neff, number.trt, (1-prob.survival.if.treated))
    number.ctrl.died<-rbinom(neff, number.ctrl, (1-prob.survival.if.not.treated))
    number.trt.lived<-number.trt-number.trt.died
    number.ctrl.lived<-number.ctrl-number.ctrl.died
    deaths.fixedrct[i,]<-number.trt.died+number.ctrl.died
    for(j in 1:neff) {  
      conting.table<-matrix(c(number.trt.died[j],number.ctrl.died[j],number.trt.lived[j],number.ctrl.lived[j]),nrow=2)
      if(	chisq.test(conting.table,correct=FALSE)$p.value <0.05 & number.trt.died[j]/number.trt[j]<number.ctrl.died[j]/number.ctrl[j]){
          adopt[j]<- adopt[j]+1
      }
    }
  }
  prob.recommendingtreatment.fixedrct<-adopt/nsim
  expected.deaths.fixedRCT<-apply(deaths.fixedrct,2, mean)
  expected.deaths.fixedRCT.distribution=apply(deaths.fixedrct,2, quantile,c(0.025,.05,.1,.5,.9,.95,.975))
  
  ##################################################
  # End of design 1: Conventional RCT (no interim)
  ##################################################
  
  # Now calculate expected deaths resulting from over the 100 days following the start of the trial
  # Assume also that total deaths are reduced by treatment in the same proportion as 14 day deaths
  
  ndays<-100
  
  # Assuming exponential growth and a doubling time of 30 days (conservative assumption from NEJM 2014 paper from WHO team)
  cases.per.day<-c(200,rep(NA,ndays-1))
  for(i in 2:ndays) cases.per.day[i]<-(2^(1/30))*cases.per.day[i-1]
  # alternatively assume cases per day is constant - effectively implying an epidemic that is at the border of being controlled
  cases.per.day2<-rep(100,ndays) # at end of October there were about 400 cases per day..so this is very optimistic 
  proportion.of.deaths.occuring.within.14days<-pnorm(14,4.2,6.4) # based on numbers reported in NEJM by WHO Ebola response team
  prob.death.if.not.treated.total<-(1/proportion.of.deaths.occuring.within.14days)*(1-prob.survival.if.not.treated) # i.e. not just over 14 days
  prob.death.if.treated.total<-prob.death.if.not.treated.total*(1-prob.survival.if.treated)/(1-prob.survival.if.not.treated)
  cumulative.study.patients.by.day<-1:(1+ndays)*patients.per.day
  
  deaths.no.intervention<-sum(cases.per.day)*prob.death.if.not.treated.total
  deaths.no.intervention2<-sum(cases.per.day2)*prob.death.if.not.treated.total
  
  expected.num.treated.in.fixedRCT.study.by.day<-matrix(rep(0,ndays*neff),nrow=ndays)
  prob.rollout.infixedRCT.by.day<-matrix(rep(0,ndays*neff),nrow=ndays)
  prob.drugiscurrentlyrolled.by.day.infixedRCT<-matrix(rep(0,ndays*neff),nrow=ndays)
  trial.duration<-14+ceiling(fixedRCT.n/patients.per.day)
  proportion.of.cases.treated.following.rollout<-1
  non.study.cases<-cases.per.day-((1:ndays)<=trial.duration)*patients.per.day
  non.study.cases2<-cases.per.day2-((1:ndays)<=trial.duration)*patients.per.day
  
  for(i in 1:neff){
    expected.num.treated.in.fixedRCT.study.by.day[1:(trial.duration-14),i] <-rep(patients.per.day*0.5,(trial.duration-14))
  }
  prob.rollout.infixedRCT.by.day[trial.duration,]<-prob.recommendingtreatment.fixedrct
  for(i in trial.duration:ndays){
    prob.drugiscurrentlyrolled.by.day.infixedRCT[i,]<-prob.recommendingtreatment.fixedrct
  }
  
  study.deaths.fixedRCT<-0.5*fixedRCT.n*prob.death.if.not.treated.total + prob.death.if.treated.total*0.5*fixedRCT.n
  non.study.deaths.fixedRCT<-rep(NA,neff)
  non.study.deaths.fixedRCT2<-rep(NA,neff)
  
  for(i in 1:neff){
    treated.non.study.cases<-prob.drugiscurrentlyrolled.by.day.infixedRCT[,i]*non.study.cases*proportion.of.cases.treated.following.rollout
    untreated.non.study.cases.postrollout<-prob.drugiscurrentlyrolled.by.day.infixedRCT[,i]*non.study.cases*(1-proportion.of.cases.treated.following.rollout)
    untreated.non.study.cases.prerollout<-(1-prob.drugiscurrentlyrolled.by.day.infixedRCT[,i])*non.study.cases
    non.study.deaths<-treated.non.study.cases*prob.death.if.treated.total[i] +(untreated.non.study.cases.prerollout+untreated.non.study.cases.postrollout)*prob.death.if.not.treated.total
    non.study.deaths.fixedRCT[i]<-sum(non.study.deaths)
    
    treated.non.study.cases2<-prob.drugiscurrentlyrolled.by.day.infixedRCT[,i]*non.study.cases2*proportion.of.cases.treated.following.rollout
    untreated.non.study.cases.postrollout2<-prob.drugiscurrentlyrolled.by.day.infixedRCT[,i]*non.study.cases2*(1-proportion.of.cases.treated.following.rollout)
    untreated.non.study.cases.prerollout2<-(1-prob.drugiscurrentlyrolled.by.day.infixedRCT[,i])*non.study.cases2
    non.study.deaths2<-treated.non.study.cases2*prob.death.if.treated.total[i] +(untreated.non.study.cases.prerollout2+untreated.non.study.cases.postrollout2)*prob.death.if.not.treated.total
    non.study.deaths.fixedRCT2[i]<-sum(non.study.deaths2)
  }
  change.fixedRCT<-100*(deaths.no.intervention-(non.study.deaths.fixedRCT +study.deaths.fixedRCT))/deaths.no.intervention
  change.fixedRCT2<-100*(deaths.no.intervention2-(non.study.deaths.fixedRCT2 +study.deaths.fixedRCT))/deaths.no.intervention2
  
  #now calc expected number starting treatment by day if survival prob if treated is 60% or 80% (corresponds to elements 26 and 36 of prob.survival.if.treated)
  number.treated.fixedRCT.survprob60<-expected.num.treated.in.fixedRCT.study.by.day[,26]+prob.drugiscurrentlyrolled.by.day.infixedRCT[,26]*non.study.cases*proportion.of.cases.treated.following.rollout
  number.treated.fixedRCT.survprob80<-expected.num.treated.in.fixedRCT.study.by.day[,36]+prob.drugiscurrentlyrolled.by.day.infixedRCT[,36]*non.study.cases*proportion.of.cases.treated.following.rollout
  
  number.treated.fixedRCT2.survprob60<-expected.num.treated.in.fixedRCT.study.by.day[,26]+prob.drugiscurrentlyrolled.by.day.infixedRCT[,26]*non.study.cases2*proportion.of.cases.treated.following.rollout
  number.treated.fixedRCT2.survprob80<-expected.num.treated.in.fixedRCT.study.by.day[,36]+prob.drugiscurrentlyrolled.by.day.infixedRCT[,36]*non.study.cases2*proportion.of.cases.treated.following.rollout
  
  
  # Design 2: Sequential RCT with interim analsyes 
  expected.num.treated.in.adaptiveRCT.by.day<-matrix(rep(0,ndays*neff),nrow=ndays)
  prob.rollout.by.day.adaptiveRCT<-matrix(rep(0,ndays*neff),nrow=ndays)
  prob.drugiscurrentlyrolled.adaptiveRCT.by.day<-matrix(rep(0,ndays*neff),nrow=ndays)
  prob.finaldecision.adaptiveRCT.by.day<-matrix(rep(0,ndays*neff),nrow=ndays)
  
  success<-outcome.phase3.rct==1 
  fail<-outcome.phase3.rct==2 
  # out_size.phase3.rct  holds number of people in sequential RCT (Design 2)
  rollout.time<- 14*success+ ceiling(out_size.phase3.rct*success/patients.per.day) 
  # rollout.time==0   means never rolled out 
  decision.day<-14+(out_size.phase3.rct)/patients.per.day
  for(i in 1:neff){
    for(j in 1:nsim){
      expected.num.treated.in.adaptiveRCT.by.day[(1:ndays)<decision.day[j,i],i]<-  expected.num.treated.in.adaptiveRCT.by.day[(1:ndays)<decision.day[j,i],i]+ 0.5*patients.per.day
      prob.rollout.by.day.adaptiveRCT[(1:ndays)==ceiling(rollout.time[j,i]),i]<-1+prob.rollout.by.day.adaptiveRCT[(1:ndays)==ceiling(rollout.time[j,i]),i]
      prob.finaldecision.adaptiveRCT.by.day[(1:ndays)>=ceiling(decision.day[j,i]),i]<-1+prob.finaldecision.adaptiveRCT.by.day[(1:ndays)>=ceiling(decision.day[j,i]),i] 
      if(rollout.time[j,i]>0){
        prob.drugiscurrentlyrolled.adaptiveRCT.by.day[(1:ndays)>=rollout.time[j,i],i]<-1+prob.drugiscurrentlyrolled.adaptiveRCT.by.day[(1:ndays)>=rollout.time[j,i],i]
      }
    }
  }
  
  expected.num.treated.in.adaptiveRCT.by.day<-expected.num.treated.in.adaptiveRCT.by.day/nsim
  prob.rollout.by.day.adaptiveRCT<-prob.rollout.by.day.adaptiveRCT/nsim
  prob.drugiscurrentlyrolled.adaptiveRCT.by.day<-prob.drugiscurrentlyrolled.adaptiveRCT.by.day/nsim
  prob.finaldecision.adaptiveRCT.by.day<-prob.finaldecision.adaptiveRCT.by.day/nsim
  rollout.time.SRCT<-rollout.time
  
  # Now calculate expected deaths resulting from over the 100 days following the start of the trial
  
  #since on average same numbers treated and not treated in sequential (adaptive) RCT
  study.deaths.adaptiveRCT<-apply(expected.num.treated.in.adaptiveRCT.by.day,2,"sum")*(prob.death.if.treated.total+prob.death.if.not.treated.total)
  
  non.study.deaths.adaptiveRCT<-rep(NA,neff)
  non.study.deaths.adaptiveRCT2<-rep(NA,neff)
  
  for(i in 1:neff){
    non.study.cases <-cases.per.day-(1-prob.finaldecision.adaptiveRCT.by.day[,i])*patients.per.day
    treated.non.study.cases<-prob.drugiscurrentlyrolled.adaptiveRCT.by.day[,i]*non.study.cases*proportion.of.cases.treated.following.rollout
    untreated.non.study.cases.postrollout<-prob.drugiscurrentlyrolled.adaptiveRCT.by.day[,i]*non.study.cases*(1-proportion.of.cases.treated.following.rollout)
    untreated.non.study.cases.prerollout<-(1-prob.drugiscurrentlyrolled.adaptiveRCT.by.day[,i])*non.study.cases
    non.study.deaths<-treated.non.study.cases*prob.death.if.treated.total[i] +(untreated.non.study.cases.prerollout+untreated.non.study.cases.postrollout)*prob.death.if.not.treated.total
    non.study.deaths.adaptiveRCT[i]<-sum(non.study.deaths)
    
    non.study.cases2 <-cases.per.day2-(1-prob.finaldecision.adaptiveRCT.by.day[,i])*patients.per.day
    treated.non.study.cases2<-prob.drugiscurrentlyrolled.adaptiveRCT.by.day[,i]*non.study.cases2*proportion.of.cases.treated.following.rollout
    untreated.non.study.cases.postrollout2<-prob.drugiscurrentlyrolled.adaptiveRCT.by.day[,i]*non.study.cases2*(1-proportion.of.cases.treated.following.rollout)
    untreated.non.study.cases.prerollout2<-(1-prob.drugiscurrentlyrolled.adaptiveRCT.by.day[,i])*non.study.cases2
    non.study.deaths2<-treated.non.study.cases2*prob.death.if.treated.total[i] +(untreated.non.study.cases.prerollout2+untreated.non.study.cases.postrollout2)*prob.death.if.not.treated.total
    non.study.deaths.adaptiveRCT2[i]<-sum(non.study.deaths2)
  }
  change.adaptiveRCT<-100*(deaths.no.intervention-(non.study.deaths.adaptiveRCT +study.deaths.adaptiveRCT))/deaths.no.intervention
  change.adaptiveRCT2<-100*(deaths.no.intervention2-(non.study.deaths.adaptiveRCT2 +study.deaths.adaptiveRCT))/deaths.no.intervention2
  
  #now calc expected number treated by day if survival prob if treated is 60% [80%] (corresponds to element 26 [36]of prob.survival.if.treated)
  number.treated.adaptiveRCT.survprob60<-expected.num.treated.in.adaptiveRCT.by.day[,26]+prob.drugiscurrentlyrolled.adaptiveRCT.by.day[,26]*non.study.cases*proportion.of.cases.treated.following.rollout
  number.treated.adaptiveRCT.survprob80<-expected.num.treated.in.adaptiveRCT.by.day[,36]+prob.drugiscurrentlyrolled.adaptiveRCT.by.day[,36]*non.study.cases*proportion.of.cases.treated.following.rollout
  
  number.treated.adaptiveRCT2.survprob60<-expected.num.treated.in.adaptiveRCT.by.day[,26]+prob.drugiscurrentlyrolled.adaptiveRCT.by.day[,26]*non.study.cases2*proportion.of.cases.treated.following.rollout
  number.treated.adaptiveRCT2.survprob80<-expected.num.treated.in.adaptiveRCT.by.day[,36]+prob.drugiscurrentlyrolled.adaptiveRCT.by.day[,36]*non.study.cases2*proportion.of.cases.treated.following.rollout
  
  
  # MSA Design 3 
  # pathway 1 is decision a from phase 2 then single arm phase II confirmed
  # pathway 2 is decision a from phase 2 then single arm phase II not confirmed but A found better in RCT
  # pathway 3 is decision a from phase 2 then single arm phase II not confirmed and A not found better in RCT (path3 is almost never taken)
  # pathway 4 is decision b from phase 2 then  A found better in RCT
  # pathway 5 is decision b from phase 2 then  A not found better in RCT
  # pathway 6 is decision c from phase 2 
  
  
  path1<-outcome.phase2==1 & outcome.phase3conf==1 # outcome - recommend
  path2<-outcome.phase2==1 & outcome.phase3conf==2 & outcome.phase3.rct==1 # outcome - recommend
  path3<-outcome.phase2==1 & outcome.phase3conf==2 & outcome.phase3.rct==2 # outcome - do not recommend
  path4<-outcome.phase2==2 & outcome.phase3.rct==1   # outcome - recommend
  path5<-outcome.phase2==2 & outcome.phase3.rct==2  # outcome - do not recommend
  path6<-outcome.phase2==3 # outcome - do not recommend
  
  # out_size.phase2 holds number of people in phase 2 
  # out_size.phase3conf holds number of people in phase 3 confirmation
  
  expected.num.treated.in.MSA.study.by.day<-matrix(rep(0,ndays*neff),nrow=ndays)
  prob.rollout.by.day<-matrix(rep(0,ndays*neff),nrow=ndays)
  prob.endofrollout.by.day<-matrix(rep(0,ndays*neff),nrow=ndays) #only applies to pathyway 3
  prob.drugiscurrentlyrolled.by.day<-matrix(rep(0,ndays*neff),nrow=ndays)
  prob.finaldecision.by.day<-matrix(rep(0,ndays*neff),nrow=ndays)
  
  
  
  # If path 1 number treated before rollout is out_size.phase2 
  #- patients.per.day  are treated so time to rollout is 14+out_size.phase2*path1/patients.per.day
  rollout.time<-14+ceiling(out_size.phase2*path1/patients.per.day) # if 
  # rollout.time==0  means never rolled out 
  decision.day<-28+(out_size.phase2+out_size.phase3conf)*path1/patients.per.day
  for(i in 1:neff){
    for(j in 1:nsim){
      if(path1[j,i]){
        expected.num.treated.in.MSA.study.by.day[(1:ndays)<decision.day[j,i],i]<-expected.num.treated.in.MSA.study.by.day[(1:ndays)<decision.day[j,i],i]+patients.per.day
        prob.rollout.by.day[(1:ndays)==ceiling(rollout.time[j,i]),i]<-1+prob.rollout.by.day[(1:ndays)==ceiling(rollout.time[j,i]),i]
        prob.finaldecision.by.day[(1:ndays)>=ceiling(decision.day[j,i]),i]<-1+prob.finaldecision.by.day[(1:ndays)>=ceiling(decision.day[j,i]),i]
        if(rollout.time[j,i]>0){
          prob.drugiscurrentlyrolled.by.day[(1:ndays)>=rollout.time[j,i],i]<-1+prob.drugiscurrentlyrolled.by.day[(1:ndays)>=rollout.time[j,i],i]
        }
      } 
    }
  }
  
  
  # If path 2 number treated before rollout is out_size.phase2  
  # but rollout ends after results of phase3conf, but then happens again after confirmation from the RCT 
  # but number treated before final decision is out_size.phase2 + out_size.phase3conf+out_size.phase3.rct
  #- patients.per.day  are treated so time to rollout is out_size.phase2*path1/patients.per.day
  
  rollout.time[path2]<-(14+ceiling((out_size.phase2)*path2/patients.per.day ))[path2] 
  endof.rollout.time<-14*2+ceiling((out_size.phase2 + out_size.phase3conf)*path2/patients.per.day ) 
  rollout2.time<-14*3+ceiling((out_size.phase2 + out_size.phase3conf + out_size.phase3.rct)*path2/patients.per.day ) 
  rct.start.time<-14*2+ceiling((out_size.phase2  + out_size.phase3conf)*path2/patients.per.day)
  #  rollout.time==0 means never rolled out 
  
  decision.day[path2]<-rollout2.time[path2]
  for(i in 1:neff){
    for(j in 1:nsim){
      if(path2[j,i]){
        expected.num.treated.in.MSA.study.by.day[(1:ndays)<rct.start.time[j,i],i]<-expected.num.treated.in.MSA.study.by.day[(1:ndays)<rct.start.time[j,i],i]+patients.per.day
        #in RCT part only half patients per day are treated
        expected.num.treated.in.MSA.study.by.day[(1:ndays)>=rct.start.time[j,i] & (1:ndays)<rollout2.time[j,i],i]<-expected.num.treated.in.MSA.study.by.day[(1:ndays)>=rct.start.time[j,i] & (1:ndays)<rollout2.time[j,i],i]+0.5*patients.per.day
        prob.rollout.by.day[(1:ndays)==ceiling(rollout.time[j,i]),i]<-1+prob.rollout.by.day[(1:ndays)==ceiling(rollout.time[j,i]),i]
        prob.finaldecision.by.day[(1:ndays)>=ceiling(decision.day[j,i]),i]<-1+prob.finaldecision.by.day[(1:ndays)>=ceiling(decision.day[j,i]),i]
        if(rollout.time[j,i]>0){
          prob.drugiscurrentlyrolled.by.day[(1:ndays)>=rollout.time[j,i] & (1:ndays)<endof.rollout.time[j,i],i]<-1+prob.drugiscurrentlyrolled.by.day[(1:ndays)>=rollout.time[j,i] & (1:ndays)<endof.rollout.time[j,i],i]
          prob.drugiscurrentlyrolled.by.day[(1:ndays)>=rollout2.time[j,i],i]<-1+prob.drugiscurrentlyrolled.by.day[(1:ndays)>=rollout2.time[j,i],i]
        }
      } 
    }
  }
  
  
  # If path 3 number treated before rollout is out_size.phase2  
  # but rollout ends after results of phase3conf and doesn't happen  again after negative RCT result 
  # but number treated before final decision is out_size.phase2 + out_size.phase3conf+out_size.phase3.rct
  #- patients.per.day  are treated so time to rollout is out_size.phase2*path1/patients.per.day
  
  rollout.time[path3]<-(14+ceiling((out_size.phase2)*path3/patients.per.day ) )[path3]
  endof.rollout.time<-14*2+ceiling((out_size.phase2 + out_size.phase3conf)*path3/patients.per.day ) 
  rct.start.time<-14*2+ceiling((out_size.phase2  + out_size.phase3conf)*path3/patients.per.day)
  decision.day[path3]<-(14*3+ceiling((out_size.phase2  + out_size.phase3conf+out_size.phase3.rct)*path3/patients.per.day))[path3]
  for(i in 1:neff){
    for(j in 1:nsim){
      if(path3[j,i]){
        expected.num.treated.in.MSA.study.by.day[(1:ndays)<rct.start.time[j,i],i]<-expected.num.treated.in.MSA.study.by.day[(1:ndays)<rct.start.time[j,i],i]+patients.per.day
        #in RCT part only half patients per day are treated
        expected.num.treated.in.MSA.study.by.day[(1:ndays)>=rct.start.time[j,i] & (1:ndays)<decision.day[j,i],i]<-expected.num.treated.in.MSA.study.by.day[(1:ndays)>=rct.start.time[j,i] & (1:ndays)<decision.day[j,i],i]+0.5*patients.per.day
        prob.rollout.by.day[(1:ndays)==ceiling(rollout.time[j,i]),i]<-1+prob.rollout.by.day[(1:ndays)==ceiling(rollout.time[j,i]),i]
        prob.finaldecision.by.day[(1:ndays)>=ceiling(decision.day[j,i]),i]<-1+prob.finaldecision.by.day[(1:ndays)>=ceiling(decision.day[j,i]),i]
        if(rollout.time[j,i]>0){
          prob.drugiscurrentlyrolled.by.day[(1:ndays)>=rollout.time[j,i] & (1:ndays)<endof.rollout.time[j,i],i]<-1+prob.drugiscurrentlyrolled.by.day[(1:ndays)>=rollout.time[j,i] & (1:ndays)<endof.rollout.time[j,i],i]
        }
      } 
    }
  }
  
  # If path 4 number enrolled  before rollout is out_size.phase2  + out_size.phase3.rct
  # and number treated is out_size.phase2  + 0.5 out_size.phase3.rct
  
  rollout.time[path4]<-(14*2+(out_size.phase2 + out_size.phase3.rct)*path4/patients.per.day)[path4] # if 
  decision.day[path4]<-(14*2+(out_size.phase2+out_size.phase3.rct)*path4/patients.per.day)[path4]
  rct.start.time<-14+ceiling((out_size.phase2)*path4/patients.per.day)
  
  for(i in 1:neff){
    for(j in 1:nsim){
      if(path4[j,i]){
        expected.num.treated.in.MSA.study.by.day[(1:ndays)<rct.start.time[j,i],i]<-expected.num.treated.in.MSA.study.by.day[(1:ndays)<rct.start.time[j,i],i]+patients.per.day
        expected.num.treated.in.MSA.study.by.day[(1:ndays)>=rct.start.time[j,i] & (1:ndays)<decision.day[j,i],i]<-expected.num.treated.in.MSA.study.by.day[(1:ndays)>=rct.start.time[j,i] & (1:ndays)<decision.day[j,i],i]+0.5*patients.per.day
        prob.rollout.by.day[(1:ndays)==ceiling(rollout.time[j,i]),i]<-1+prob.rollout.by.day[(1:ndays)==ceiling(rollout.time[j,i]),i]
        prob.finaldecision.by.day[(1:ndays)>=ceiling(decision.day[j,i]),i]<-1+prob.finaldecision.by.day[(1:ndays)>=ceiling(decision.day[j,i]),i] 
        if(rollout.time[j,i]>0){
          prob.drugiscurrentlyrolled.by.day[(1:ndays)>=rollout.time[j,i],i]<-1+prob.drugiscurrentlyrolled.by.day[(1:ndays)>=rollout.time[j,i],i]
        }
      }
    }
  }
  
  
  # If path 5 number enrolled  before decision is out_size.phase2  + out_size.phase3.rct
  # and number treated is out_size.phase2  + 0.5 out_size.phase3.rct
  # There is no rollout.
  rollout.time[path5]<-0
  decision.day[path5]<-(14*2+(out_size.phase2+out_size.phase3.rct)*path5/patients.per.day)[path5]
  rct.start.time<-14+ceiling((out_size.phase2)*path5/patients.per.day)
  
  for(i in 1:neff){
    for(j in 1:nsim){
      if(path5[j,i]){
        expected.num.treated.in.MSA.study.by.day[(1:ndays)<rct.start.time[j,i],i]<-expected.num.treated.in.MSA.study.by.day[(1:ndays)<rct.start.time[j,i],i]+patients.per.day
        expected.num.treated.in.MSA.study.by.day[(1:ndays)>=rct.start.time[j,i] & (1:ndays)<decision.day[j,i],i]<-expected.num.treated.in.MSA.study.by.day[(1:ndays)>=rct.start.time[j,i] & (1:ndays)<decision.day[j,i],i]+0.5*patients.per.day
        prob.finaldecision.by.day[(1:ndays)>=ceiling(decision.day[j,i]),i]<-1+prob.finaldecision.by.day[(1:ndays)>=ceiling(decision.day[j,i]),i]
      } 
    }
  }
  
  # If path 6 number enrolled  before decision is out_size.phase2  
  # and number treated is out_size.phase2  
  # There is no rollout.
  
  rollout.time[path6]<-0
  decision.day[path6]<-(14+(out_size.phase2)*path6/patients.per.day)[path6]
  
  for(i in 1:neff){
    for(j in 1:nsim){
      if(path6[j,i]){
        expected.num.treated.in.MSA.study.by.day[(1:ndays)<decision.day[j,i],i]<-expected.num.treated.in.MSA.study.by.day[(1:ndays)<decision.day[j,i],i]+patients.per.day
        prob.finaldecision.by.day[(1:ndays)>=ceiling(decision.day[j,i]),i]<-1+prob.finaldecision.by.day[(1:ndays)>=ceiling(decision.day[j,i]),i]
      } 
    }
  }
  
  expected.num.treated.in.MSA.study.by.day<-expected.num.treated.in.MSA.study.by.day/nsim
  prob.rollout.by.day<-prob.rollout.by.day/nsim
  prob.drugiscurrentlyrolled.by.day<-prob.drugiscurrentlyrolled.by.day/nsim
  prob.finaldecision.by.day<-prob.finaldecision.by.day/nsim
  prob.study.ongoing<-1-prob.finaldecision.by.day
  expected.num.treated.in.MSA.study.<-apply(expected.num.treated.in.MSA.study.by.day,2,sum)
  expected.num.intotal.MSA.study<-apply(prob.study.ongoing*patients.per.day,2,sum)
  expected.num.untreated.in.MSA.study<-expected.num.intotal.MSA.study-expected.num.treated.in.MSA.study.
  mean.time.to.decision.MSA<-apply(decision.day,2,mean)
  mean.time.to.rollout.MSA<-apply(rollout.time,2,mean)
  time.of.rollout.or.rejection<-rollout.time
  time.of.rollout.or.rejection[rollout.time==0]<-decision.day[rollout.time==0]
  mean.time.to.rolloutorejection.MSA<- apply(time.of.rollout.or.rejection,2,mean)
  
  time.to.decision.MSA.distribution<-apply(decision.day,2, quantile,c(0.025,.05,.1,.5,.9,.95,.975))
  time.to.rollout.MSA.distribution<-apply(rollout.time,2, quantile,c(0.025,.05,.1,.5,.9,.95,.975))
  time.to.rolloutorejection.MSA.distribution<-apply(time.of.rollout.or.rejection,2, quantile,c(0.025,.05,.1,.5,.9,.95,.975))
  
  # Now calculate expected deaths resulting from over the 100 days following the start of the trial
  
  study.deaths.MSA.treated<-apply(expected.num.treated.in.MSA.study.by.day,2,"sum")*(prob.death.if.treated.total) 
  study.deaths.MSA.untreated<-expected.num.untreated.in.MSA.study*(prob.death.if.not.treated.total) 
  study.deaths.MSA<-study.deaths.MSA.treated+study.deaths.MSA.untreated
  non.study.deaths.MSA<-rep(NA,neff)
  non.study.deaths.MSA2<-rep(NA,neff)
  
  for(i in 1:neff){
    non.study.cases <-cases.per.day-(1-prob.finaldecision.by.day[,i])*patients.per.day
    treated.non.study.cases<-prob.drugiscurrentlyrolled.by.day[,i]*non.study.cases*proportion.of.cases.treated.following.rollout
    untreated.non.study.cases.postrollout<-prob.drugiscurrentlyrolled.by.day[,i]*non.study.cases*(1-proportion.of.cases.treated.following.rollout)
    untreated.non.study.cases.prerollout<-(1-prob.drugiscurrentlyrolled.by.day[,i])*non.study.cases
    non.study.deaths<-treated.non.study.cases*prob.death.if.treated.total[i] +(untreated.non.study.cases.prerollout+untreated.non.study.cases.postrollout)*prob.death.if.not.treated.total
    non.study.deaths.MSA[i]<-sum(non.study.deaths)
    
    non.study.cases2 <-cases.per.day2-(1-prob.finaldecision.by.day[,i])*patients.per.day
    treated.non.study.cases2<-prob.drugiscurrentlyrolled.by.day[,i]*non.study.cases2*proportion.of.cases.treated.following.rollout
    untreated.non.study.cases.postrollout2<-prob.drugiscurrentlyrolled.by.day[,i]*non.study.cases2*(1-proportion.of.cases.treated.following.rollout)
    untreated.non.study.cases.prerollout2<-(1-prob.drugiscurrentlyrolled.by.day[,i])*non.study.cases2
    non.study.deaths2<-treated.non.study.cases2*prob.death.if.treated.total[i] +(untreated.non.study.cases.prerollout2+untreated.non.study.cases.postrollout2)*prob.death.if.not.treated.total
    non.study.deaths.MSA2[i]<-sum(non.study.deaths2)
  }
  change.MSA<-100*(deaths.no.intervention-(non.study.deaths.MSA +study.deaths.MSA))/deaths.no.intervention
  change.MSA2<-100*(deaths.no.intervention2-(non.study.deaths.MSA2 +study.deaths.MSA))/deaths.no.intervention2
  
  #now calc expected number treated by day if survival prob if treated is 60% or 80% (corresponds to elements 26 and 36 of prob.survival.if.treated)
  number.treated.adaptive.survprob60<- expected.num.treated.in.MSA.study.by.day[,26]+prob.drugiscurrentlyrolled.by.day[,26]*non.study.cases*proportion.of.cases.treated.following.rollout
  number.treated.adaptive.survprob80<- expected.num.treated.in.MSA.study.by.day[,36]+prob.drugiscurrentlyrolled.by.day[,36]*non.study.cases*proportion.of.cases.treated.following.rollout
  
  number.treated.adaptive2.survprob60<- expected.num.treated.in.MSA.study.by.day[,26]+prob.drugiscurrentlyrolled.by.day[,26]*non.study.cases2*proportion.of.cases.treated.following.rollout
  number.treated.adaptive2.survprob80<- expected.num.treated.in.MSA.study.by.day[,36]+prob.drugiscurrentlyrolled.by.day[,36]*non.study.cases2*proportion.of.cases.treated.following.rollout
  
  change.MSA.100cases.per.day<-change.MSA2
  change.adaptiveRCT.100cases.per.day<-change.adaptiveRCT2
  change.fixedRCT.100cases.per.day<-change.fixedRCT2
  
  number.treated.fixedRCT.survprob80.100cases.per.day <- number.treated.fixedRCT2.survprob80
  number.treated.adaptiveRCT.survprob80.100cases.per.day <- number.treated.adaptiveRCT2.survprob80
  number.treated.adaptive.survprob80.100cases.per.day<-number.treated.adaptive2.survprob80
  
  #  Now plot the results
  figname<-paste("Three evaluation strategies compared (",100*prob.survival.if.not.treated," per cent untreated mortailty).pdf",sep="")
  pdf(figname, height=8, width=7, pointsize=12)
  par(mfcol=c(3,2),mar=c(4,4,1,4),oma=c(1,1,1,1),cex.lab=1)
  # prob of recommending treatment
  plot(prob.survival.if.treated,Prob.A.confirmed, xlab='Probability of surviving to day 14 if treated', xaxt='n',ylab="Probability of adopting treatment",col="blue" ,type='l',bty='n')
  #rect(0.1,0,0.5,1,col=rgb(50,50,50,50,maxColorValue=255),border=NA)
  
  lines(prob.survival.if.treated,Prob.A.confirmed,col="blue")
  lines(prob.survival.if.treated,prob.phase3rct.recommend.treatment,col="green4")  # group sequential RCT
  #lines(prob.survival.if.treated,prob.recommendingtreatment.fixedrct,col="red",lty=2,lwd=2)
  points(prob.survival.if.treated,prob.recommendingtreatment.fixedrct,col="red",pch=20,cex=.5)
  
  axis(side=1,at=(1:9)/10)
  #legend(0.08,0.85,legend=c("Conventional RCT", "Sequential RCT", "Multi-stage approach"),col=c("red","green4","blue"),bty='n',lty=1,seg.len = 1)
  legend(0.08,0.85,legend=c("Conventional RCT", "Sequential RCT", "Multi-stage \napproach"),col=c("red","green4","blue"),bty='n',pch=c(20,NA,NA),lty=c(NA,1,1),seg.len = 1,cex=1)
  
  #abline(v=0.333,col=grey(.30),lwd=1.5)
  
  text(0.1,1,"A")
  
  #below assumes 5 patients can be recruited per day throughout entire period
  patients.per.day<-5
  
  #Design 3: MSA
  
  plot(prob.survival.if.treated,mean.time.to.decision.MSA, xlab="Probability of surviving to day 14 if treated", ylab="Time to roll-out or rejection (days)",lty=2,col="blue" ,type='l',xlim=c(0.1,0.9),ylim=c(0,125),xaxt="n",bty='n')
  lines(prob.survival.if.treated,mean.time.to.rolloutorejection.MSA,lty=1,col="blue")
  
  lower<-time.to.rolloutorejection.MSA.distribution[2,]# 5th percentile
  upper<-time.to.rolloutorejection.MSA.distribution[6,] # 95th percentile
  polygon(c(prob.survival.if.treated, rev(prob.survival.if.treated)), c(lower, rev(upper)),col = rgb(0,0,100,50,maxColorValue=255), border = NA)
  
  
  
  
  #plot(prob.survival.if.treated,14+Mean.recruits.to.outcome/patients.per.day, xlab="", ylab="Time to roll-out or rejection (days)",lty=2,col="blue" ,type='l',xlim=c(0.1,0.9),ylim=c(0,500/5),xaxt="n",bty='n')
  #lines(prob.survival.if.treated,14+Mean.recruits.to.rollout.or.abandonment/patients.per.day,lty=1,col="blue")
  text(0.1,120,"B")
  
  #lower<-14+num.patients.distribution[2,]/patients.per.day # 5th percentile
  #upper<-14+num.patients.distribution[6,]/patients.per.day # 95th percentile
  #polygon(c(prob.survival.if.treated, rev(prob.survival.if.treated)), c(lower, rev(upper)),col = rgb(0,0,100,50,maxColorValue=255), border = NA)
  #lines(prob.survival.if.treated,14+Mean.recruits.to.outcome/patients.per.day,lty=2,col="blue")
  #lines(prob.survival.if.treated,14+Mean.recruits.to.rollout.or.abandonment/patients.per.day,lty=1,col="blue")
  
  # Design 1: conventional RCT
  
  # axis(side=4, at=14+c(0,100,200,300,400,500)/patients.per.day, labels=c(0,100,200,300,400,500))
  axis(side=1,at=(1:9)/10)
  #mtext(side=4,"Number of patients",line=2,cex=.7)  # because of the e variable number 14 day periods after last recruitment  (because of variable number of study components in MSA)
  # we don't have a 1-1 mapping from number of paitents to time in the MSA design.  We do in the other two designs.
  # lines(prob.survival.if.treated,14+rep(360,neff)/patients.per.day,col="red")  # fixed RCT
  points(prob.survival.if.treated,14+rep(360,neff)/patients.per.day,col="red",pch=20,cex=.5)  # fixed RCT
  
  # Design 2: group sequential RCT
  lower<-14+num.patients.distribution.RCTalone[2,]/patients.per.day # 5th percentile
  upper<-14+num.patients.distribution.RCTalone[6,]/patients.per.day # 95th percentile
  polygon(c(prob.survival.if.treated, rev(prob.survival.if.treated)), c(lower, rev(upper)),col = rgb(0,100,0,50,maxColorValue=255), border = NA)
  lines(prob.survival.if.treated,14+mean.sample.size.phase3rct/patients.per.day,col="green4")
  
  
  
  par(mar=c(4,4,1,4))
  
  plot(prob.survival.if.treated,Prob.RCT.performed, ylab="Probability RCT is performed",col="blue" ,type='l', xlab='Probability of surviving to day 14 if treated', xaxt='n',bty='n')
  axis(side=1,at=(1:9)/10)
  text(0.1,1,"C")
  
  #Design 3: MSA
  plot(prob.survival.if.treated,num.deaths.mean.MSA, ylab="Deaths in study",col="blue" ,type='l',ylim=c(0,300), xlab="Probability of surviving to day 14 if treated", xaxt='n',bty='n')
  lower<-num.deaths.distribution.MSA[2,] # 5th percentile
  upper<-num.deaths.distribution.MSA[6,] # 95th percentile
  polygon(c(prob.survival.if.treated, rev(prob.survival.if.treated)), c(lower, rev(upper)),col = rgb(0,0,100,50,maxColorValue=255), border = NA)
  axis(side=1,at=(1:9)/10)
  text(0.9,300,"D")
  # Design 1: conventional RCT
  #  lines(prob.survival.if.treated,expected.deaths.fixedRCT,col="red")
  points(prob.survival.if.treated,expected.deaths.fixedRCT,col="red",pch=20,cex=.5)
  lower<-expected.deaths.fixedRCT.distribution[2,] # 5th percentile
  upper<-expected.deaths.fixedRCT.distribution[6,] # 95th percentile
  polygon(c(prob.survival.if.treated, rev(prob.survival.if.treated)), c(lower, rev(upper)),col = rgb(100,0,0,50,maxColorValue=255), border = NA)
  
  
  lines(prob.survival.if.treated,num.deaths.mean.MSA, col="blue")
  # Design 2: group sequential RCT
  lines(prob.survival.if.treated,num.deaths.mean.group.sequentialRCT, col="green4")
  lower<-num.deaths.distribution.group.sequentialRCT[2,] # 5th percentile
  upper<-num.deaths.distribution.group.sequentialRCT[6,] # 95th percentile
  polygon(c(prob.survival.if.treated, rev(prob.survival.if.treated)), c(lower, rev(upper)),col = rgb(0,100,0,50,maxColorValue=255), border = NA)
  
  #Design 3: MSA
  plot(prob.survival.if.treated,change.MSA.100cases.per.day, ylab="Mortality reduction (%)",col="blue" ,type='l',lty=1,ylim=c(-1,100), xlab='Probability of surviving to day 14 if treated', xaxt='n',bty='n')
  # Design 2: group sequential RCT
  lines(prob.survival.if.treated,change.adaptiveRCT.100cases.per.day, col="green4",lty=1 )
  # Design 1: conventional RCT
  #lines(prob.survival.if.treated,change.fixedRCT.100cases.per.day, col="red",lty=1 )
  points(prob.survival.if.treated,change.fixedRCT.100cases.per.day, col="red",pch=20,cex=.5)
  
  #Design 3: MSA
  lines(prob.survival.if.treated,change.MSA, col="blue" ,lty=2)
  # Design 2: group sequential RCT
  lines(prob.survival.if.treated,change.adaptiveRCT, col="green4",lty=2 )
  # Design 1: conventional RCT
  #lines(prob.survival.if.treated,change.fixedRCT, col="red",lty=2 )
  points(prob.survival.if.treated,change.fixedRCT, col="red",pch=1,cex=.5)
  
  axis(side=1,at=(1:9)/10)
  text(0.9,100,"E")
  
  # Design 2: group sequential RCT
  plot(1:100,number.treated.adaptiveRCT.survprob80,log="y",col="green4",xlab="Day" ,ylim=c(1,20000), ylab="Patients starting treatment per day",type="l",bty='n',lty=2)
  #Design 3: MSA
  lines(number.treated.adaptive.survprob80,col="blue",lty=2)
  # Design 1: conventional RCT
  #lines(number.treated.fixedRCT.survprob80,col="red",lty=2,lwd = 1.2)
  points(number.treated.fixedRCT.survprob80,col="red",pch=1,cex=0.5)
  
  # Design 2: group sequential RCT
  lines(number.treated.adaptiveRCT.survprob80.100cases.per.day,col="green4",lty=1)
  #Design 3: MSA
  lines(number.treated.adaptive.survprob80.100cases.per.day,col="blue",lty=1)
  # Design 1: conventional RCT
  #lines(number.treated.fixedRCT.survprob80.100cases.per.day,col="red",lty=1,lwd = 1.2)
  points(number.treated.fixedRCT.survprob80.100cases.per.day,col="red",pch=20,cex=.5)
  
  text(100,20000,"F")
  
  
  dev.off()
  
  #  numbers needed for table one approximations (actual numbers are based on analytial results - but we can confirm these results by simulation) 
  prob.survival.if.not.treated
  prob.survival.if.treated
  #MSA
  Prob.A.confirmed 
  Mean.recruits.to.outcome
  #SRCT
  prob.phase3rct.recommend.treatment
  mean.sample.size.phase3rct
  
  #fixed RCT
  prob.recommendingtreatment.fixedrct
  # 360 patients
  
