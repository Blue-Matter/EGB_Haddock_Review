# Sim testing functions

doSimTest<-function(dir,useinitsmod=F,nsim=4,seed=0,catch=1,ind=1,ecov=1,NAA_p=0,M_p=0,sel_p=0,ecov_p=0,q_p=0,inits=NA){

  fit<-readRDS(paste0(dir,"/Fit.rda"))

  #if(!is.na(dirmod))mod<-readRDS(paste0(dirmod,"/Fit.rda"))
  #There are elements input$data$simulate_data and input$data$simulate_state, input$data$simulate_period which are all 0/1 indicator vectors.
  #simulate_data is length 3: catch, indices, Ecov
  #simulate_state is length 5 for process errors in NAA, M, sel, Ecov, q (only when those process errors are being used)
  #simulate_period is 2 for the data period and projection period which affects both simulate_data and simulate_state.

  fit$env$data$simulate_state = c(NAA_p, M_p, sel_p, ecov_p, q_p)# keep process errors as estimated
  fit$env$data$simulate_data = c(catch, ind, ecov)

  indata<-list()
  for(i in 1:nsim){

    set.seed(i+seed)
    simdata = fit$simulate(complete = TRUE)

    if(useinitsmod){ # do you want to use the model of inits for the sim testing (different models for sim and test will need different inits)
      siminput = inits$input # estimates using the model in the inits object inits
    }else{
      siminput = fit$input # estimates using same model as the simulator
      if(!is.na(inits[1])) siminput$par = inits$par  # use inits object for initialization
    }

    nams<-c("agg_catch","catch_paa","agg_indices","index_paa") # select simulated data to be copied
    # nams<-names(siminput$data)
    for(j in 1:length(nams)){ 
      siminput$data[[nams[j]]]<-simdata[[nams[j]]]
    }

    indata[[i]] = siminput
    
  }

  # will provide (by default) simulated data (using estimated/input variances). And a fit to the simulated data can be accomplished by:

  WHAMcluster()
  Fits <- sfLapply(1:nsim, fit_wham_parallel, datlist = indata, do.retro=F)
  # sapply(1,fit_wham_parallel, datlist = indata, do.retro=F)
  sfStop()
  Fitsmall <- lapply(Fits,extract_WHAM)
  simout<-list(extract_WHAM(fit),Fitsmall)
  simout

}

plotSimSum<-function(simout){
  
  par(mfcol=c(2,3),mai=c(0.4,0.7,0.1,0.01),omi=c(0.3,0.1,0.1,0.01))
  conv<-as.vector(lapply(simout[[2]],function(x)x$conv))==0
  LHFs<-sapply(simout[[2]],function(x)x$nll[1])
  conv<-conv&!(LHFs>(mean(LHFs)+2*sd(LHFs)))
  
  SSB_T<-simout[[1]]$SSB/1000
  SSB_E<-sapply(simout[[2]],function(x)x$SSB/1000)[,conv]
  years<-simout[[1]]$years
  simcomp2(VT=SSB_T,VE=SSB_E,lab="SSB (kt)",xlab="",xtick=years)
  legend('topright',legend=paste0("n = ",sum(conv)), bty="n")
  mtext("Year",1,outer=F,line=2.7)
  
  F_T<-simout[[1]]$Fbar
  F_E<-sapply(simout[[2]],function(x)x$Fbar)[,conv]
  simcomp2(VT=F_T,VE=F_E,lab="Fbar",xlab="",xtick=years)
  mtext("Year",1,outer=F,line=2.7)
  
  ny<-nrow(simout[[1]]$NAA)
  N_T<-simout[[1]]$NAA[ny,]
  ages<-1:ncol(simout[[1]]$NAA)
  N_E<-sapply(simout[[2]],function(x)x$NAA[ny,])[,conv]
  simcomp2(VT=N_T,VE=N_E,lab="Numbers 2019",xlab="",xtick=ages)
  mtext("Age",1,outer=F,line=2.7)
  
}

plotSSB<-function(simout){
  conv<-as.vector(lapply(simout[[2]],function(x)x$conv))==0
  SSB_T<-simout[[1]]$SSB
  years<-simout[[1]]$years
  SSB_E<-sapply(simout[[2]],function(x)x$SSB)[,conv]
  simcomp(VT=SSB_T,VE=SSB_E,lab="SSB",xlab="Year",xtick=years)
}

plotfN<-function(simout){
  conv<-as.vector(lapply(simout[[2]],function(x)x$conv))==0
  ny<-nrow(simout[[1]]$NAA)
  N_T<-simout[[1]]$NAA[ny,]
  ages<-1:ncol(simout[[1]]$NAA)
  N_E<-sapply(simout[[2]],function(x)x$NAA[ny,])[,conv]
  simcomp(VT=N_T,VE=N_E,lab="Numbers 2019",xlab="Age",xtick=ages)
}

plotRec<-function(simout){
  conv<-as.vector(lapply(simout[[2]],function(x)x$conv))==0
  ny<-nrow(simout[[1]]$NAA)
  years<-simout[[1]]$years
  N_T<-simout[[1]]$NAA[,1]
  N_E<-sapply(simout[[2]],function(x)x$NAA[,1])[,conv]
  simcomp(VT=N_T,VE=N_E,lab="Recruitment",xlab="Year",xtick=years)
}

plotVB<-function(simout){
  conv<-as.vector(lapply(simout[[2]],function(x)x$conv))==0
  VB_T<-getVB(simout[[1]])
 years<-simout[[1]]$years
 VB_E<-sapply(simout[[2]],getVB)[,conv]
 simcomp(VT=VB_T,VE=VB_E,lab="Vuln. Bio.",xlab="Year",xtick=years)
}

plotFXSPR<-function(simout){
  conv<-as.vector(lapply(simout[[2]],function(x)x$conv))==0
  FXSPR_T<-simout[[1]]$FXSPR
  years<-simout[[1]]$years
  FXSPR_E<-sapply(simout[[2]],function(x)x$FXSPR)[,conv]
  simcomp(VT=FXSPR_T,VE=FXSPR_E,lab="FXSPR",xlab="Year",xtick=years)
}

getsel<-function(x)x$FAA[,1,]/apply(x$FAA[,1,],1,max)
getVB<-function(x)apply(x$BAA*getsel(x),1,sum)

plotMQ<-function(simout){
  conv<-as.vector(lapply(simout[[2]],function(x)x$conv))==0
  
  years<-simout[[1]]$years
  ny<-length(years)
  FXSPR_T<-simout[[1]]$FXSPR[ny]
  FXSPR_E<-sapply(simout[[2]],function(x)x$FXSPR[ny])[conv]
  VB_T<-getVB(simout[[1]])[ny]/1000
  VB_E<-sapply(simout[[2]],function(x){getVB(x)[ny]})[conv]/1000
  MR_T<-FXSPR_T*VB_T
  MR_E<-FXSPR_E*VB_E

  VT=FXSPR_T
  VE=FXSPR_E
  lab="FXSPR"
  par(mfcol=c(2,3),mai=c(0.05,0.65,0.05,0.1),omi=rep(0,4))
  qcomp(FXSPR_T,FXSPR_E,"FXSPR")
  qcomp(VB_T,VB_E,"Vuln. Biom. (kt)")
  qcomp(MR_T,MR_E,"FXSPR x VB (kt)")
}

doaxis2<-function(){
  axis(1,c(-1E10,1E10),rep("",2))
  axis(2)
  axis(3,c(-1E10,1E10),rep("",2))
  axis(4,c(-1E10,1E10),rep("",2))
}

qcomp<-function(VT,VE,lab="",xlab="",qupper=0.999,ecol="#ff000050"){

  rng<-c(0,max(c(VT,VE)))
  plot(VT,ylim=rng,pch=19,cex=1.1,ylab=lab,yaxs="i",axes=F)
  doaxis2()
  qs<-matrix(quantile(VE,c(0.05,0.25,0.5,0.75,0.95)),ncol=1)
  qplot(qs,1,ecol=ecol)
  points(VT,pch=19,cex=1.1)

  err<-(VE-VT)/VT*100


  qs<-matrix(quantile(err,p=c(0.05,0.25,0.5,0.75,0.95)),ncol=1)
  rng<-range(0,qs)
  plot(0,ylim=rng,pch=19,cex=1.1,ylab=paste(lab,"(Est-Sim)/Sim"),axes=F)
  doaxis2()
  abline(h=(-100:100)*10,col='light grey')
  qplot(qs,1,ecol=ecol)
  points(0,pch=19,cex=1.1)

}

simcomp<-function(VT,VE,lab="",xlab="",qupper=0.999,ecol="#ff000050",xtick){

  par(mfrow=c(2,2),mai=c(0.4,0.7,0.1,0.01),omi=c(0.6,0.1,0.5,0.01))
  VT[VT==Inf]<-NA
  VE[VE==Inf]<-NA
  rng<-c(0,quantile(c(VT,VE),qupper,na.rm=T))
  plot(xtick,VT,ylim=rng,type="l",lwd=2,ylab=lab,yaxs="i")
  matplot(xtick,VE,type="l",col=ecol,add=T,lty=1,lwd=2)
  lines(xtick,VT,lwd=2)

  qs<-apply(VE,1,quantile,p=c(0.05,0.25,0.5,0.75,0.95),na.rm=T)
  rng<-c(0,max(qs,na.rm=T))
  plot(xtick,VT,ylim=rng,type="l",lwd=2,ylab=lab,yaxs="i")
  qplot(qs,xtick,ecol=ecol)
  lines(xtick,VT,lwd=2)

  err<-(VE-VT)/VT*100
  rng<-quantile(err,c(1-qupper,qupper),na.rm=T)
  plot(xtick,rep(0,length(xtick)),ylim=rng,type="l",lwd=2,ylab=paste(lab,"(Est-Sim)/Sim"))
  matplot(xtick,err,type="l",col=ecol,add=T,lty=1,lwd=1)

  qs<-apply(err,1,quantile,p=c(0.05,0.25,0.5,0.75,0.95),na.rm=T)
  rng<-range(qs)
  plot(xtick,rep(0,length(xtick)),ylim=rng,type="l",lwd=2,ylab=paste(lab,"(Est-Sim)/Sim"))
  abline(h=(-100:100)*10,col='light grey')
  abline(h=(-100:100)*50,col='dark grey')
  qplot(qs,xtick,ecol=ecol)

  mtext(xlab,1,outer=T,line=0.7)

}

simcomp2<-function(VT,VE,lab="",xlab="",qupper=0.999,ecol="#ff000050",xtick,zlab=""){
  
  #par(mfrow=c(2,2),mai=c(0.4,0.7,0.1,0.01),omi=c(0.6,0.1,0.5,0.01))
  VT[VT==Inf]<-NA
  VE[VE==Inf]<-NA
  rng<-c(0,quantile(c(VT,VE),qupper,na.rm=T))
  #plot(xtick,VT,ylim=rng,type="l",lwd=2,ylab=lab,yaxs="i")
  #matplot(xtick,VE,type="l",col=ecol,add=T,lty=1,lwd=2)
  #lines(xtick,VT,lwd=2)
  
  qs<-apply(VE,1,quantile,p=c(0.05,0.25,0.5,0.75,0.95),na.rm=T)
  rng<-c(0,max(qs,na.rm=T))
  plot(xtick,VT,ylim=rng,type="l",lwd=2,ylab=lab,yaxs="i")
  qplot(qs,xtick,ecol=ecol)
  lines(xtick,VT,lwd=2)
  mtext(zlab,3,0.4)
  
  err<-(VE-VT)/VT*100
  rng<-quantile(err,c(1-qupper,qupper),na.rm=T)
  #plot(xtick,rep(0,length(xtick)),ylim=rng,type="l",lwd=2,ylab=paste(lab,"(Est-Sim)/Sim"))
  #matplot(xtick,err,type="l",col=ecol,add=T,lty=1,lwd=1)
  
  qs<-apply(err,1,quantile,p=c(0.05,0.25,0.5,0.75,0.95),na.rm=T)
  rng<-range(qs,na.rm=T)
  if(rng[1]<(-150))rng[1]<-(-150)
  if(rng[2]>(150))rng[2]<-150
  plot(xtick,rep(0,length(xtick)),ylim=rng,type="l",lwd=2,ylab=paste(lab,"(Est-Sim)/Sim"))
  abline(h=(-100:100)*10,col='light grey')
  abline(h=(-100:100)*50,col='dark grey')
  qplot(qs,xtick,ecol=ecol)
  
  mtext(xlab,1,outer=T,line=0.7)
  
}


qplot<-function(qs,xs,ecol){

  for(i in 1:ncol(qs)){
    lines(rep(xs[i],2),qs[c(1,5),i],col=ecol,lwd=2)
    lines(rep(xs[i],2),qs[c(2,4),i],col=ecol,lwd=4)
  }
  points(xs,qs[3,],col=ecol,pch=19,cex=1.1)

}
