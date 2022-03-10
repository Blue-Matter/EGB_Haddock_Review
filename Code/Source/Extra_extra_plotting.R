# new plots

Recent_Ind<-function(listo,YLabs=NULL,XLabs=c("NMFS Spr","NMFS Fall","DFO Spr")){
  
  cols<-c("red","green","blue")
  labcex=0.85
  ni<-ncol(listo[[1]]$agg_indices)
  
  get_obs<-function(x){
    ny<-nrow(x$agg_indices)
    yind<-ny-(9:0)
    x$agg_indices[yind,]
  }
  
  get_pred<-function(x){
    ny<-nrow(x$pred_indices)
    yind<-ny-(9:0)
    x$pred_indices[yind,]
  }
  
  obs<-lapply(listo,get_obs)
  pred<-lapply(listo,get_pred)
  ny<-length(listo[[1]]$years)
  yrs<-listo[[1]]$years[ny-(9:0)]
  xlimALL<-range(yrs)
  
  obsmax<-t(sapply(obs,function(x)apply(x,2,max,na.rm=T)))
  predmax<-t(sapply(pred,function(x)apply(x,2,max,na.rm=T)))
  maxi<-apply(rbind(obsmax,predmax),2,max)
  
  ncomp<-length(listo)
  if(is.null(YLabs))YLabs = paste("Fit",1:ncomp)
  if(is.null(XLabs))XLabs = c(paste("Index",1:ni))
   
  par(mfrow=c(ncomp,ni+1),mai=c(0.1,0.25,0.05,0.05),omi=c(0.5,0.6,0.25,0.01))
  for(i in 1:ncomp){
     for(j in 1:ni){
      yline<-pretty(seq(0,maxi[j]*2,length.out=18))
      obso<-obs[[i]][,j]
      obso[obso==0]<-NA
      plot(yrs,obso,ylim=c(0,maxi[j]*1.02),xlim=xlimALL,col='white',yaxs='i',axes=F)
      abline(v=yrs[which.max(obso)],lty=2)
      abline(h=yline,col="lightgrey")
      points(yrs,obso,col="black")
      if(i==ncomp)axis(1)
      if(i!=ncomp)axis(1,c(-10000,100000))
      if(j==1)axis(2)
      if(j!=1)axis(1,c(-10000,100000))
      #axis(3,c(-10000,100000))
      #axis(4,c(-10000,100000))
      lines(yrs,pred[[i]][,j],col=cols[j])
      if(i ==1)mtext(XLabs[j],col=cols[j],3,line=0.5,cex=labcex)
      if(j==1)mtext(YLabs[i],2,line=2.5,cex=labcex)
    }
    np<-pred[[i]]/rep(apply(pred[[i]],2,mean),each=length(yrs))
    matplot(yrs,np,col=cols,lty=1,type="l",ylim=c(0,max(np)*1.01),yaxs="i",axes=F)
    if(i!=ncomp)axis(1,c(-10000,100000))
    if(i==ncomp)axis(1)
    if(j==1)axis(2)
    if(j!=1)axis(1,c(-10000,100000))
   
    abline(v=yrs[apply(np,2,which.max)],col=cols,lty=2)
    
  }
  mtext("Index (N/tow)",2,outer=T,line=2.5,cex=labcex)
  mtext("Year",1,outer=T,line=1.8,cex=labcex)
  
}

Recent_Comp<-function(sums,runnams,ino=1){
  
  yrs<-sums[[1]]$years[length(sums[[1]]$years)-(14:0)]
  na<-dim(sums[[1]]$pred_catch_paa)[3]
  ncomp<-length(sums)
  
  getres<-function(x,ino){
    
    if(ino==1){
      res=log(x$catch_paa[1,,]/x$pred_catch_paa[,1,])#/x$pred_catch_paa[,1,]
    }else{
      res=log(x$ind_paa[ino-1,,]/x$pred_index_paa[,ino-1,])#/x$pred_catch_paa[,1,]
    }  
    ind<-length(x$years)-(14:0)
    res[ind,]
 
  }  
  
  res<-lapply(sums,getres,ino=ino)
  
  dat<-data.frame(Year=rep(yrs,na*ncomp),
                  Age=rep(rep(1:na,each=length(yrs)),ncomp),
                  Names=rep(runnams,each=length(yrs)*na),
                  log_rat=unlist(res),
                  Txt = round(exp(unlist(res)),2))
  
  
  gplot<-ggplot(dat) +
    geom_tile(aes(y=Year,x=Age,fill=log_rat))+
    scale_y_reverse()+
    scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))+
    facet_wrap(~Names)+
    geom_text(aes(y=Year,x=Age,label=Txt),size=2.5)+
    colorspace::scale_fill_continuous_diverging()
  
 gplot
  
}


TEG<-function(arr){ # make index for list calculation
  dim<-new('list')
  ndims<-length(arr)
  for(i in 1:ndims)dim[[i]]<-1:arr[i]
  as.matrix(expand.grid(dim))
}

Ucomp<-function(sums, runnams,inams=c("NMFS Spr","NFMS Fall","DFO Spr")){
  
  getobscat<-function(x)x$catch
  
  getobsindW<-function(x){
    nind<-dim(x$ind_paa)[1]
    WAA<-array(NA,dim(x$ind_paa))
    for(i in 1:nind)WAA[i,,]<-x$waa[i,,]
    
    #WAA<-x$BAA/x$NAA
    pind<-TEG(dim(x$ind_paa))
    totind<-pind[,2:1]
    Wind<-pind[,2:3]
    obs_IAA<-array(WAA[pind]*x$ind_paa[pind]*x$agg_indices[totind],dim(x$ind_paa))
    apply(obs_IAA,2:1,sum)
  }
  
  getpredcat<-function(x)x$pred_catch[,1]

  getpredindW<-function(x){
    nind<-dim(x$ind_paa)[1]
    WAA<-array(NA,dim(x$QAA))
    for(i in 1:nind)WAA[,i,]<-x$waa[i,,]
    qind<-TEG(dim(x$QAA))
    NFMind<-qind[,c(1,3)]
    Tind<-qind[,c(1:2)]
    pred_IAA<-array(WAA[qind]*x$QAA[qind]*x$NAA[NFMind]*exp((-x$MAA[NFMind]-x$FAA_tot[NFMind])*x$fracyr_indices[Tind]),dim(x$QAA))
    apply(pred_IAA,1:2,sum)
  }
 
  ncomp<-length(sums)
  nind<-dim(sums[[1]]$agg_indices)[2]
  ny<-dim(sums[[1]]$agg_indices)[1]
  
  obscat<-sapply(sums,getobscat)
  obsind<-array(sapply(sums,getobsindW),dim=c(ny,nind,ncomp))
  
  predcat<-sapply(sums,getpredcat)
  predind<-array(sapply(sums,getpredindW),dim=c(ny,nind,ncomp))
  
  ind<-TEG(dim(obsind))
  Cind<-ind[,c(1,3)]
  Uobs<-Upred<-array(NA,dim(obsind))
  Uobs[ind]<-obscat[Cind]/obsind[ind]
  Upred[ind]<-predcat[Cind]/predind[ind]
  Uobs[Uobs==Inf]<-NA
  Upred[Upred==Inf]<-NA
  
  UobsN <- Uobs / array(rep(apply(Uobs,2:3,mean,na.rm=T),each=ny),dim(Uobs)) 
  UpredN <- Upred / array(rep(apply(Upred,2:3,mean,na.rm=T),each=ny),dim(Upred)) 
  
  nc<-nind
  nr<-ncomp
  par(mfrow=c(nr,nc),mai=c(0.05,0.1,0.05,0.05),omi=c(0.4,0.8,0.25,0.01))
  cols<-c("black","red")
  x<-sums[[1]]
  for(cc in 1:ncomp){
    for(i in 1:nind){
    
     matplot(x$years,cbind(UobsN[,i,cc],UpredN[,i,cc]),type="l",ylab="",xlab="",col="white",lty=1,axes=F)
     abline(v=1:1000*10,col='light grey')
     abline(h=0:20,col="light grey")
     points(x$years,UobsN[,i,cc],col=cols[1])
     lines(x$years,UpredN[,i,cc],col=cols[2])
     if(cc==ncomp)axis(1)
     if(cc!=ncomp)axis(1,c(-10000,100000))
     if(i==1)axis(2)
     if(i!=1)axis(2,c(-10000,100000))
     if(cc==1 & i ==nind)legend('topright',legend=c("Observed","Predicted"),text.col=cols,bty='n')
     if(cc==1)mtext(inams[i],line=0.5)
     if(i==1)mtext(runnams[cc],2,line=2.8)
     
    }
  }
  
  mtext("Normalized implied harvest rate U = C(w) / I(w)",2,outer=T,line=4)
  
}


plotinds<-function(listy){
  par(mfrow=c(3,2),omi=c(0.5,0.5,0,0),mai=c(0.5,0.5,0.1,0.1))
  ind<-1:3
  yrs<-listy[[1]][,1]
  mus<-sapply(listy,function(x)x[,2])[,ind]
  nmus<-mus/rep(apply(mus,2,mean,na.rm=T),each=length(yrs))
  cvs<-sapply(listy,function(x)x[,3])[,ind]
  mus[mus==0]<-NA #  convert zeros to NAs
  
  cols<-c('black','green','red','blue')
  
  matplot(yrs,mus,ylab='Index',xlab='Year',type="l",col=cols,lty=1)
  abline(h=0,col='grey')
  matpoints(yrs,mus,col=cols,type="p",pch=1)
  legend('topleft',legend=c("NMFS_spring" , "NMFS_fall" , "DFO_spring" , "NMFS_fall"),bty='n',text.col=cols)
  
  matplot(yrs,nmus,ylab='Index',xlab='Year',type="l",col=cols,lty=1)
  abline(h=0,col='grey')
  matpoints(yrs,nmus,col=cols,type="p",pch=1)
  
  matplot(yrs,log(mus),ylab='Index',xlab='Year',type="l",col=cols,lty=1)
  abline(h=0,col='grey')
  matpoints(yrs,log(mus),col=cols,type="p",pch=1)
  
  for(i in 1:3){
    
    LB<-qnorm(0.025)*listy[[i]][,3]+log(listy[[i]][,2])
    UB<-qnorm(0.975)*listy[[i]][,3]+log(listy[[i]][,2])
    UB[UB==-Inf]<-NA
    LB[LB==-Inf]<-NA
    plot(yrs,log(listy[[i]][,2]),col=cols[i],ylim=c(min(LB,na.rm=T),max(UB,na.rm=T)),ylab='Log Index',xlab='Year',type='l',lwd=2)
    
    for(j in 1:length(LB)) lines(rep(listy[[i]][j,1],2),c(LB[j],UB[j]),col=cols[i])
    legend('topleft',legend=c("NMFS_spring" , "NMFS_fall" , "DFO_spring" , "NMFS_fall")[i],bty='n',text.col=cols[i])
    
  }
}

Indcomp2<-function(ilist, ind=1:3, 
                  nams=NULL, 
                  inams=c("NMFS Spr","NMFS Fall","DFO Spr"),
                  cols=c("black","red","blue","green")){
  

  ni<-length(ind)
  ncomp<-length(ilist)
  if(is.null(nams))nams=paste("Model",1:ncomp)
  
  par(mfrow=c(ni,3),omi=c(0.3,0.3,0.3,0),mai=c(0.3,0.3,0.1,0.1))
  
  getys<-function(x){
    ind1<-x[[1]]
    not0<-ind1[,2]!=0
    firsty<-ind1[match(T,not0),1]
    ys<-firsty:max(ind1[,1])
  }
  ys<-lapply(ilist,getys)
  yrs<-unlist(ys[which.max(sapply(ys,length))])
  yind<-lapply(ys,function(x,yrs)match(x,yrs),yrs=yrs)
  
  getis<-function(x,yrs){
    ind1<-x[[1]]
    not0<-ind1[,2]!=0
    firsty<-ind1[match(T,not0),1]
    ys<-firsty:max(ind1[,1])
    match(ys,ind1[,1])
  }
  iind<-lapply(ilist,getis,yrs=yrs)
 
 
  for(i in ind){
    mi<-list()
    for(j in 1:ncomp)mi[[j]]<-ilist[[j]][[i]][,2]
    indmat<-array(NA,c(length(yrs),ncomp))
    for(j in 1:ncomp)indmat[yind[[j]],j]<-mi[[j]][iind[[j]]]
    indmat[indmat==0]<-NA
    indmat_n<-indmat/rep(apply(indmat,2,mean,na.rm=T),each=nrow(indmat))
    indmat_l<-log(indmat)
    
    matplot(yrs,indmat,col="white",type="l",lty=1)
    abline(h=(0:100)*100,col="lightgrey")
    matplot(yrs,indmat,col=cols,type="l",lty=1,add=T)
    matplot(yrs,indmat,col=cols,type="p",pch=1,add=T)
    if(i==ind[1])legend('topleft',legend=nams,text.col=cols,bty='n')
    mtext(inams[i],2,line=2.2)
    
    matplot(yrs,indmat_n,col="white",type="l",lty=1)
    abline(h=(0:20),col="lightgrey")
    matplot(yrs,indmat_n,col=cols,type="l",lty=1,add=T)
    matplot(yrs,indmat_n,col=cols,type="p",pch=1,add=T)
    if(i==ind[1])mtext("Normalized",line=0.3)
    
    matplot(yrs,indmat_l,col="white",type="l",lty=1)
    abline(h=(-10:10),col="lightgrey")
    matplot(yrs,indmat_l,col=cols,type="l",lty=1,add=T)
    matplot(yrs,indmat_l,col=cols,type="p",pch=1,add=T)
    if(i==ind[1])mtext("Log(index)",line=0.3)
    
  }
   
 
}

