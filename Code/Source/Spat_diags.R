# Diagnostic tools


# ==== Correlation plots of abundance (D11a, D11b) ===============================================

getinds<-function(listo){
  
  inds<-sapply(listo,function(x)length(x$catch))
  minds<-min(inds)
  outy<-list()
  for(i in 1:length(inds))outy[[i]]<-(inds[i]-minds+1):inds[i]
  outy
}


betaplot <- function(listo,Inames=NULL,type="Obs"){
  
  runnam<-names(listo)
  inds<-getinds(listo)
  ni<-ncol(listo[[1]]$agg_indices)
  if(is.null(Inames))Inames<-paste("Index",1:ni)
  if(type=="Obs")Inames<-paste(Inames,"(observed)")
  if(type!="Obs")Inames<-paste(Inames,"(model predicted)")
  ncol<-ceiling(ni^0.5)
  nrow<-ceiling(ni/ncol)
  
  par(mfrow=c(nrow,ncol),omi=c(0.5,0.5,0,0),mai=c(0.5,0.5,0.5,0.01))
  
  for(i in 1:ni){
    
    if(type=="Obs"){
      x<-listo[[2]]$agg_indices[inds[[2]],i]
      y<-listo[[1]]$agg_indices[inds[[1]],i]
    }else{
      x<-listo[[2]]$pred_indices[inds[[2]],i]
      y<-listo[[1]]$pred_indices[inds[[1]],i]
      
    }
    keep<-x!=0 & y!=0
    x<-x[keep]
    y<-y[keep]
    #x<-(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))
    #y<-(y-min(y,na.rm=T))/(max(y,na.rm=T)-min(y,na.rm=T))
    x<-x/max(x,na.rm=T)
    y<-y/max(y,na.rm=T)
    
    plot(x,y,pch=19,col="white",xlab="",ylab="")
    mod<-function(par)sum((x-(y^par))^2)
    
    ysim<-xsim<-seq(0,1,length.out=100)
    betas<-c(1/4,1/2,1,2,4)
    yhyp<-sapply(betas,function(x)xsim^x)
    matplot(xsim,yhyp,col=c('red','orange','green','orange','red'),type="l",add=T,lty=c(1,1,1,1,1),lwd=c(1,1,2,1,1))
    
    
    opt<-optimize(mod,c(0.1,10))
    abline(v=pretty(seq(min(x),max(x),length.out=20)),col='grey')
    abline(h=pretty(seq(min(y),max(y),length.out=20)),col='grey')
    text(0.1,0.97,"Hyperdeplete",col="grey",cex=0.8)
    text(0.9,0.03,"Hyperstable",col="grey",cex=0.8)
    text(c(0.24,0.35,0.64,0.75),c(0.75,0.65,0.34,0.24),c("4","2","1/2","1/4"),col=c('red','orange','orange','red'),cex=0.8)
    
    xsim2<-ysim^opt$minimum
    lines(xsim2,ysim,col="blue")
    points(x,y,pch=19,col="black")
    legend('top',legend=paste("Est beta =",round(opt$minimum,4)),text.col="blue",bty='n')
    mtext(Inames[i],line=0.2,font=2,cex=0.9)
    
  }
  
  mtext(paste(runnam[1],"(larger area)"),2,line=0.5,outer=T,font=2)
  mtext(paste(runnam[2],"(smaller area)"),1,line=0.5,outer=T,font=2)
  
}

D1a <- function(listo,Inames=NULL,type="Obs")betaplot(listo, Inames,type)  # Observed
D1b <- function(listo,Inames=NULL,type="Pred")betaplot(listo, Inames,type) # Model predicted


D1c <- function(listo){
  
  inds<-getinds(listo)
  runnam<-names(listo)
  ny<-dim(listo[[2]]$BAA)[1]
  na<-dim(listo[[2]]$BAA)[2]
  VB<-array(NA,c(ny,length(runnam)))
  
  for(i in 1:length(listo)){
    
    FAAtot<-listo[[i]]$FAA_tot[inds[[i]],]
    FAAnorm<-FAAtot/apply(FAAtot,1,max)
    VB[,i]<-apply(FAAnorm*listo[[i]]$BAA[inds[[i]],],1,sum)
    
  }
  
  VBrat<-VB[,2]/VB[,1]
  B<-VB[,1]
  
  years<-listo[[2]]$years
  
  opt<-lm(y~x,data.frame(y=VBrat,x=B))
  coefs<-opt$coefficients
  xsim<-range(B)
  ysim<-coefs[1]+xsim*coefs[2]
  
  test<-anova(opt)
  Pval<-test$`Pr(>F)`[1]
  
  runnam<-names(listo)
  ylim=range(VBrat)
  if(ylim[1]>0.9)ylim[1]<-0.9
  if(ylim[2]<1.1)ylim[2]<-1.1
  
  htick<-pretty(seq(ylim[1],ylim[2],length.out=20))
  vtick<-pretty(seq(min(B),max(B),length.out=20))
  
  par(mfrow=c(1,2),mai=c(0.5,0.5,0.5,0.01),omi=c(0.5,0.5,0.01,0.01))
  plot(B,VBrat,col='white',xlab="",ylab="",ylim=ylim)
  
  abline(v=vtick,col='grey')
  abline(h=1,lwd=2,lty=1,col='green')
  abline(h=htick,col='grey')
  
  lines(xsim,ysim,col="blue")
  points(B,VBrat,pch=19)
  legend('topright',legend=c(paste("Intercept =", round(coefs[1],3)),
                             paste("Pval (slope) =", round(Pval,4))), text.col="blue",bty="n")
  
  plot(B,VBrat,col='white',xlab="",ylab="",ylim=ylim)
  text(B,VBrat,years,col=rainbow(length(years)),cex=0.9,font=1)
  mtext(paste("Vuln. Biom. for",runnam[1],"(larger area)"),1,line=1,font=2,outer=T)
  mtext(paste("Vuln. Biom. ratio ",runnam[2],"(smaller area) / ", runnam[1], " (larger area)"),2,line=1,outer=T,font=2,cex=0.7)
  
  
}


  
# ==== Ontogenetic selectivity (D12) =======================

  
D2 <- function(listo){
  
  
  years<-listo[[2]]$years
  yind<-floor(seq(1,length(years),length.out=4))
  ylab<-years[yind]
  ny<-length(years)
  na<-dim(listo[[2]]$BAA)[2]

  par(mfrow=c(4,2),omi=c(0.5,0.2,0.8,0),mai=c(0.3,0.4,0.1,0.01))
 
  for(i in 1:4){
    
    Fs<-sapply(listo,function(x)x$FAA[yind[i],1,])
    sel<-Fs/array(rep(apply(Fs,2,max),each=na),dim(Fs))
    matplot(1:na,sel,col="white",xlab="",ylab="")
    abline(h=pretty(seq(min(sel),max(sel),length.out=20)),col="grey")
    abline(v=0:(na+1),col="grey")
    matplot(1:na,sel,col=c("black","orange"),lty=1,lwd=2,type="l",add=T)
    if(i==1)legend('topleft',legend=names(listo),text.col=c("black","orange"),bty='n',text.font=2,cex=1.2)
    mtext(ylab[i],2,line=2.5,font=2)
    
    paa<-sapply(listo,function(x)x$pred_catch_paa[yind[i],1,])
    matplot(1:na,paa,col="white",xlab="",ylab="")
    abline(h=pretty(seq(min(paa),max(paa),length.out=20)),col="grey")
    abline(v=0:(na+1),col="grey")
    matplot(1:na,paa,col=c("black","orange"),lty=1,lwd=2,type="l",add=T)
    
  }
  
  mtext(c("Selectivity","Predicted catch fraction"),adj=c(0.25,0.78),outer=T,line=0.8,font=2)
  mtext("Age",1,line=1,outer=T,font=T)
  
}


# ==== F corr with biomass (D13) ===========================

D3a<-function(listo){
  
  inds<-getinds(listo)
  apF  <- lapply(listo,function(x)apply(x$FAA,1,max))
  B    <- apply(listo[[1]]$BAA[inds[[1]],],1,sum)
  Frat <- apF[[2]][inds[[2]]]/apF[[1]][inds[[1]]]
  years<-listo[[2]]$years
  runnam<-names(listo)
  
  opt<-lm(y~x,data.frame(y=Frat,x=B))
  coefs<-opt$coefficients
  xsim<-range(B)
  ysim<-coefs[1]+xsim*coefs[2]
  
  test<-anova(opt)
  Pval<-test$`Pr(>F)`[1]
  
  runnam<-names(listo)
  ylim=range(Frat)
  if(ylim[1]>0.9)ylim[1]<-0.9
  if(ylim[2]<1.1)ylim[2]<-1.1
  
  htick<-pretty(seq(ylim[1],ylim[2],length.out=20))
  vtick<-pretty(seq(min(B),max(B),length.out=20))
  
  par(mfrow=c(1,2),mai=c(0.5,0.5,0.5,0.01),omi=c(0.5,0.5,0.01,0.01))
  plot(B,Frat,col='white',xlab="",ylab="",ylim=ylim)
  
  abline(v=vtick,col='grey')
  abline(h=1,lwd=2,lty=1,col='green')
  abline(h=htick,col='grey')
  
  lines(xsim,ysim,col="blue")
  points(B,Frat,pch=19)
   legend('topright',legend=c(paste("Intercept =", round(coefs[1],3)),
                             paste("Pval (slope) =", round(Pval,3))), text.col="blue",bty="n")
 
  plot(B,Frat,col='white',xlab="",ylab="",ylim=ylim)
  text(B,Frat,years,col=rainbow(length(years)),cex=0.9,font=1)
  mtext(paste("Biomass for",runnam[1],"(larger area)"),1,line=1,font=2,outer=T)
  mtext(paste("ApF ratio ",runnam[2],"(smaller area) / ", runnam[1], " (larger area)"),2,line=1,outer=T,font=2,cex=0.7)
  
  
}


D3b<-function(listo){
  
  inds<-getinds(listo)
  apFt  <- lapply(listo,function(x)apply(x$FAA,1,max))
  apF<-NULL
  for(i in 1:length(listo))apF<-cbind(apF,apFt[[i]][inds[[i]]])
  runnam<-names(listo)
  years<-listo[[2]]$years
  par(mfrow=c(1,2),mai=c(0.4,1,0.5,0.01),omi=c(0.01,0.01,0.01,0.01))
  
  matplot(years,apF,col="white",xlab="",ylab="")
  abline(h=pretty(seq(min(apF),max(apF),length.out=20)),col="grey")
  abline(v=pretty(seq(min(years),max(years),length.out=20)),col="grey")
  
  matplot(years,apF,col=c("black","orange"),lty=1,lwd=2,type="l",add=T)
  legend('topright',legend=runnam,text.col=c("black","orange"),bty='n',text.font=2)
  mtext("Apical F",2,line=2.0,font=2)
  
  Frat <- apF[,2]/apF[,1]
  ylim=range(Frat)
  if(ylim[1]>0.9)ylim[1]<-0.9
  if(ylim[2]<1.1)ylim[2]<-1.1
  plot(years, Frat,xlab="",ylab="",col='white',ylim=ylim)
  abline(v=pretty(seq(min(years),max(years),length.out=20)),col="grey")
  
  htick<-pretty(seq(ylim[1],ylim[2],length.out=20))
  abline(h=htick,col='grey')
  abline(h=1,col="green",lwd=2)
  points(years, Frat,pch=19)
  lines(years,Frat)
  abline(h=mean(Frat),col="blue")
  legend('bottomright',legend=paste("Mean ratio =",round(mean(Frat),3)),text.col='blue',bty='n')
  mtext(paste0("Ratio apical F (",runnam[2]," / ",runnam[1],")"),2,line=2.0,font=2)
  
  
}  


D4<-function(listo, Inames=NULL){
  
  
  runnam<-names(listo)
  ni<-ncol(listo[[1]]$agg_indices)
  if(is.null(Inames))Inames<-paste("Index",1:ni)
  dat<-data.frame(run = rep(runnam,each=ni),index=rep(Inames,2),q=as.vector(sapply(listo,function(x)x$q[1,])))
  
  ggplot(dat,aes(fill=run, y=q, x=run)) + 
    geom_bar(position="dodge", stat="identity") +
    #ggtitle("Studying 4 species..") +
    facet_wrap(~index, scales = "free") +
    scale_fill_manual(values=c("orange","black"))+
    theme(legend.position="none") +
    xlab("")
  
}
  

D5<-function(listo){
  
  inds<-getinds(listo)
  runnam<-names(listo)
  nr<-length(runnam)
  ny<-dim(listo[[1]]$BAA)[1]
  na<-dim(listo[[1]]$BAA)[2]
  FAAlast<-sapply(listo,function(x)x$FAA_tot[nrow(x$FAA_tot),])
  FAAnorm<-FAAlast/array(rep(apply(FAAlast,2,max),each=na),dim(FAAlast))
  Blast    <- sapply(listo,function(x)x$BAA[nrow(x$BAA),])  
  VB<-apply(Blast*FAAnorm,2,sum)
  Ftarg<-sapply(listo,function(x)x$FXSPR[length(x$FXSPR)])
  OFL<-VB*Ftarg
  
  Cat<-rep(NA,length(listo))
  for(i in 1:length(listo))Cat[i]<-sum(listo[[i]]$catch[inds[[i]]])
  
  dat<-data.frame(quant=rep(c("Tot_Catch","VB","FXSPR","VB x FXSPR"),each=nr),run=rep(runnam,4),val=c(Cat,VB,Ftarg,OFL))
  colvals<-c("black","red","green","blue","grey","orange")[1:nr]
  ggplot(dat,aes(fill=run, y=val, x=run)) + 
    geom_bar(position="dodge", stat="identity") +
    #ggtitle("Studying 4 species..") +
    facet_wrap(~quant, scales = "free") +
    scale_fill_manual(values=colvals)+
    theme(legend.position="none") +
    xlab("")
  
}
  
  
D6<-function(listo){
  
  inds<-getinds(listo)
  runnam<-names(listo)
  ny<-dim(listo[[2]]$BAA)[1]
  na<-dim(listo[[2]]$BAA)[2]
  VB<-array(NA,c(ny,length(runnam)))
  
  for(i in 1:length(listo)){
    
    FAAtot<-listo[[i]]$FAA_tot[inds[[i]],]
    FAAnorm<-FAAtot/apply(FAAtot,1,max)
    VB[,i]<-apply(FAAnorm*listo[[i]]$BAA[inds[[i]],],1,sum)
   
  }
  
  years<-listo[[2]]$years
  par(mfrow=c(1,2),mai=c(0.4,1,0.5,0.01),omi=c(0.01,0.01,0.01,0.01))
  
  matplot(years,VB,col="white",xlab="",ylab="")
  abline(h=pretty(seq(min(VB),max(VB),length.out=20)),col="grey")
  abline(v=pretty(seq(min(years),max(years),length.out=20)),col="grey")
  
  matplot(years,VB,col=c("black","orange"),lty=1,lwd=2,type="l",add=T)
  legend('topright',legend=runnam,text.col=c("black","orange"),bty='n',text.font=2)
  mtext("Vulnerable Biomass",2,line=2.0,font=2)
  
  Vrat <- VB[,2]/VB[,1]
  ylim=range(Vrat)
  if(ylim[1]>0.9)ylim[1]<-0.9
  if(ylim[2]<1.1)ylim[2]<-1.1
  plot(years, Vrat,xlab="",ylab="",col='white',ylim=ylim)
  abline(v=pretty(seq(min(years),max(years),length.out=20)),col="grey")
  
  htick<-pretty(seq(ylim[1],ylim[2],length.out=20))
  abline(h=htick,col='grey')
  abline(h=1,col="green",lwd=2)
  points(years, Vrat,pch=19)
  lines(years,Vrat)
  abline(h=mean(Vrat),col="blue")
  legend('topright',legend=paste("Mean ratio =",round(mean(Vrat),3)),text.col='blue',bty='n')
  mtext(paste0("Ratio vulnerable biomass (",runnam[2]," / ",runnam[1],")"),2,line=2.0,font=2)
  
}



#matp(sapply(listo,function(x)log(x$N[,1])),col=cols, ylab='log Recruitment (N age 0)')

D7a<-function(listo){
  
  inds<-getinds(listo)
  runnam<-names(listo)
  ny<-dim(listo[[2]]$BAA)[1]
  na<-dim(listo[[2]]$BAA)[2]
  Rec<-array(NA,c(ny,length(listo)))
  for(i in 1:length(listo))    Rec[,i]<-listo[[i]]$NAA[inds[[i]],1]
  
  years<-listo[[2]]$years
  par(mfrow=c(2,2),mai=c(0.4,1,0.5,0.01),omi=c(0.01,0.01,0.01,0.01))
  
  # recruitment TS
  matplot(years,Rec,col="white",xlab="",ylab="")
  abline(h=pretty(seq(min(Rec),max(Rec),length.out=20)),col="grey")
  abline(v=pretty(seq(min(years),max(years),length.out=20)),col="grey")
  
  matplot(years,Rec,col=c("black","orange"),lty=1,lwd=2,type="l",add=T)
  legend('topleft',legend=runnam,text.col=c("black","orange"),bty='n',text.font=2)
  mtext("Est. Recruitment",2,line=2.0,font=2)
  
  # log recruitment TS
  matplot(years,log(Rec),col="white",xlab="",ylab="")
  abline(h=pretty(seq(min(log(Rec)),max(log(Rec)),length.out=20)),col="grey")
  abline(v=pretty(seq(min(years),max(years),length.out=20)),col="grey")
  
  matplot(years,log(Rec),col=c("black","orange"),lty=1,lwd=2,type="l",add=T)
  legend('topleft',legend=runnam,text.col=c("black","orange"),bty='n',text.font=2)
  
  mtext("Log Est. Recruitment",2,line=2.0,font=2)
 
  Vrat <-  Rec[,1]- Rec[,2]
  ylim=range(Vrat)
  if(ylim[1]>0.9)ylim[1]<-0.9
  if(ylim[2]<1.1)ylim[2]<-1.1
  plot(years, Vrat,xlab="",ylab="",col='white',ylim=ylim)
  abline(v=pretty(seq(min(years),max(years),length.out=20)),col="grey")
  
  htick<-pretty(seq(ylim[1],ylim[2],length.out=20))
  abline(h=htick,col='grey')
  abline(h=1,col="green",lwd=2)
  points(years, Vrat,pch=19)
  lines(years,Vrat)
  abline(h=mean(Vrat),col="blue")
  legend('top',legend=paste("Mean dif. =",round(mean(Vrat),1)),text.col='blue',bty='n')
  mtext(paste0("Dif recruitment (",runnam[1]," - ",runnam[2],")"),2,line=2.0,font=2)
 
  Vrat <- Rec[,1]/Rec[,2]
  ylim=range(Vrat)
  if(ylim[1]>0.9)ylim[1]<-0.9
  if(ylim[2]<1.1)ylim[2]<-1.1
  plot(years, Vrat,xlab="",ylab="",col='white',ylim=ylim)
  abline(v=pretty(seq(min(years),max(years),length.out=20)),col="grey")
  
  htick<-pretty(seq(ylim[1],ylim[2],length.out=20))
  abline(h=htick,col='grey')
  abline(h=1,col="green",lwd=2)
  points(years, Vrat,pch=19)
  lines(years,Vrat)
  abline(h=mean(Vrat),col="blue")
  legend('top',legend=paste("Mean ratio =",round(mean(Vrat),2)),text.col='blue',bty='n')
  mtext(paste0("Ratio recruitment (",runnam[1]," / ",runnam[2],")"),2,line=2.0,font=2)

}

D7b<-function(listo){
  
  inds<-getinds(listo)
  runnam<-names(listo)
  ny<-dim(listo[[2]]$BAA)[1]
  na<-dim(listo[[2]]$BAA)[2]
  Rec<-array(NA,c(ny,length(runnam)))
  for(i in 1:length(listo))    Rec[,i]<-listo[[i]]$NAA[inds[[i]],1]
  
  years<-listo[[2]]$years
  par(mfrow=c(1,2),mai=c(0.6,1,0.5,0.01),omi=c(0.01,0.01,0.01,0.01))
  
  # recruitment TS
  
   
  plot(Rec[,1],Rec[,2],col="white",xlab="",ylab="")
  lns<-pretty(seq(min(Rec),max(Rec),length.out=30))
  abline(h=lns,col="grey")
  abline(v=lns,col="grey")
  
  opt<-lm(y~x,data.frame(y=Rec[,2],x=Rec[,1]))
  coefs<-opt$coefficients
  xs<-seq(0,max(Rec[,1]),length.out=20)
  ys<-predict(opt,newdata=data.frame(x=xs))
  lines(xs,ys,col="blue")
  
  points(Rec[,1],Rec[,2],col=c("black"),pch=19)
  mtext(paste("Est. Rec. ",runnam[1],"(larger area)"),1,line=2,font=2)
  mtext(paste("Est. Rec. ",runnam[2],"(smaller area)"),2,line=2,font=2)
  
  #legend('topleft',legend=c(paste("Slope =", round(coefs[2],3)),
  #                           paste("Pval (slope) =", round(Pval,3))), text.col="blue",bty="n")
  
  Rec<-log(Rec)
  plot(Rec[,1],Rec[,2],col="white",xlab="",ylab="")
  lns<-pretty(seq(min(Rec),max(Rec),length.out=20))
  abline(h=lns,col="grey")
  abline(v=lns,col="grey")
  
  opt<-lm(y~x,data.frame(y=Rec[,2],x=Rec[,1]))
  coefs<-opt$coefficients
  xs<-seq(0,max(Rec[,1]),length.out=20)
  ys<-predict(opt,newdata=data.frame(x=xs))
  lines(xs,ys,col="blue")
  
  points(Rec[,1],Rec[,2],col=c("black"),pch=19)
  mtext(paste("Est. log Rec. ",runnam[1],"(larger area)"),1,line=2.,font=2)
  mtext(paste("Est. log Rec. ",runnam[2],"(smaller area)"),2,line=2,font=2)
 
  #legend('topleft',legend=c(paste("Slope =", round(coefs[2],3)),
   #                         paste("Pval (slope) =", round(Pval,3))), text.col="blue",bty="n")
  
 
}


D8<-function(listo){
  AAResComp(listo, YLabs=names(listo),XLabs=c("Catch","NMFS Spr","NMFS Fall","DFO Spr"),Nscaled=F)
}



