# some plotting stuff
makeTransparent<-function (someColor, alpha = 100){
  newColor <- col2rgb(someColor)
  apply(newColor, 2, function(curcoldata) {
    rgb(red = curcoldata[1], green = curcoldata[2], blue = curcoldata[3],
        alpha = alpha, maxColorValue = 255)
  })
}


plotallsel<-function(dirs,runnams=NULL){
  nn<-length(dirs)
  if(is.null(runnams))runnams=paste("Model",1:nn)
  ncol<-floor(nn^0.5)
  nrow<-ceiling(nn/ncol)
  par(mfrow=c(nrow,ncol),mai=c(0.4,0.1,0.3,0.05),omi=c(0.4,0.6,0.25,0.01))
  for(i in 1:nn){
    outfile=dirs[i]#paste0(dirs[i],"/Fit.rda")
    Fit<-readRDS(outfile)
    plot_wham_sel(Fit)
    mtext(runnams[i],line=0.3,adj=0.25)
  }
}


plot_wham_sel<-function(m1,indnam =c("Catch block 1","Catch block 2","NMFS_s","NMFS_f","DFO_s")){

  cols<-c('black','red','green','blue','orange','grey','pink','brown')

  sels<-sapply(m1$rep$selAA,function(x)x[1,])

  na<-dim(sels)[1]
  yadj<-array(rep(seq(0,-0.005,length.out=dim(sels)[2]),each=na),dim(sels))
  ##fixage<-matrix(0:(na-1),nrow=nf,ncol=na, byrow=T)[is.na(sels)]
  #fixage<-c(7,5,7,2)-1

  matplot(1:na,sels+yadj,col=cols,type="l",lty=1,lwd=2,xlab="",ylab="")
  abline(h=seq(0,1,by=0.2),col='light grey')
  abline(v=seq(0,100,by=1),col='light grey')
  #abline(v=fixage,col=cols,lty=2)
  matplot(1:na,sels+yadj,col=cols,type="l",lty=1,lwd=2,add=T)
  abline(h=seq(0,1,by=0.2),col='light grey')
  legend('bottom',legend=indnam,text.col=cols,bty='n')
  mtext("Age",1,line=1.8)
  mtext("Selectivity",2,line=2.2)

}

matp<-function(mat,xs = 1969:2019,xlab="",ylab="",cols,yint=T,ylim=NULL,lwd=2,lty=1){
  rng<-range(mat,na.rm=T)
  if(is.null(ylim))ylim=rng
  if(yint)ylim[1]<-0
  matplot(xs,mat,type='l',col='white', yaxs='i',xlab="",ylab="",ylim=ylim)

  dif<-(ylim[2]-ylim[1])/5
  lys<-pretty(seq(ylim[1]-dif,ylim[2]+dif,length.out=20))
  abline(h=lys,col='light grey')
  abline(v=seq(1900,2100,by=10),col='light grey')

  matplot(xs,mat,type='l',col=cols, add=T, lwd=lwd,lty=lty)
  mtext(xlab,1,line=2.2)
  mtext(ylab,2,line=2.5)
}

Comp1<-function(listo,namey,cols=c("black","red","green","blue","orange","grey","purple","pink","darkred","darkblue"),lwd=2,include_nll=T){

  lcols<-makeTransparent(cols,90)
  isretro<-!is.null(listo[[1]]$rho)
  if(isretro){
    rho_SSB<- sapply(listo,function(x)x$rho[1])
    rho_Fbar<-sapply(listo,function(x)x$rho[2])
  }
  AIC<-sapply(listo,function(x)x$AIC)
  dAIC<-AIC-min(AIC)
  totLH<-sapply(listo,function(x)x$nll[1])

  layout(matrix(c(1,2,5,3,4,5),byrow=T,nrow=2),widths=c(1,1,0.5))
  par(mai=c(0.4,0.55,0.05,0.05), omi=c(0.3,0.05,0.05,0.05))
  matp(sapply(listo,function(x)x$SSB),xs=listo[[1]]$years, col=lcols, ylab="SSB",yint=T,lwd=lwd)
  if(isretro)legend('topleft',legend=paste("rho =",round(rho_SSB,3)),text.col=cols,bty='n')
  if(include_nll)legend('top',legend=paste("nll =",round(totLH,1)),text.col=cols,bty='n')
  matp(sapply(listo,function(x)log(x$N[,1])),xs=listo[[1]]$years,col=lcols, ylab='log Recruitment (N age 1)',lwd=lwd)
  if(include_nll)legend('bottomright',legend=paste("dAIC =",round(dAIC,1)),text.col=cols,bty='n')
  matp(sapply(listo,function(x)apply(x$FAA,1,max)),xs=listo[[1]]$years,col=lcols,ylab="Fully selected F",yint=T,lwd=lwd)
  if(isretro)legend('topright',legend=c("NConv","Conv")[as.integer(sapply(listo,function(x){x$conv==0}))+1],text.col=cols,bty='n')
  
  matp(sapply(listo,function(x)x$Fbar),xs=listo[[1]]$years,col=lcols, ylab="Fbar",yint=T,lwd=lwd)
  if(isretro)legend('topright',legend=paste("rho =",round(rho_Fbar,3)),text.col=cols,bty='n')
  par(mai=rep(0.1,4))
  plot(1,1,axes=F,xlab="",ylab="",col="white"); legend('left',legend=namey, text.col=cols,bty='n',cex=1.2)

}

Comp2<-function(listo,namey=NULL,
                cols=c("black","red","green","blue","orange","grey","purple","pink","darkred","darkblue"),
                lty=1, lwd=2, include_nll=F, include_conv=F,Frhopos='right',alpha=90){
  
  
  lcols<-makeTransparent(cols,alpha)
  isretro<-!is.null(listo[[1]]$rho)
  if(isretro){
    rho_SSB<- sapply(listo,function(x)x$rho[1])
    rho_Fbar<-sapply(listo,function(x)x$rho[2])
    for(i in 1:length(rho_SSB))if(is.null(rho_SSB[[i]]))rho_SSB[[i]]<-NA
    for(i in 1:length(rho_Fbar))if(is.null(rho_Fbar[[i]]))rho_Fbar[[i]]<-NA
    rho_SSB<-unlist(rho_SSB)
    rho_Fbar<-unlist(rho_Fbar)
  }
  
  
  AIC<-sapply(listo,function(x)x$AIC)
  dAIC<-AIC-min(AIC,na.rm=T)
  totLH<-sapply(listo,function(x)x$nll[1])
  exc<-!(!is.nan(totLH)|!is.na(totLH))
  
  layout(matrix(c(1,2,5,3,4,5),byrow=T,nrow=2),widths=c(1,1,0.5))
  par(mai=c(0.4,0.55,0.05,0.05), omi=c(0.3,0.05,0.05,0.05))
  
  ncomp<-length(listo)
  if(is.null(namey))namey=paste("Model",1:ncomp)
  
  yrs<-lapply(listo,function(x)x$years)
  lens<-sapply(yrs,function(x)length(x))
  years<-yrs[[which.max(lens)]]
  nt<-length(years)
  SSB <- N1 <- FAA <- Fbar <- array(NA,c(nt,ncomp))
  for(i in 1:ncomp){
    ny<-length(listo[[i]]$years)
    ind<-nt-((ny-1):0)
    SSB[ind,i]<-listo[[i]]$SSB
    N1[ind,i]<-listo[[i]]$NAA[,1]
    FAA[ind,i]<-apply(listo[[i]]$FAA,1,max)
    Fbar[ind,i]<-listo[[i]]$Fbar
  }  
  SSB[,exc]<-NA
  N1[,exc]<-NA
  FAA[,exc]<-NA
  Fbar[,exc]<-NA
    
  matp(SSB/1000,xs=years, col=lcols, ylab="SSB (kt)",yint=T,lwd=lwd,lty=lty)
  if(isretro)legend('topleft',legend=paste("rho =",round(rho_SSB,3)),text.col=cols,bty='n')
  if(include_nll)legend('top',legend=paste("nll =",round(totLH,1)),text.col=cols,bty='n')
  
  matp(N1,xs=years,col=lcols, ylab='Recruitment (N age 1)',lwd=lwd,,lty=lty)
  if(include_nll)legend('topleft',legend=paste("dAIC =",round(dAIC,1)),text.col=cols,bty='n')
  
  matp(FAA,xs=years,col=lcols,ylab="Fully selected F",yint=T,lwd=lwd,,lty=lty)
  if(include_conv)legend('topright',legend=c("NConv","Conv")[as.integer(sapply(listo,function(x){x$conv==0}))+1],text.col=cols,bty='n')
  
  matp(Fbar,xs=years,col=lcols, ylab="Mean F",yint=T,lwd=lwd,,lty=lty)
  if(isretro)legend(Frhopos,legend=paste("rho =",round(rho_Fbar,3)),text.col=cols,bty='n')
  par(mai=rep(0.1,4))
  plot(1,1,axes=F,xlab="",ylab="",col="white"); legend('left',legend=namey, text.col=cols,bty='n',cex=1.2)
  
}


Comp3<-function(listo,namey=NULL,
                cols=c("black","red","green","blue","orange","grey","purple","pink","darkred","darkblue"),
                lty=1, lwd=2, include_nll=F, include_conv=F){
  
  
  lcols<-makeTransparent(cols,90)
  isretro<-!is.null(listo[[1]]$rho)
  if(isretro){
    rho_SSB<- sapply(listo,function(x)x$rho[1])
    rho_Fbar<-sapply(listo,function(x)x$rho[2])
  }
  for(i in 1:length(rho_SSB))if(is.null(rho_SSB[[i]]))rho_SSB[[i]]<-NA
  for(i in 1:length(rho_Fbar))if(is.null(rho_Fbar[[i]]))rho_Fbar[[i]]<-NA
  rho_SSB<-unlist(rho_SSB)
  rho_Fbar<-unlist(rho_Fbar)
  
  AIC<-sapply(listo,function(x)x$AIC)
  dAIC<-AIC-min(AIC,na.rm=T)
  totLH<-sapply(listo,function(x)x$nll[1])
  exc<-!(!is.nan(totLH)|!is.na(totLH))
  
  layout(matrix(c(1,2,3),byrow=T,nrow=1),widths=c(1,1,0.5))
  par(mai=c(0.4,0.55,0.05,0.05), omi=c(0.3,0.05,0.05,0.05))
  
  ncomp<-length(listo)
  if(is.null(namey))namey=paste("Model",1:ncomp)
  
  yrs<-lapply(listo,function(x)x$years)
  lens<-sapply(yrs,function(x)length(x))
  years<-yrs[[which.max(lens)]]
  nt<-length(years)
  SSB <- N1 <- FAA <- Fbar <- array(NA,c(nt,ncomp))
  for(i in 1:ncomp){
    ny<-length(listo[[i]]$years)
    ind<-nt-((ny-1):0)
    SSB[ind,i]<-listo[[i]]$SSB
    N1[ind,i]<-listo[[i]]$NAA[,1]
    FAA[ind,i]<-apply(listo[[i]]$FAA,1,max)
    Fbar[ind,i]<-listo[[i]]$Fbar
  }  
  SSB[,exc]<-NA
  N1[,exc]<-NA
  FAA[,exc]<-NA
  Fbar[,exc]<-NA
  
  matp(SSB/1000,xs=years, col=lcols, ylab="Spawning Stock Biomass (kt)",yint=T,lwd=lwd,lty=lty)
  if(isretro)legend('topleft',legend=paste("rho =",round(rho_SSB,3)),text.col=cols,bty='n')
  if(include_nll)legend('top',legend=paste("nll =",round(totLH,1)),text.col=cols,bty='n')
  
  matp(Fbar,xs=years,col=lcols, ylab="Mean Exploitation Rate (Fbar)",yint=T,lwd=lwd,,lty=lty)
  if(isretro)legend('topright',legend=paste("rho =",round(rho_Fbar,3)),text.col=cols,bty='n')
  par(mai=rep(0.1,4))
  plot(1,1,axes=F,xlab="",ylab="",col="white"); legend('left',legend=namey, text.col=cols,bty='n',cex=1.2)
  
}


Prof<-function(listo,
               cols=c("black","grey","green","green","green","blue","blue","blue","red","orange","pink","purple"),
               ltys=c(1,1,1,2,3,1,2,3,1,1,1,1),parnam=NULL,parlab="",keep=NULL){

  if(is.null(parnam))parnam=1:length(listo)
  
  par(mai=c(0.8,0.8,0.05,0.05))
  LHmat<-sapply(listo,function(x)x$nll[c(1:11,14)]) #[c(1,2,3,5)])
  LHmat<-LHmat-apply(LHmat,1,min)
  if(is.null(keep))keep<-1:nrow(LHmat)
  LHmat<-LHmat[keep,]
  matp(mat=t(LHmat),xs=parnam,cols=cols[keep], lty=ltys[keep],ylab="delta nll",yint=T)
  legend('topleft',legend=row.names(LHmat), lty=ltys[keep],col=cols[keep],text.col=cols[keep],bty='n',cex=1.6)
  mtext(parlab,1,outer=F,line=1.5,cex=0.75)

}

Prof2<-function(listo,
               cols=c("black","green","green","green","blue","blue","blue","red","orange"),
               leg=c("NMFS Spr Index","NMFS Fall Index","DFO Spr Index",
                     "NMFS Spr Comp","NMFS Fall Comp","DFO Spr Comp",
                     "Fleet Catch","Fleet Comp"),
               def=0.2,
               ltys=c(1,1,2,3,1,2,3,1,1),parnam=NULL,parlab=""){
  
  if(is.null(parnam))parnam=1:length(listo)
  
  par(mai=c(0.8,0.8,0.05,0.05))
  LHmat<-sapply(listo,function(x)x$nll[c(1:11,14)]) #[c(1,2,3,5)])
  LHmat<-LHmat[3:10,]
  LHmat<-rbind(apply(LHmat,2,sum),LHmat)
  val<-LHmat[1,]
  sdv<-sd(val)
  muv<-mean(val)
  keep <- !(val>muv+sdv)
  LHmat<-LHmat[,keep]
  LHmat<-LHmat-apply(LHmat,1,min)
  row.names(LHmat)<-c("Total",leg)#"Index 1",    "Index 2",    "Index 3",    "Index comp 1", "Index comp 2", "Index comp 3", "Catch", "Catch comp")
  
  matp(mat=t(LHmat),xs=parnam[keep],cols=cols, lty=ltys,ylab="Delta nll",yint=T)
  abline(v=def,col='grey')
  legend('topleft',legend=row.names(LHmat), lty=ltys,col=cols,text.col=cols,bty='n',cex=1)
  mtext(parlab,1,outer=F,line=2.1,cex=0.95)
  
}




CompFAA<-function(listo,namey, cols=c("black","red","green","blue","orange","grey"),curYr=2019,yind=c(1970,1980, 1990, 2000, 2010, 2015, 2019)){


  par(mfrow=c(ceiling(length(yind)/2),2),mai=c(0.35,0.25,0.02,0.05), omi=c(0.25,0.3,0.05,0.05))
  ny<-dim(listo[[1]]$FAA)[1]
  na<-dim(listo[[1]]$FAA)[3]
  syr<-curYr-ny+1

  ymax<-max(sapply(listo,function(x)max(x$FAA[ny+(yind-curYr),1,])))

  for(y in 1:length(yind)){

    yi<-ny+(yind[y]-curYr)
    matp(sapply(listo,function(x)x$FAA[yi,1,]),xs=0:(na-1),col=cols,ylab="",yint=T,ylim=c(0,ymax))
    legend('topleft',legend=yind[y],bty='n')

  }
  plot(1,1,axes=F,xlab="",ylab="",col="white"); legend('left',legend=namey, text.col=cols,bty='n',cex=1.6)
  mtext("Fishing mortality rate",2,outer=T,line=0.6,cex=1)
  mtext('Age',1,outer=T,line=0.6,cex=1)

}

CompVPA<-function(listo,namey,cols=c("black","red","green","blue","orange","grey")){

  layout(matrix(c(1,2,5,3,4,5),byrow=T,nrow=2),widths=c(1,1,0.5))
  par(mai=c(0.4,0.55,0.05,0.05), omi=c(0.3,0.05,0.05,0.05))
  matp(sapply(listo,function(x)x$SSB),cols=cols, ylab="SSB")
  matp(sapply(listo,function(x)log(x$N[,2])),col=cols, ylab='log N age 1')
  matp(sapply(listo,function(x)apply(x$FAA,1,max)),col=cols,ylab="Fully selected F")
  matp(sapply(listo,function(x)x$Fbar),col=cols, ylab="Fbar")
  plot(1,1,axes=F,xlab="",ylab="",col="white"); legend('left',legend=namey, text.col=cols,bty='n',cex=1.6)
  mtext('Year',1,outer=T,line=0.8,cex=1)

}

docomp<-function(files,mnams=NULL){

  if(is.null(mnams))mnams=paste("m",1:length(files))
  listy<-list()
  for(i in 1:length(files))listy[[i]]=readRDS(files[i])
  names(listy)=mnams
  compare_wham_models(listy)

}

NAA_RE_plot<-function(Fit){

  RE<-Fit$rep$NAA_devs
  mods<-list(Fit)

  years = Fit$input$data$year1_model+(2:Fit$input$data$n_years_model)-1
  n_years = length(years)
  n_ages = Fit$env$data$n_ages
  ages <- 1:n_ages
  df.plot<-data.frame(RE=as.vector(RE),Year=rep(years,n_ages),Age=rep(ages,each=n_years))

  ggplot(df.plot, ggplot2::aes(x=Year, y=Age)) +
          geom_tile(aes(fill=RE)) +
    scale_fill_gradient2(name = "", low = scales::muted("blue"), mid = "white", high = scales::muted("red"))

 # +
   #       geom_label(aes(x=Year, y=Age, label=lab), size=5, alpha=1, #fontface = "bold",
  #                   data=data.frame(Year=1976.5, Age=5.8, lab=df.mods$Model[plot.mods], NAA_lab=factor(NAA_lab, levels=NAA_lab_levels), GSI_how=factor(GSI_how, levels=GSI_lab_levels))) +
    #      scale_x_continuous(expand=c(0,0)) +
     #     scale_y_continuous(expand=c(0,0)) +
      #    theme_bw() +
       #   facet_grid(rows=vars(NAA_lab), cols=vars(GSI_how), drop=F) +


}


dolines<-function(xs,ys,col='black',lwd=2){
  for(i in 1:nrow(xs))lines(xs[i,],ys[i,],col=col,lwd=lwd)  
}

IndComp<-function(listo,YLabs=NULL,XLabs=NULL,ylimobs=F){
  
  labcex=0.85
  ni<-ncol(listo[[1]]$agg_indices)
  obs<-lapply(listo,function(x)x$agg_indices)
  pred<-lapply(listo,function(x)x$pred_indices)
  xlimALL<-range(lapply(listo,function(x)x$years))
  
  obsmax<-t(sapply(obs,function(x)apply(x,2,max,na.rm=T)))
  predmax<-t(sapply(pred,function(x)apply(x,2,max,na.rm=T)))
  maxi<-apply(rbind(obsmax,predmax),2,max)
  if(ylimobs)maxi=apply(obsmax,2,function(x)max(x)*1.01)
  ncomp<-length(listo)
  if(is.null(YLabs))YLabs = paste("Fit",1:ncomp)
  if(is.null(XLabs))XLabs = c(paste("Index",1:ni))
   
  
  par(mfrow=c(ncomp,ni),mai=c(0.02,0.25,0.00,0.025),omi=c(0.5,0.6,0.25,0.01))
  for(i in 1:ncomp){
    
    yrs<-listo[[i]]$years
    
    for(j in 1:ni){
      yline<-pretty(seq(0,maxi[j]*2,length.out=18))
      obso<-obs[[i]][,j]
      obso[obso==0]<-NA
      plot(yrs,obso,ylim=c(0,maxi[j]*1.02),xlim=xlimALL,col='white',yaxs='i',axes=F)
      abline(h=yline,col="lightgrey")
      points(yrs,obso,pch=19,col="#ff000080",cex=0.9)
      if(i==ncomp)axis(1)
      if(i!=ncomp)axis(1,c(-10000,100000))
      if(j==1)axis(2)
      if(j!=1)axis(2,c(-10000,100000))
      #axis(3,c(-10000,100000))
      #axis(4,c(-10000,100000))
      lines(yrs,pred[[i]][,j],col="blue")
       if(i ==1)mtext(XLabs[j],3,line=0.5,cex=labcex)
      if(j==1)mtext(YLabs[i],2,line=2.5,cex=labcex)
    }
    
  }
  mtext("Index (N/tow)",2,outer=T,line=2.5,cex=labcex)
  mtext("Year",1,outer=T,line=2.3,cex=labcex)
  
}


IndResComp<-function(listo,YLabs=NULL,XLabs=NULL,resrange=c(-2,2)){
  
  labcex=0.85
  ni<-ncol(listo[[1]]$agg_indices)
  res<-lapply(listo,function(x)log(x$agg_indices/x$pred_indices))
  res<-lapply(res,function(x){x[x==-Inf]<-NA; x})
  xlimALL<-range(lapply(sums,function(x)x$years))
  #resrange<-quantile(unlist(res),c(0.01,0.99),na.rm=T)
  
   
  ncomp<-length(listo)
  if(is.null(YLabs))YLabs = paste("Fit",1:ncomp)
  if(is.null(XLabs))XLabs = c(paste("Index",1:ni))
  
  par(mfrow=c(ncomp,ni),mai=c(0.3,0.25,0.1,0.05),omi=c(0.4,0.6,0.25,0.01))
  for(i in 1:ncomp){
    
    yrs<-listo[[i]]$years
    
    for(j in 1:ni){
      
      plot(yrs,res[[i]][,j],ylim=resrange,xlim=xlimALL,col='white')
      abline(h=seq(-10,10,by=0.5),col='lightgrey')
      dolines(cbind(yrs,yrs),cbind(0,res[[i]][,j]))
      if(i ==1)mtext(XLabs[j],3,line=0.5,cex=labcex)
      if(j==1)mtext(YLabs[i],2,line=3.5,cex=labcex)
    }
    
  }
  mtext("log residual (obs/pred)",2,outer=T,line=2.5,cex=labcex)
  mtext("Year",1,outer=T,line=0.8,cex=labcex)
 
}

AAResComp<-function(sums,Nscaled=T,YLabs=NULL,XLabs=NULL,Nscaler=50){
  labcex=0.85
  ncomp<-length(sums)
  out<-sums[[1]]
  ni<-dim(out$pred_index_paa)[2]
  if(is.null(YLabs))YLabs = paste("Fit",1:ncomp)
  if(is.null(XLabs))XLabs = c("Catch",paste("Index",1:ni))
  par(mfrow=c(ncomp,1+ni),mai=c(0.05,0.1,0.05,0.05),omi=c(0.4,0.6,0.25,0.01))

  ylimALL<-rev(range(lapply(sums,function(x)x$years)))
 
  for(i in 1:ncomp){
    
    out<-sums[[i]]
    ytick<-out$years#[length(out$years):1]
    pred=out$pred_catch_paa[,1,]
    obs=out$catch_paa[1,,]
    N=out$catch_Neff[,1]
    AAResPlot(obs,pred,N,Nscaled,ytick=ytick,ylab=T,xlab=(i==ncomp),ylim=ylimALL,Nscaler=Nscaler)
    
    if(i ==1)mtext(XLabs[1],3,line=0.5,cex=labcex)
    mtext(YLabs[i],2,line=3.5,cex=labcex)
    for(j in 1:ni){

      pred=out$pred_index_paa[,j,]
      obs=out$ind_paa[j,,]
      N=out$index_Neff[,j]
      AAResPlot(obs,pred,N,Nscaled,ytick=ytick,xlab=(i==ncomp),ylim=ylimALL,Nscaler=Nscaler,leg=(j==ni&i==1))
      if(i ==1)mtext(XLabs[j+1],3,line=0.5,cex=labcex)
    }
  }
  mtext("Age",1,outer=T,line=1.5,cex=labcex)
}

AAResPlot<-function(obs,pred,N,Nscaled=F,ytick,xlab=F,ylab=F,ylim=NULL,Nscaler,leg=F){

  
  N2<-N3<-N*obs
  if(!Nscaled) N3=1
  scaler<-2.5
  res = N3*(obs - pred)/sqrt(N2 * pred * (1 - pred)) # N2 = 50; obs=0.1; pred=0.05; N2 * (obs - pred)/sqrt(N2 * pred * (1 - pred))
  #res = (obs - pred)/sqrt(pred * (1 - pred))
  bubwidth = abs(res * scaler)^0.5
  pos<-res>0
  na <- ncol(obs)
  ny <- nrow(obs)
  xs<-rep(1:na,each=ny)
  ys<-ytick[rep(1:ny,na)]
  if(is.null(ylim))ylim=rev(range(ytick))
  plot(xs[pos],ys[pos],ylim=ylim,cex=bubwidth[pos],xlab="",ylab="",axes=F)
  if(xlab)axis(1,0:9,0:9)
  if(!xlab)axis(1,c(-1000,1000),c(-1000,1000))
  if(ylab)axis(2,(191:220)*10,(191:220)*10,las=2)
  if(!ylab)axis(2,c(-1000,10000),c(-1000,10000))
  axis(3,c(-1000,1000),c(-1000,1000))
  axis(4,c(-1000,10000),c(-1000,10000))
  points(xs[!pos],ys[!pos],pch=19,cex=bubwidth[!pos])
  if(leg){
    legend('topleft',legend=c("neg","pos"),pch=c(19,1),cex=0.8,bty='n')
    levs<-c(1,2)
    legend('topright',legend=levs,pch=rep(1,2),pt.cex=(levs*scaler)^0.5,cex=0.8,bty='n')
    
  }
}


Myplot<-function(filey,col="black"){
  Fit<-readRDS(filey)
  M<-Fit$rep$MAA[,1]
  plot(Fit$years,M,ylim=c(0,max(M)*1.01),xlab="Year",yaxs='i',type='l',col=col)
  abline(h=(0:10)/10,col='grey')
  abline(v=seq(1900,2050,by=5),col='grey')
  lines(Fit$years,M,col=col)
}



REMyplot<-function(filey,col="black"){
  Fit<-readRDS(filey)
  M<-exp(Fit$rep$M_a+Fit$rep$M_re[,1])
  plot(Fit$years,M,ylim=c(0,max(M)*1.01),xlab="Year",yaxs='i',type='l')
  abline(h=(0:10)/10,col='grey')
  abline(v=seq(1900,2050,by=5),col='grey')
  lines(Fit$years,M,col=col)
}


REMaplot<-function(filey,col="black"){
  Fit<-readRDS(filey)
  M<-exp(Fit$rep$M_a+Fit$rep$M_re[1,])
  plot(1:length(M),M,ylim=c(0,max(M)*1.01),xlab="Age",yaxs='i',type='l')
  abline(h=(0:10)/10,col='grey')
  abline(v=seq(1900,2050,by=5),col='grey')
  lines(1:length(M),M,col=col)
}


REMy_comp_plot<-function(mdirs,runnams=NULL,REapprox=NULL,cols=c("black","red","green","blue","orange","grey","purple","pink","darkred","darkblue","darkgreen","azure","darkgrey")){
  
  if(is.null(REapprox))REapprox<-rep(F,length(mdirs))
  if(is.null(runnams))runnams=paste("Model",1:length(mdirs))
  coly<-makeTransparent(cols)
  Fit<-readRDS(mdirs[2])
  
  getMy<-function(x,mdirs,REapprox){
    Fit<-readRDS(mdirs[x])
    if(REapprox[x]){
      RE<-Fit$rep$NAA_devs
      impliedM <- (-RE)
      annualized<-c(NA,apply(impliedM[,3:9],1,mean))
      return(annualized+Fit$rep$MAA[,2])
    }else{  
      return(Fit$rep$MAA[,2])
    }
  }
  
  Ms<-sapply(1:length(mdirs),getMy,mdirs=mdirs,REapprox=REapprox)
  layout(matrix(c(1,2),byrow=T,nrow=1),widths=c(1,0.7))
  par(mai=c(0.4,0.65,0.05,0.05), omi=c(0.3,0.05,0.05,0.05))
  
  matplot(Fit$years,Ms,ylim=c(0,max(Ms,na.rm=T)*1.01),ylab="",xlab="Year",yaxs='i',type='l',col=cols,lty=1,lwd=2)
  mtext("Annual Instantaneous Natural Mortality Rate (M)",2,line=2.5)
  plot(1,1,axes=F,xlab="",ylab="",col="white"); legend('left',legend=runnams, text.col=cols,bty='n',cex=1.2)
  
  
}



SelDens<-function(Fit,nams=c("Fleet","NMFS Spr","NMFS Fall","DFO Spr")){
  
  sel0<-Fit$rep$selAA
  pt<-Fit$input$data$selblock_pointer_fleets[,1]
  sel0[[1]][pt==2,]<-sel0[[2]][pt==2,]
  sel<-sel0[c(1,3,4,5)]
  
  nsel<-length(sel)
  ny<-nrow(sel[[1]])
  na<-ncol(sel[[1]])
  years<-Fit$years
  dat<-data.frame(Year=rep(years,na*nsel),
             Age=rep(rep(1:na,each=ny),nsel),
             Selnam=rep(nams,each=na*ny),
             Selectivity=unlist(sel))
             
  ggplot(dat,aes(Age,Year,fill=Selectivity)) +
    geom_tile()+
    facet_wrap(~Selnam)+
    scale_y_reverse()
 
}



SelDenslist<-function(Fit,Fit2,runnams=c("iid","2DAR1"),nams=c("Fleet","NMFS Spr","NMFS Fall","DFO Spr"),size=2.5){
  
  sel0<-Fit$rep$selAA
  pt<-Fit$input$data$selblock_pointer_fleets[,1]
  sel0[[1]][pt==2,]<-sel0[[2]][pt==2,]
  sel<-sel0[c(1,3,4,5)]
  
  nsel<-length(sel)
  ny<-nrow(sel[[1]])
  na<-ncol(sel[[1]])
  years<-Fit$years
  dat<-data.frame(Year=rep(years,na*nsel),
                  Age=rep(rep(1:na,each=ny),nsel),
                  Selnam=rep(nams,each=na*ny),
                  Selectivity=unlist(sel),
                  Run=rep(runnams[1],ny*na*nsel),
                  text=round(unlist(sel),2))
  
  Fit<-Fit2
  sel0<-Fit$rep$selAA
  pt<-Fit$input$data$selblock_pointer_fleets[,1]
  sel0[[1]][pt==2,]<-sel0[[2]][pt==2,]
  sel<-sel0[c(1,3,4,5)]
  
  nsel<-length(sel)
  ny<-nrow(sel[[1]])
  na<-ncol(sel[[1]])
  years<-Fit$years
  dat2<-data.frame(Year=rep(years,na*nsel),
                  Age=rep(rep(1:na,each=ny),nsel),
                  Selnam=rep(nams,each=na*ny),
                  Selectivity=unlist(sel),
                  Run=rep(runnams[2],ny*na*nsel),
                  text=round(unlist(sel),2))
  
  dat3<-rbind(dat,dat2)
  dat3<-dat3[dat3$Selnam=="Fleet",]
  
  ggplot(dat3,aes(Age,Year,fill=Selectivity)) +
    geom_tile()+
    facet_wrap(~Run)+
    scale_y_reverse()+
    geom_text(aes(y=Year,x=Age,label=text),size=size)
  
}


PlotNAA_RE<-function(Fit,as_mult=F,inc_age_1=F){
  
  res<-Fit$rep$NAA_devs

  years<-Fit$years[2:length(Fit$years)]
  ny<-length(years)
  na<-ncol(res)
  ages<-1:na
  if(!inc_age_1){
    res=res[,2:na]
    na<-na-1
    ages<-2:(na+1)
  }
  
  res<-as.vector(res)
  
  if(as_mult)res=exp(res)-1
  
  dat<-data.frame(Year=rep(years,na),
                  Age=rep(ages,each=ny),
                  NAA_RE=res)
  
  ggplot(dat) +
    geom_tile(aes(x=Year,y=Age,fill=NAA_RE))+
    #scale_y_reverse()+
    colorspace::scale_fill_continuous_diverging()
  
   
  
}

dosuv<-function(resmat){
  
  na<-ncol(resmat)
  ny<-nrow(resmat)
  surv<-array(NA,dim=dim(resmat))

  for(a in na:1)    surv[1,a]<-exp(resmat[1,a])
  for(y in 2:ny)    surv[y,1]<-exp(resmat[y,1])
  
  for(y in 2:ny){
    
    for(a in na:2){
      surv[y,a]<-exp(resmat[y,a])*surv[y-1,a-1]
    }
  }
  
  surv
  
}

PlotNAA_RElist<-function(Fitlist,nams=NULL,inc_age_1=F,addmult=F,is.surv=F){
  
  res<-lapply(Fitlist,function(x)x$NAA_devs)
 
  na<-ncol(res[[1]])
  
  if(!inc_age_1){
    res<-lapply(res,function(x)x[,2:ncol(x)])
    ages<-2:na
    na<-na-1
  }else{
    ages<-1:na
  }
  surv<-lapply(res,dosuv)
  
  ncomp<-length(Fitlist)
  if(is.null(nams))nams<-paste("Model",1:ncomp)
   
  getyr<-function(Fit){
     years<-Fit$years[2:length(Fit$years)]
     rep(years,na)
  } 
   
  geta<-function(Fit){
     years<-Fit$years[2:length(Fit$years)]
     rep(ages,each=length(years))
  }  
  
  Year1<-lapply(Fitlist,getyr)
  Names<-rep(nams,unlist(lapply(Year1,function(x)length(x))))
  Year<-unlist(Year1)
  Age<-unlist(lapply(Fitlist,geta))
  res<-unlist(res)
  nr<-ceiling(ncomp^0.5)
  multi<-round(exp(res),2)
  if(is.surv)multi<-round(unlist(surv),2)
  
  getstats<-function(x){
    inds<-log(x$ind_paa/aperm(x$pred_index_paa,c(2,1,3)))
    inds[inds==Inf|inds==-Inf]<-NA
    cs<-log(x$catch_paa[1,,]/x$pred_catch_paa[,1,])
    cs[cs==Inf|cs==-Inf]<-NA
    paste(c("R","C","S","F","D"),round(c(sd(x$NAA_devs[,2:9]),sd(cs,na.rm=T),apply(inds,1,sd,na.rm=T))*100,0),sep=": ")
  }
  
  tdat<-data.frame(Age=rep(seq(2.2,8.8,length.out=5),length(Fitlist)),
                   Year=rep(min(Year)-2,length(Fitlist)*5),
                   Names=rep(nams,each=5))
  
  txt=unlist(lapply(Fitlist,getstats))
  
  dat<-data.frame(Year=Year,
                  Age=Age,
                  Names=Names,
                  NAA_RE=res,
                  multi=multi)
  
  if(addmult){
  ploty<-ggplot(dat) +
    geom_tile(aes(y=Year,x=Age,fill=NAA_RE))+
      scale_y_reverse()+
      scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))+
      facet_wrap(~Names, nrow=1)+
      geom_text(aes(y=Year,x=Age,label=multi),size=2.5)+
      colorspace::scale_fill_continuous_diverging()+
      geom_text(data = tdat,label = txt, aes(y=Year,x=Age),size=3)
  }else{
    ploty<-ggplot(dat) +
      geom_tile(aes(y=Year,x=Age,fill=NAA_RE))+
      scale_y_reverse()+
      scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))+
      facet_wrap(~Names, nrow=1)+
      colorspace::scale_fill_continuous_diverging()+
      geom_text(data = tdat,label = txt, aes(y=Year,x=Age),size=3)
    
  }
  #if(tempy){
   # ploty<-ggplot(dat) +
    #  geom_tile(aes(x=Year,y=Age,fill=NAA_RE))+
      #scale_y_reverse()+
     # scale_y_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))+
      #facet_wrap(~Names, nrow=nr)+
      #colorspace::scale_fill_continuous_diverging()+
      #geom_text(data = tdat,label = txt, aes(y=Year,x=Age),size=3)
    
  #}
  ploty
  
}

PlotNAA_both<-function(Fit){
  
  res<-Fit$rep$NAA_devs[,2:ncol(rep$NAA_devs)]
  
  na<-ncol(res)
  ages<-(1:na)+1
  
  surv<-round(dosuv(res),2)
  
  years<-Fit$years[2:length(Fit$years)]
  Year<-rep(years,na)
  Age<-rep(ages,each=length(years))
  
  multi<-round(exp(unlist(res)),2)
  
 # tdat<-data.frame(Age=rep(seq(2.2,8.8,length.out=5),length(Fitlist)),
  #                 Year=rep(min(Year)-2,length(Fitlist)*5))
                   
  
 # txt=getstats(rep)
  
  dat1<-data.frame(Year=Year,
                  Age=Age,
                  NAA_RE=as.vector(res),
                  multi=as.vector(multi))
  
  
  p1<-ggplot(dat1) +
    geom_tile(aes(y=Year,x=Age,fill=NAA_RE))+
    scale_y_reverse()+
    scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))+
    geom_text(aes(y=Year,x=Age,label=multi),size=2.5)+
    colorspace::scale_fill_continuous_diverging()
  
  dat2<-data.frame(Year=Year,
                   Age=Age,
                   NAA_RE=as.vector(res),
                   multi=as.vector(surv))
  
  p2<-ggplot(dat2) +
    geom_tile(aes(y=Year,x=Age,fill=NAA_RE))+
    scale_y_reverse()+
    scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))+
    geom_text(aes(y=Year,x=Age,label=multi),size=2.5)+
    colorspace::scale_fill_continuous_diverging()
  
     
  cowplot::plot_grid(p1,p2,labels=c("(a) Multiplier","(b) Natural Survival (Down Cohort, Excl. M)"),label_size = 9) 

  
}





cohortrat<-function(mat,inds=1:4){
  
  mat<-mat[,inds]
  nc<-ncol(mat)
  ny<-nrow(mat)
  nyc<-ny-nc+1
  rat<-array(NA,c(nyc,nc))
  for(y in 1:nyc){
    for(c in 1:nc){
      rat[y,c]<-mat[y-1+c,c]/mat[y,1]
    }
    
  }
  rat 
  
}



age1fracplot<-function(dat,Ilabs=c("NMFS Spr","NMFS Fall","DFO Spr"),ystart=10,log=F){
  
  cols=c('red','green','blue')
  dats<-dat
  years<-dats$year1+(0:(dats$n_years-1))
  ni<-length(dats$IAA_mats)-1
  nplot<-ni+1
  ncol<-ceiling(nplot^0.5)
  nrow=ceiling(nplot/ncol)
  
  par(mfrow=c(nrow,ncol),mai=c(0.4,0.55,0.3,0.05), omi=c(0.3,0.25,0.05,0.05))

  Cfracs<-cohortrat(dats$CAA_mats[[1]])
  Cfracs[Cfracs=='Inf']<-NA
  Ifracs<-lapply(dats$IAA_mats,cohortrat,inds=4:7)
  
  ylims=c(0,max(unlist(Cfracs),unlist(Ifracs),na.rm=T))
  
  
  yind<-c(ystart:nrow(Cfracs))
  if(!log){
    matplot(years[yind],Cfracs[yind,2:4],col=cols,type="l",lty=1,ylab="")
  }else{
    matplot(years[yind],log(Cfracs[yind,2:4]),col=cols,type="l",lty=1,ylab="")
  }
  abline(v=2014,col='grey')
  mtext("Fleet",line=0.5)
  
  for(i in 1:ni){
     if(!log){
       matplot(years[yind],Ifracs[[i]][yind,2:4],col=cols,type="l",lty=1,ylab="")
     }else{
       matplot(years[yind],log(Ifracs[[i]][yind,2:4]),col=cols,type="l",lty=1,ylab="")
     }
       mtext(Ilabs[i],line=0.5)
     abline(v=2014,col='grey')
  }
  
  legend('topleft',legend=paste("age",2:4,"/ age 1"),bty='n',text.col=cols)
  mtext('Ratio of age composition, within cohort',2,outer=T)
  mtext("Year (cohort is year - 1)",1,outer=T)

}


cohortrat2<-function(mat,inds=1:4){
  
  mat<-mat[,inds]
  nc<-ncol(mat)
  ny<-nrow(mat)
  nyc<-ny-nc+1
  rat<-array(NA,c(nyc,nc))
  for(y in 1:nyc){
    for(c in 1:nc){
      rat[y,c]<-mat[y-1+c,c]/mat[y,2]
    }
    
  }
  rat 
  
}



age2fracplot<-function(dat,Ilabs=c("NMFS Spr","NMFS Fall","DFO Spr"),ystart=10,log=F){
  
  cols=c('green','blue')
  dats<-dat$dat
  years<-dats$year1+(0:(dats$n_years-1))
  ni<-length(dats$IAA_mats)-1
  nplot<-ni+1
  ncol<-ceiling(nplot^0.5)
  nrow=ceiling(nplot/ncol)
  
  par(mfrow=c(nrow,ncol),mai=c(0.4,0.55,0.3,0.05), omi=c(0.3,0.25,0.05,0.05))
  
  Cfracs<-cohortrat2(dats$CAA_mats[[1]])
  Cfracs[Cfracs=='Inf']<-NA
  Ifracs<-lapply(dats$IAA_mats,cohortrat2,inds=4:7)
  
  ylims=c(0,max(unlist(Cfracs),unlist(Ifracs),na.rm=T))
  
  yind<-c(ystart:nrow(Cfracs))
  if(!log){
    matplot(years[yind],Cfracs[yind,3:4],col=cols,type="l",lty=1,ylab="")
  }else{
    matplot(years[yind],log(Cfracs[yind,3:4]),col=cols,type="l",lty=1,ylab="")
  }
  
  abline(v=2014,col='grey')
  mtext("Fleet",line=0.5)
  
  for(i in 1:ni){
    if(!log){
      matplot(years[yind],Ifracs[[i]][yind,3:4],col=cols,type="l",lty=1,ylab="")
    }else{
      matplot(years[yind],log(Ifracs[[i]][yind,3:4]),col=cols,type="l",lty=1,ylab="")
    }
    
    mtext(Ilabs[i],line=0.5)
    abline(v=2014,col='grey')
  }
  
  legend('topleft',legend=paste("age",3:4,"/ age 2"),bty='n',text.col=cols)
  mtext('Ratio of age composition, within cohort',2,outer=T)
  mtext("Year (cohort is year - 1)",1,outer=T)
  
  
}





print("WHAM/Source/Extra_plotting.R loaded")
