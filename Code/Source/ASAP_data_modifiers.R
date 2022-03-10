# Input object manipulation function

IndW_adj<-function(asap,fleetno=2,defCV=0.2,plot=T){
  
  nams=c("NMFS Spr","NMFS Fall","DFO Spr")
  years<-asap$dat$year1+(1:asap$dat$n_years)-1
  asap2<-asap # make copy for debugging
  ni<-length(asap$dat$IAA_mats)
  na<-asap$dat$n_ages
  aind<-3+1:na
  if(plot)par(mfrow=c(2,2),mai=c(0.5,0.5,0.1,0.1))
  
  for(i in 1:length(nams)){
    mat<-asap$dat$IAA_mats[[i]]
    index1<-mat[,2]
    mu1<-mean(index1[index1!=0])
    WAA<-asap$dat$WAA_mats[[asap$dat$WAA_pointers[fleetno+i]]]
    index2<-apply(mat[,aind]*WAA,1,sum)
    mu2<-mean(index2[index2!=0])
    index2<-index2*mu1/mu2 # give same scale as before
    mat[,2]<-index2 #index
    mat[index1!=0,3]<-defCV
    if(plot){
      ind1<-index1
      ind1[ind1==0]<-NA
      ind2<-index2
      ind2[ind2==0]<-NA
      hl<-(0:100)*100# pretty(seq(0,max(ind1,ind2,na.rm=T)*2,length.out=40))
      matplot(years,cbind(ind1,ind2),ylim=c(0,max(ind1,ind2,na.rm=T)*1.02),type='l',col=c("black","red"),lty=1,yaxs='i')
      abline(h=hl,col='light grey')
      matplot(years,cbind(ind1,ind2),type='l',col=c("black","red"),lty=1,add=T)
      
      matplot(years,cbind(ind1,ind2),col=c("black","red"),pch=1,cex=0.9,add=T)
      if(i==1)legend('left',legend=c("Numbers","Weight (rescaled)"),text.col=c("black","red"),bty='n')
      legend('top',legend=nams[i],bty="n")
      }
    asap2$dat$IAA_mats[[i]]<-mat
  }
  asap2$dat$index_units[]<-1
  asap2
  
}

selb1_adj<-function(asap){
  asap$comments<-paste0(asap$comments,"_data_1selb")
  
  asap$dat$sel_block_assign[[1]][]<-1
  asap$dat$sel_block_option<-2  # 1 age specific 2 logistic
  asap$dat$n_fleet_sel_blocks<-1
  asap
}
FALLig_adj<-function(asap){
  asap$comments<-paste0(asap$comments,"_ignore_Fall_index")
  ESScol<-ncol(asap$dat$IAA_mats[[2]])

  asap$dat$IAA_mats[[2]][,3]<-10
  asap$dat$IAA_mats[[2]][,ESScol]<-0.01
  asap
}

ind_prec_adj<-function(asap,ind=1,adj=2){
  asap$comments<-paste0(asap$comments,"index ", ind," precinc = ",adj)
  
  for(i in ind){
    CV<-asap$dat$IAA_mats[[i]][,3]
    prec<-1/CV^2
    nuCV=sqrt(1/(prec*adj))
    asap$dat$IAA_mats[[i]][,3]<-nuCV
  }
  
  asap
}


IESS_adj<-function(asap,iwt=1){

  ni<-asap$dat$n_indices

  if(length(iwt)==1)Inds<-rep(iwt,ni)

  for(i in 1:ni){
    mat<-asap$dat$IAA_mats[[i]]
    ESScol<-ncol(mat)
    ESS<-mat[,ESScol]
    asap$dat$IAA_mats[[i]][,ESScol]<-ESS*Inds[i]
  }

  asap$comments<-paste0(asap$comments,"_IESS",iwt)
  asap

}

FlatInd_adj<-function(asap,frac=0){
  
  ni<-asap$dat$n_indices
  
  for(i in 1:ni){
    mat<-asap$dat$IAA_mats[[i]]
    Ind<-mat[,2]
    muI<-mean(Ind,na.rm=T)
    res<-log(Ind/muI)
    nInd<-muI*exp(res*frac)
    nInd[is.na(nInd)]<-0
    
    asap$dat$IAA_mats[[i]][,2]<-nInd
    #plot(Ind,type='l',col='red');lines(nInd,col='blue')
  }
  
  asap$comments<-paste0(asap$comments,"_flattenInd_frac=",frac)
  asap
  
}

CESS_adj<-function(asap,wt=1){

  asap$dat$catch_Neff[]<-asap$dat$catch_Neff[]*wt
  asap$comments<-paste0(asap$comments,"_CESS",wt)
  asap

}


ESS_adj<-function(asap,iwt=1){

  ni<-asap$dat$n_indices

  if(length(iwt)==1)Inds<-rep(iwt,ni)

  for(i in 1:ni){
    mat<-asap$dat$IAA_mats[[i]]
    ESScol<-ncol(mat)
    ESS<-mat[,ESScol]
    asap$dat$IAA_mats[[i]][,ESScol]<-ESS*Inds[i]
  }
  asap$dat$catch_Neff[]<-asap$dat$catch_Neff[]*wt
  asap$comments<-paste0(asap$comments,"_ESS",iwt)
  asap

}


DropIndF_adj<-function(asap){

  dims<-dim(asap$dat$IAA_mats[[2]])
  asap$dat$IAA_mats[[2]][dims[1],c(2:3,13)]<-0
  asap$comments<-paste0(asap$comments,"_nolastFall")
  asap

}

M_adj<-function(asap, M=0.2){
  asap$dat$M[]<-M
  asap$comments<-paste0(asap$comments,"_M = ",M)
  asap
}


MLorz<-function(asap){
  Mmu<-mean(asap$dat$M)
  Mvec<-1/apply(asap$dat$WAA_mats[[1]],2,mean)
  M0<-Mvec/mean(Mvec)*Mmu
  asap$dat$M[]<-rep(M0,each=asap$dat$n_years)
  asap$comments<-paste0(asap$comments,"_LorenzM")
  asap
}

FSBlockadj<-function(asap,yrs=c(1976,1994),types=c(1,1,1),fno=1){

  if(length(yrs)>0){

    if(length(types)==1) types=rep(types,length(yrs)+1)
    years<-asap$dat$year1+(0:(asap$dat$n_years-1))
    asap$dat$n_fleet_sel_blocks<-nb<-length(yrs)+1
    bpoints<-c(min(years),yrs,max(years))

    intl<-NULL
    for(bb in 1:nb){

      if(bb<nb)intl<-c(intl, rep(bb,bpoints[bb+1]-bpoints[bb]))
      if(bb==nb)intl<-c(intl, rep(bb,bpoints[bb+1]-bpoints[bb]+1))

    }

    asap$dat$sel_block_assign[[fno]]<-intl
    asap$dat$sel_block_option <- types # 1 appears age specific (non parametric) 2 appears logistic, 3 is double logistic

  }else{

    asap$dat$n_fleet_sel_blocks<-1
    asap$dat$sel_block_assign[[fno]][]<-1
    asap$dat$sel_block_option <- types[1]

  }

  asap

}

sR_adj<-function(asap,sR=1.0){

  asap$dat$recruit_cv<-rep(sR,asap$dat$n_years)
  asap$comments<-paste0(asap$comments,"_sigR",sR)
  asap

}

Mrecent<-function(asap,ny=7,M=0.3){
  
  asap$comments<-paste0(asap$comments,"_Mrecent_M=",M,"recyrs=",ny)
  ind<-asap$dat$n_years-((ny-1):0)
  asap$dat$M[ind,]<-M
  asap

}  

no_age_1<-function(asap,aind=4:12){
  
  dat<-dat2<-asap$dat
  #aind<-grep("Age",colnames(dat2$IAA_mats[[1]]))
  na<-length(aind)
  for(i in 1:length(dat2$IAA_mats)){
    dat$IAA_mats[[i]][,aind[1]]<-0
    tot<-apply(dat2$IAA_mats[[i]][,aind],1,sum)
    p1<-dat2$IAA_mats[[i]][,aind[1]]
    div<-(tot-p1)/tot
    div[is.na(div)]<-1
    #dat$IAA_mats[[i]][,aind[2:na]]<-dat2$IAA_mats[[i]][,aind[2:na]]/div
    dat$IAA_mats[[i]][,2]<-dat$IAA_mats[[i]][,2]*div
  }
  
  dat$CAA_mats[[1]][,1]<-0
  
  asap$dat<-dat
  asap
}
 
# turn off age 1s for select surveys and the fleet
no1sel<-function(asap,asap_no1,surv=NULL,fleet=F){
  new_asap<-asap
  if(fleet)new_asap$dat$CAA_mats[[1]][,1]<-0
  if(!is.na(surv[1])){ # if an survey number is provided
    for(i in surv)new_asap$dat$IAA_mats[[i]]<-asap_no1$dat$IAA_mats[[i]]
  }
  new_asap
}

# turn off age 1s for select surveys and the fleet
no1sel_yspec<-function(asap,asap_no1,yspec=2014,surv=1:3,fleet=TRUE){
  years<-asap$dat$year1+(0:(asap$dat$n_years-1))
  yind<-match(yspec,years)
  new_asap<-asap
  
  if(fleet)new_asap$dat$CAA_mats[[1]][yind,1]<-0
  if(!is.na(surv[1])){ # if an survey number is provided
    for(i in surv)new_asap$dat$IAA_mats[[i]][yind,]<-asap_no1$dat$IAA_mats[[i]][yind,]
  }
  new_asap
}

  
print("Input modifier functions loaded")
