# ==============================================================================================================================================
# Tentative, Preliminary, Exploratory TRAC Stock Assessment of Haddock in EGB area using
# Woods Hole Assessment Model (WHAM) (Miller and Stock 2021, https://timjmiller.github.io/wham/articles)
# ==============================================================================================================================================

# Misc functions for extracting and comparing WHAM models

# Tom Carruthers
# March 2021


sim_WHAM<-function(Fit){

  out<-list()
  est<-Fit$simulate
  out$NAA <- est()$NAA
  out$BAA <- out$NAA*Fit$input$data$waa[2,,]
  out$FAA <- est()$FAA_tot
  out$SSB <- est()$SSB
  out$Fbar <- est()$Fbar

  out

}


extract_WHAM<-function(Fit){

  out<-list()
  out$model_name<-Fit$model_name
  # Estimates
  est<-Fit$rep
  out$NAA_devs<-est$NAA_devs
  out$NAA <- est$NAA
  out$BAA <- out$NAA*Fit$input$data$waa[2,,]
  out$FAA <- est$FAA
  out$QAA <- est$QAA
  out$MAA <- est$MAA
  out$selAA <- est$selAA
  out$FAA_tot <- est$FAA_tot
  out$SSB <- est$SSB
  out$Fbar <- est$Fbar
  out$pred_NAA <- est$pred_NAA
  out$pred_indices <- est$pred_indices
  out$q<-est$q
  out$pred_catch <- est$pred_catch
  out$pred_catch_paa <- est$pred_catch_paa
  out$pred_index_paa <-  est$pred_index_paa
  out$nll <- c(Obj = Fit$opt$obj, Total = est$nll, indices = apply(est$nll_agg_indices,2,sum),
               index_comp = apply(est$nll_index_acomp,2,sum), catch = sum(est$nll_agg_catch),
               catch_comp = sum(est$nll_catch_acomp),
               k= length(Fit$opt$par),
               nll_sel = est$nll_sel, nll_M = est$nll_M, nll_Ecov = est$nll_Ecov, nll_NAA = est$nll_NAA)

  out$npar <- length(Fit$opt$par)
  out$AIC  <- 2 * (Fit$opt$obj + out$npar) # WHAM: k = length(x$opt$par);  2*(x$opt$obj + k) # AIC
  
  # Retro
  if(!is.null(Fit$peels))out$rho<-mohns_rho(Fit)

  # Data
  dat <- Fit$input$data
  out$agg_indices <- dat$agg_indices
  out$catch_paa <- dat$catch_paa
  out$ind_paa <- dat$index_paa
  out$catch_Neff <- dat$catch_Neff
  out$index_Neff <- dat$index_Neff
  out$fracyr_indices <- dat$fracyr_indices
  out$waa <- dat$waa
  out$waa_pointer_indices <- dat$waa_pointer_indices
  out$waa_pointer_fleets <- dat$waa_pointer_fleets
  
  out$conv <- Fit$opt$convergence
  #if(!is.null(dat$agg_catch)){
    out$catch <- dat$agg_catch[,1]
  #}else{
  #  out$catch <- dat$CAA

  #}

  # Derived quantities
  out$FXSPR<-exp(Fit$rep$log_FXSPR)

  # Indexing and labels
  out$years <- Fit$years
  out

}

extract_WHAM_list<-function(mdirs){

  out<-list()
  #mods<-list()
  for(i in 1:length(mdirs)){
    Fit<-readRDS(mdirs[[i]])
    out[[i]]<-extract_WHAM(Fit)
  }
  #res<-compare_wham_models(mods)
  names(out)<-names(mdirs)

  out

}

j1<-function(filenam)jpeg(filenam,res=400,units='in',width=4,height=4)
j2<-function(filenam)jpeg(filenam,res=400,units='in',width=8,height=4)
j4<-function(filenam)jpeg(filenam,res=400,units='in',width=8,height=8)
j8<-function(filenam)jpeg(filenam,res=400,units='in',width=8,height=10)


countRE<-function(Fit){
  n=0
  is_N_RE<-Fit$input$data$n_NAA_sigma > 1# !all(Fit$parList$log_NAA[,2:ncol(Fit$parList$log_NAA)]==10) # are there N_RE
  if(is_N_RE) n=n+sum(Fit$parList$log_NAA!=10)-1
  if(!all(Fit$parList$selpars_re==0))n=n+sum(Fit$parList$selpars_re!=0)-1
  if(!all(Fit$parList$M_re==0))n=n+length(unique(as.vector(Fit$parList$M_re)))-1
  n
}

getNpar<-function(mdirs,nams=NULL){
  if(is.null(nams))nams=paste("Model",1:length(mdirs))
  Fits<-lapply(mdirs,function(x)readRDS(x))
  Fixed=sapply(Fits,function(x)length(x$opt$par))
  RE=sapply(Fits,function(x)countRE(x))
  data.frame(Model=nams,Fixed=Fixed,RE=RE,Total=Fixed+RE)
}

getDiags<-function(mdirs,nams=NULL){
  
  sums<-extract_WHAM_list(mdirs)
  isretro<-!is.null(listo[[1]]$rho)
  if(isretro){
    rho_SSB<- sapply(listo,function(x)x$rho[1])
    rho_Fbar<-sapply(listo,function(x)x$rho[2])
  }
  AIC<-sapply(listo,function(x)x$AIC)
  dAIC<-AIC-min(AIC)
  totLH<-sapply(listo,function(x)x$nll[1])
  if(is.null(nams))nams=paste("Model",1:length(mdirs))
  if(isretro)return(data.frame(Model=nams,AIC=AIC,dAIC=dAIC,nll=totLH,rho_SSB=rho_SSB,rho_Fbar=rho_Fbar))
  if(!isretro)return(data.frame(Model=nams,AIC=AIC,dAIC=dAIC,nll=totLH))
  
}

print('Misc functions loaded /WHAM/Source/Misc.R')
