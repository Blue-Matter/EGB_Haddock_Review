
dofit<-function(in_obj,run_name=NULL,do.sdrep=T, do.retro=T, do.osa=F, do.proj=F, save.sdrep=F, do.check=T, outdir=NULL, silent=T,out.type="html",ext=NULL){

  if(is.null(ext))ext="/WHAM/EGB_model_fits/"
  set.seed(1)
  if(is.null(outdir)){
    outfile=paste0(getwd(),ext,run_name,"/Fit.rda")
    outdir=paste0(getwd(),ext,run_name)
  }else{
    outfile=paste0(outdir,"/Fit.rda")
  }

  out <- fit_wham(in_obj, do.sdrep=do.sdrep, do.retro=do.retro, do.osa=do.osa, do.proj=do.proj, save.sdrep=save.sdrep,
                  do.check=do.check, MakeADFun.silent = silent, retro.silent = silent)

  if(grepl("_RESCALE",in_obj$model_name)){
    out$check_convergence <- conv_info(mod=out)
    # rescale selectivities and q
    out<- rescale.out(out)
  }
  
  check_convergence(out)
  out

}

# OM needs to do rescale after sampling
dofitOM<-function(in_obj,run_name=NULL,do.sdrep=T, do.retro=T, do.osa=F, do.proj=F, save.sdrep=F, do.check=T, outdir=NULL, silent=T,out.type="html",ext=NULL){
  
  
  if(is.null(ext))ext="/WHAM/EGB_model_fits/"
  set.seed(1)
  if(is.null(outdir)){
    outfile=paste0(getwd(),ext,run_name,"/Fit.rda")
    outdir=paste0(getwd(),ext,run_name)
  }else{
    outfile=paste0(outdir,"/Fit.rda")
  }
  
  out<-fit_wham(in_obj, do.sdrep=do.sdrep, do.retro=do.retro, do.osa=do.osa, do.proj=do.proj, save.sdrep=save.sdrep,
                  do.check=do.check, MakeADFun.silent = silent, retro.silent = silent)
  
  out$check_convergence <- conv_info(mod=out)
  check_convergence(out)
  out
  
}



doRep<-function(dirs,out.type="html",res=72){

  for(i in 1:length(dirs)){
    outfile=paste0(dirs[i],"/Fit.rda")
    Fit<-readRDS(outfile)
    plot_wham_output(mod=Fit, out.type=out.type,res=res,dir.main=dirs[i])
    print(paste("Done:",dirs[i]))
  }
 
}

plot_wham_list<-function(whamdirs){
  for(i in 1:length(whamdirs)){
    out<-readRDS(paste0(whamdirs[i],"/Fit.rda"))
    plot_wham_output(mod=out, out.type='html',res=72,dir.main="C:/temp")
  }
}
