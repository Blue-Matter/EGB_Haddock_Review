# parallel computing scripts

WHAMcluster<-function(){
  
  library(snowfall)
  library(parallel)
  sfInit(parallel=T,cpus=detectCores()/2)
  sfLibrary("wham", character.only = TRUE, verbose = FALSE)
  sfExport(list=list("conv_info","rescale.out"))
  
}


fit_wham_parallel<-function(x ,datlist, dirs=NULL, do.osa=F, do.retro=F, do.sdrep=F, do.proj=F, save.sdrep=F){
  
  out<- fit_wham(datlist[[x]], do.osa = do.osa, do.retro = do.retro, 
                  do.sdrep = do.sdrep, do.proj = do.proj, 
                  save.sdrep = save.sdrep,  MakeADFun.silent = T,
                  retro.silent = T, osa.opts = list(method = "oneStepGeneric",
                                                     parallel = FALSE)) 
  
  if(grepl("_RESCALE",datlist[[x]]$model_name)){
    out$check_convergence <- conv_info(mod=out)
    
    # rescale selectivities and q
    out<- rescale.out(out)
    # out.rescale$model_name = "f1e.free.rescale"
  }
  
  if(!is.null(dirs))  saveRDS(out,dirs[x])
  
  out
  
}