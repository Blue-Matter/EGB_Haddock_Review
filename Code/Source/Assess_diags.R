# Diagnostics for assessment

make_jit_list<-function(dat,n=5,cv=0.1,seed=1){

  datlist<-list()
  set.seed(seed)
  for(jj in 1:n){
    in_J<-dat
    in_J$par$log_N1_pars[] <-in_J$par$log_N1_pars[] + rnorm(length(in_J$par$log_N1_pars),0,cv)
    in_J$par$log_NAA[,1] <-in_J$par$log_NAA[,1] + rnorm(length( in_J$par$log_NAA[,1]),0,cv)
    datlist[[jj]]<-in_J
  }

  datlist

}



make_M_list<-function(dat,n=11){

  Mults<- seq(-0.4,0.4,length.out=n)
  datlist<-list()
  for(jj in 1:n){
    in_J<-dat
    in_J$par$M_a <- in_J$par$M_a+Mults[jj]
    datlist[[jj]]<-in_J
  }

  datlist

}



make_Mstep_list<-function(Mod,n=11,ny=7){
  
  Mults<- seq(-0.7,0.7,length.out=n)
  nyears<-length(Mod$years)
  ind<-nyears-((ny-1):0)
  datlist<-list()
  for(jj in 1:n){
    in_J<-Mod
    in_J$par$M_re[ind,] <- Mults[jj]
    datlist[[jj]]<-in_J
  }

  datlist
  
}





