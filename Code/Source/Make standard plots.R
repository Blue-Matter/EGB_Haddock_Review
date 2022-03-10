
# standard plots

doplots<-function(ext,sums, mdirs,runnams,ind=NA){
  pfnams<-paste0(getwd(),ext,"_",c("Comparison","ACompRes","Ind","IndRes","My","Sel","Rec_Ind","Rec_Comp","U"),".jpg")
  if(is.na(ind[1]))ind<-1:9
  
  if(1 %in%ind){
  jpeg(pfnams[1],res=500,units='in',width=9,height=7)
    Comp2(sums, namey=runnams,include_nll=T, include_conv=T)
  dev.off()
  }
  
  if(2 %in%ind){
  jpeg(pfnams[2],res=400,units='in',width=8,height=length(sums)*1.5)
    AAResComp(sums, YLabs=runnams,XLabs=c("Catch","NMFS Spr","NMFS Fall","DFO Spr"),Nscaled=F)
  dev.off()
  }

  if(3 %in%ind){
  jpeg(pfnams[3],res=400,units='in',width=7,height=8)
    IndComp(sums, YLabs=runnams,XLabs=c("NMFS Spr","NMFS Fall","DFO Spr"))
  dev.off()
  }

  if(4 %in%ind){
  jpeg(pfnams[4],res=400,units='in',width=7,height=8)
    IndResComp(sums, YLabs=runnams,XLabs=c("NMFS Spr","NMFS Fall","DFO Spr"))
  dev.off()
  }

  if(5 %in%ind){
  jpeg(pfnams[5],res=400,units='in',width=6,height=4)
    REMy_comp_plot(mdirs,runnams)
  dev.off()
  }

  if(6 %in%ind){
  jpeg(pfnams[6],res=400,units='in',width=6,height=ceiling(length(mdirs)^0.5)*2.5)
    plotallsel(mdirs,runnams)
  dev.off()
  }

  if(7 %in%ind){
  jpeg(pfnams[7],res=400,units='in',width=6,height=ceiling(length(mdirs)^0.5)*2.5)
    Recent_Ind(sums,runnams)
  dev.off()
  }

  if(8 %in%ind){
  temp<-Recent_Comp(sums,runnams)
  ggsave(pfnams[8],temp)
  }

  if(9 %in%ind){ 
  jpeg(pfnams[9],res=400,units='in',width=6.5,height=length(sums)*1.5)
    Ucomp(sums, runnams)
  dev.off()
  }

}