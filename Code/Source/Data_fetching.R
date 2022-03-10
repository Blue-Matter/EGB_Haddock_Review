    

get_dat <- function(type="base"){
  
  if(type=="base"){
    asaptemp <- readRDS("Data/EGB_2021/Formatted/EGB_base.rda")
    asaptemp$comments<-"EGB Base"
  }else if(type=="drop_DFO"){
    asaptemp <- readRDS("Data/EGB_2021/Formatted/EGB_dropDFO.rda")
    asaptemp$comments<-"EGB drop_DFO"
  }else if(type=="no_deep"){
    asaptemp <- readRDS("Data/EGB_2021/Formatted/EGB_nodeep.rda")
    asaptemp$comments<-"EGB no_deep"
  }else if(type=="no_fall"){
    asaptemp <- readRDS("Data/EGB_2021/Formatted/EGB_nofall.rda")
    asaptemp$comments<-"EGB no_fall"
  }else if(type=="unblocked"){
    asaptemp <- readRDS("Data/EGB_2021/Formatted/EGB_unblocked.rda")
    asaptemp$comments<-"EGB unblocked_sel"
  }else if(type=="noAge1"){
    asaptemp <- readRDS("Data/EGB_2021/Formatted/EGB_noAge1.rda")
    asaptemp$comments<-"EGB noAge1"
  }
  
  asaptemp
  
}


get_datW <- function(type="base"){
  
  if(type=="base"){
    asap3 <- readRDS("Data/EGB+WGB/Formatted/EGB+WGB_base.rda")
  }else if(type=="unblocked"){
    asap3 <- readRDS("Data/EGB+WGB/Formatted/EGB+WGB_unblocked.rda")
  }else if(type=="drop_DFO"){
    asap3 <- readRDS("Data/EGB+WGB/Formatted/EGB+WGB_dropDFO.rda")
  }else if(type=="no_deep"){
    asap3 <- readRDS("Data/EGB+WGB/Formatted/EGB+WGB_nodeep.rda")
  }else if(type=="no_fall"){
    asap3 <- readRDS("Data/EGB+WGB/Formatted/EGB+WGB_nofall.rda")
  }
  asap3$comments<-"EGB+WGB Base"
  
  asap3
  
}


get_datGB <- function(type="base"){
  
  if(type=="base"){
    asap3 <- readRDS("Data/GB/Formatted/GB_base.rda")
    asap3$comments<-"GB Base"
  }else if(type=="no_fall"){
    asap3 <- readRDS("Data/GB/Formatted/GB_nofall.rda")
    asap3$comments<-"GB no_fall"
  }
  
  asap3
  
}
