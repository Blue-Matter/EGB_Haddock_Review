# set of functions for preparing wham input objects given varying specifications of natural mortality rate, random effects and likelihood functions

logit<-function(x)exp(x)/(1+exp(x))
invlogit<-function(x)log(x/(1-x))


# makes inits for selectivity parameters in the 'free' rescaling approach
make_Base_sel_init<-function(asap){
  
  years = asap$dat$year1 - 1 + 1:asap$dat$n_years 
  
  #fleet specs: logistic for both blocks
  fix_fleet_sel = lapply(1:2, function(x) NA)
  init_fleet_sel = lapply(1:2, function(x) c(2, 0.3))
  
  # get ASAP index selectivity specs (just take [1:3], since we are dropping index4=Y41) ====
  init_index_sel <- lapply(asap$dat$index_sel_ini, function(x) x[1:9,1]) [1:3]
  fix_index_sel <- lapply(asap$dat$index_sel_ini, function(x) which(x[1:9,2] < 0)) [1:3]
  
  # # fix spring index ages that hit boundary
  fix_index_sel[[1]] <- c(2, 4,5, 8)
  init_index_sel[[1]] <- rep(0.50, 9)
  init_index_sel[[1]][fix_index_sel[[1]]] <- 1.0
  
  # # fix dfo index ages that hit boundary
  fix_index_sel[[3]] <- c(3,9)
  init_index_sel[[3]] <- rep(0.50, 9)
  init_index_sel[[3]][fix_index_sel[[3]]] <- 1.0
 
  list( model=c(rep("logistic", 2), rep("age-specific", 3) ),
                    re=c( "none", "iid", rep("none", 3) ) ,
                    initial_pars=c(init_fleet_sel, init_index_sel) ,
                    fix_pars=rep(list(NA),5)) # 
  
}

make_noRE_sel_init<-function(asap){
  
  years = asap$dat$year1 - 1 + 1:asap$dat$n_years 
  
  #fleet specs: logistic for both blocks
  fix_fleet_sel = lapply(1:2, function(x) NA)
  init_fleet_sel = lapply(1:2, function(x) c(2, 0.3))
  
  # get ASAP index selectivity specs (just take [1:3], since we are dropping index4=Y41) ====
  init_index_sel <- lapply(asap$dat$index_sel_ini, function(x) x[1:9,1]) [1:3]
  fix_index_sel <- lapply(asap$dat$index_sel_ini, function(x) which(x[1:9,2] < 0)) [1:3]
  # # fix spring index ages that hit boundary
  fix_index_sel[[1]] <- c(2, 4,5, 8)
  init_index_sel[[1]] <- rep(0.50, 9)
  init_index_sel[[1]][fix_index_sel[[1]]] <- 1.0
  
  # # fix dfo index ages that hit boundary
  fix_index_sel[[3]] <- c(3,9)
  init_index_sel[[3]] <- rep(0.50, 9)
  init_index_sel[[3]][fix_index_sel[[3]]] <- 1.0
  
  list( model=c(rep("logistic", 2), rep("age-specific", 3) ),
        re=rep("none", 5),
        initial_pars=c(init_fleet_sel, init_index_sel) ,
        fix_pars=rep(list(NA),5)) # 
}

# ensures fixed parameters as NMFS have specified - this may change in make_Base_sel_init()
make_NMFS_sel_init<-function(asap){
  
  sel_list<-make_Base_sel_init(asap)
  fix_fleet_sel = lapply(1:2, function(x) NA)
  fix_index_sel <- lapply(asap$dat$index_sel_ini, function(x) which(x[1:9,2] < 0)) [1:3]
  fix_index_sel[[1]] <- c(2, 4, 5, 8)
  fix_index_sel[[3]] <- c(3,9)
  sel_list$fix_pars<-c(fix_fleet_sel, fix_index_sel)  
  sel_list
}


make_Mstep_ecov<-function(asap,ny){
  
  years = asap$dat$year1 - 1 + 1:asap$dat$n_years 
  
  ecov = list(
    label = "dummy",
    year = years,
    mean = matrix(rep(0:1, c(length(years)-ny, ny))), #dummy variable matrix
    logsigma = matrix(log(0.001), length(years), 1), #tiny observation error variance
    use_obs = matrix(1,length(years), 1), #no missing observations of the covariate
    where = "M", #affects M rather than recruitment, or catchability pars.
    how = 1, #turns on estimation of effect on M (rather than just model the environmental process without effects on the population)
    lag = 0, #covariate affects M in the same year
    process_model = "rw" #random walk model for the latent covariate (necessary, but not relevant)
  )
  ecov
}

prep_rescale<-function(input){
  dummy<-readRDS("Data/Dummy_inits/dummy.rds")
  input$map$logit_q = factor(rep(NA,input$data$n_indices)) # don't estimate q
  input$par$logit_q[] = gen.logit(2,0,1000)  # fix it at 2
  input$par$logit_selpars[,1:9] = dummy$parList$logit_selpars[,1:9]  #start sel pars at decent place
  input
}

prep_rescale_absel<-function(input){
  dummy<-readRDS("Data/Dummy_inits/dummy.rds")
  input$map$logit_q = factor(rep(NA,input$data$n_indices)) # don't estimate q
  input$par$logit_q[] = gen.logit(2,0,1000)  # fix it at 2
  input$par$logit_selpars[,1:9] = dummy$parList$logit_selpars[c(1,5,3,4,5),1:9]  #start sel pars at decent place
  input
}

# The - random walk in M
M_RE_prep<-function(asap,mod="logistic-normal-miss0"){
  asap$comments<-paste0(asap$comments,"_M_RE_RESCALE")
  input <- prepare_wham_input(asap, recruit_model = 2, model_name=asap$comments, 
                              M=list(model="constant",re="ar1_y",est_ages=NULL),
                              age_comp = mod,
                              selectivity = make_Base_sel_init(asap))     
  prep_rescale(input)
}



M_RE_fixSD_prep<-function(asap,mod="logistic-normal-miss0", sigma = 0.1791195,rho_y=0.8354769, fix_rho_y = TRUE){
  
  asap$comments<-paste0(asap$comments,"_M_RE_fixSD_RESCALE")
  input <- prepare_wham_input(asap, recruit_model = 2, model_name=asap$comments, 
                              M=list(model="constant",re="ar1_y",est_ages=NULL),
                              age_comp = mod,
                              selectivity = make_Base_sel_init(asap)) 
  
  input$par$M_repars <- c(log(sigma),  0.000000,  invlogit(rho_y)) # inits for log(sigma), logit(rho_a), logit(rho_y)
  
  if(fix_rho_y) {
    input$map$M_repars <- as.factor(c(1,NA,NA)) # fix rho parameter
  }else{
    input$map$M_repars <- as.factor(c(NA,NA,1)) # fix sigma parameter
  }
  prep_rescale(input)
  
}


#  - random walk in M ar1y
M_RE_iid_prep<-function(asap,mod="logistic-normal-miss0"){
  asap$comments<-paste0(asap$comments,"_M_RE")
  
  selectivity = make_Base_sel_init(asap)
  selectivity$model[2]<-"age-specific"
  selectivity$initial_pars[[2]]<-c(0.1,0.3,0.5,0.7,0.9,1,1,1,1)
  selectivity$initial_pars[[4]][8]<-1
  selectivity$fix_pars[[3]]<-8
  selectivity$fix_pars[[4]]<-8
  selectivity$fix_pars[[5]]<-9
  
  prepare_wham_input(asap, recruit_model = 2, model_name=asap$comments, 
                              M=list(model="constant",re="ar1_y",est_ages=NULL),
                              age_comp = mod,
                              selectivity = selectivity)     
  
}


# Base model without any fixed selectivity parameters
M_RE_fix_prep<-function(asap,mod="logistic-normal-miss0"){
  asap$comments<-paste0(asap$comments,"_M_RE_RESCALE")
  input <- prepare_wham_input(asap, recruit_model = 2, model_name=asap$comments, 
                              M=list(model="constant",re="ar1_y",est_ages=NULL),
                              age_comp = mod,
                              selectivity = make_NMFS_sel_init(asap))     
  prep_rescale(input)
}

# An alternative approach in which M is estimated as a step up in recent ny years
M_est_prep<-function(asap,ny=10,mod="logistic-normal-miss0"){
  asap$comments<-paste0(asap$comments,"_M_est_RESCALE",ny)
  input <- prepare_wham_input(asap, recruit_model = 2, model_name=asap$comments,  
                              ecov=make_Mstep_ecov(asap=asap,ny=ny),
                              age_comp = mod,
                              selectivity = make_Base_sel_init(asap)) 
  prep_rescale(input)
}

M_est_prep_OM<-function(asap,ny=10,mod="logistic-normal-miss0"){
  asap$comments<-paste0(asap$comments,"_M_est",ny)
  input <- prepare_wham_input(asap, recruit_model = 2, model_name=asap$comments,  
                              ecov=make_Mstep_ecov(asap=asap,ny=ny),
                              age_comp = mod,
                              selectivity = make_Base_sel_init(asap)) 
  #prep_rescale(input)
}




# M is not estimated and a fixed value in the asap file
No_M_prep<-function(asap,mod="logistic-normal-miss0"){
  asap$comments<-paste0(asap$comments,"_Mfix_RESCALE")
  input <- prepare_wham_input(asap, recruit_model = 2, model_name=asap$comments,   
                              age_comp = mod,
                              selectivity = make_Base_sel_init(asap))     
  prep_rescale(input)
 
}

No_M_noselRE_prep<-function(asap,mod="logistic-normal-miss0"){
  asap$comments<-paste0(asap$comments,"_Mfix_RESCALE")
  input <- prepare_wham_input(asap, recruit_model = 2, model_name=asap$comments,   
                              age_comp = mod,
                              selectivity = make_noRE_sel_init(asap))     
  prep_rescale(input)
  
}



No_M_iid_prep<-function(asap,mod="logistic-normal-miss0"){
  asap$comments<-paste0(asap$comments,"_Mfix_seliid")
  
  selectivity = make_Base_sel_init(asap)
  selectivity$model[2]<-"age-specific"
  selectivity$re[2]<-"2dar1" # iid ar1 ar1_y 2dar1
  selectivity$initial_pars[[4]][8]<-1
  selectivity$initial_pars[[2]]<-c(0.1,0.3,0.5,0.7,0.9,1,1,1,1)
  
  selectivity$fix_pars[[3]]<-8
  selectivity$fix_pars[[4]]<-8
  selectivity$fix_pars[[5]]<-9  
  prepare_wham_input(asap, recruit_model = 2, model_name=asap$comments,   
                              age_comp = mod,
                              selectivity = selectivity) 
  
  
}  

# no fleet blocking - requires asap data file to have no blocking
Unblocked_prep<-function(asap,mod="logistic-normal-miss0"){
  asap$comments<-paste0(asap$comments,"_Unblocked_RESCALE")
  
  sel_list<-make_Base_sel_init(asap)
  ind<-2:5
  sel_list$re<-sel_list$re[ind]
  sel_list$initial_pars<-sel_list$initial_pars[ind]
  sel_list$fix_pars<-sel_list$fix_pars[ind]

  input<-prepare_wham_input(asap, recruit_model=2, model_name=asap$comments,
                     M=list(model="constant",re="ar1_y",est_ages=NULL),
                     age_comp = mod,
                     selectivity=sel_list)
  prep_rescale(input)
  
}

# no fall index - requires asap data file to have fall index turned off
Nofall_prep<-function(asap,mod="logistic-normal-miss0"){
  
  asap$comments<-paste0(asap$comments,"_NoFallInd_RESCALE")
  # Some init par vectors (maybe overwritten by inits below)
  sel_list<-make_Base_sel_init(asap)
  ind<-c(1,2,3,5)
  sel_list$re<-sel_list$re[ind]
  sel_list$initial_pars<-sel_list$initial_pars[ind]
  sel_list$fix_pars<-sel_list$fix_pars[ind]
  
  input<-prepare_wham_input(asap, recruit_model=2, model_name=asap$comments,
                     M=list(model="constant",re="ar1_y",est_ages=NULL),
                     age_comp = mod,
                     selectivity=sel_list)
  prep_rescale(input)
}


# Swaps M AR-1 yr for NAA by age and year
NAA_RE_prep<-function(asap,mod="logistic-normal-miss0"){
  
  asap$comments<-paste0(asap$comments,"_NAAre_RESCALE")
  input<-prepare_wham_input(asap, recruit_model=2, model_name=asap$comments,
                     NAA_re=list(sigma="rec+1",cor="iid"),
                     age_comp = mod,
                     selectivity=make_Base_sel_init(asap))
  prep_rescale(input)
  
}

# Selectivity is a 1D ar1-y process
Sel1D_prep<-function(asap,mod="logistic-normal-miss0"){
  
  asap$comments<-paste0(asap$comments,"_selAR1yr_RESCALE")
  sel_list<-make_Base_sel_init(asap)
  sel_list$re<-c("none", "ar1_y",  "none", "none", "none")
  input <- prepare_wham_input(asap, recruit_model = 2, model_name=asap$comments, 
                              M=list(model="constant",re="ar1_y",est_ages=NULL),
                              age_comp = mod,
                              selectivity = sel_list)     
  prep_rescale(input)
  
}

# Selectivity is a 2D ar1 process
Sel2D_prep<-function(asap,ny=10,mod="logistic-normal-miss0"){
  
  asap$comments<-paste0(asap$comments,"_selAR12D_RESCALE")
  sel_list<-make_Base_sel_init(asap)
  sel_list$re<-c("none", "2dar1",  "none", "none", "none")
  input <- prepare_wham_input(asap, recruit_model = 2, model_name=asap$comments, 
                              ecov=make_Mstep_ecov(asap=asap,ny=ny),
                              age_comp = mod,
                              selectivity = sel_list)     
  prep_rescale(input)
 
}



# Selectivity is a 2D ar1 process
Sel2D_noNAA_prep<-function(asap,mod="logistic-normal-miss0"){
  
  asap$comments<-paste0(asap$comments,"_selAR12D_RESCALE")
  sel_list<-make_Base_sel_init(asap)
  sel_list$re<-c("none", "2dar1",  "none", "none", "none")
  input <- prepare_wham_input(asap, recruit_model = 2, model_name=asap$comments, 
                              #M=list(model="constant",re="ar1_y",est_ages=NULL),
                              age_comp = mod,
                              selectivity = sel_list)     
  prep_rescale(input)
  
}

# Density dependent M - based on indices
DDM_prep<-function(asap,mod="logistic-normal-miss0"){
  
  asap$comments<-paste0(asap$comments,"_DDM_RESCALE")
  
  ecov <- list(
    label = "NMFS_FAll",
    mean = as.matrix(asap$dat$IAA_mats[[1]][,2]),
    logsigma = 'est_1', # estimate obs sigma, 1 value shared across years
    year = asap$dat$IAA_mats[[1]][,1],
    use_obs = matrix(1, ncol=1, nrow=asap$dat$n_years), # use all obs (=1)
    lag = 0, # NMFS_Fall in year t affects M in same year
    process_model = "ar1", # NMFS_Fall modeled as AR1 (random walk would be "rw")
    where = "M", # NMFS_Fall affects natural mortality
    how = 1, # include NMFS_Fall effect on M
    link_model = "linear")

  input<-prepare_wham_input(asap, recruit_model=2, model_name=asap$comments,
                     ecov=ecov,
                     age_comp = mod,
                     selectivity=make_Base_sel_init(asap))
  prep_rescale(input)
  
}

# RE on M but age not year
REMage_prep<-function(asap,mod="logistic-normal-miss0"){

  asap$comments<-paste0(asap$comments,"_M_RE_age_RESCALE")
  input <- prepare_wham_input(asap, recruit_model = 2, model_name=asap$comments, 
                              M=list(model="constant",re="ar1_a",est_ages=NULL),
                              age_comp = mod,
                              selectivity = make_Base_sel_init(asap))     
  prep_rescale(input)
  
}


# change fleet selectivity to age-specific in recent years  (can't estimate M then)
SAAlate_prep<-function(asap,ny=10,mod="logistic-normal-miss0"){
  
  asap$comments<-paste0(asap$comments,"_SAAlate_noMRE_RESCALE")
  sel_list<-make_Base_sel_init(asap)
  sel_list$model<-c("logistic","age-specific","age-specific","age-specific","age-specific")
  sel_list$re<-rep("none",5)
  sel_list$initial_pars[[2]]<-c(0.1,0.2,0.4, 0.6,0.8,0.9, 1,1,1)
  input<-prepare_wham_input(asap, recruit_model=2, model_name=asap$comments,
                            ecov=make_Mstep_ecov(asap=asap,ny=ny),
                            age_comp = mod,
                            selectivity=sel_list)
  prep_rescale(input)
  
}



SS_2dar1_prep<-function(asap,mod="logistic-normal-miss0"){
  
  asap$comments<-paste0(asap$comments,"_SS_prep_RESCALE")
  
  input<-prepare_wham_input(asap, recruit_model = 2, model_name = asap$comments,                         
                            NAA_re = list(cor='2dar1', sigma='rec+1'),
                            age_comp = mod,
                            selectivity = make_Base_sel_init(asap))  
  prep_rescale(input)
  
}

SS_iid_prep<-function(asap,mod="logistic-normal-miss0"){
  
  asap$comments<-paste0(asap$comments,"_SS_prep_RESCALE")
  
  input<-prepare_wham_input(asap, recruit_model = 2, model_name = asap$comments,                         
                            NAA_re = list(cor='iid', sigma='rec+1'),
                            age_comp = mod,
                            selectivity = make_Base_sel_init(asap))  
  prep_rescale(input)
  
}


# the f1e.free config of Liz Brooks designed for the GB

NMFS_prep<-function(asap,mod="logistic-normal-miss0"){
  
  asap$comments<-paste0(asap$comments,"_NMFS_GB_RESCALE")
  
  input<-prepare_wham_input(asap, recruit_model = 2, model_name = asap$comments,                         
                     NAA_re = list(cor='2dar1', sigma='rec+1'),
                     age_comp = mod,
                     selectivity = make_NMFS_sel_init(asap))  
  prep_rescale(input)
  
}


NMFS_MRE_prep<-function(asap,mod="logistic-normal-miss0"){
  
  asap$comments<-paste0(asap$comments,"_NMFS_MRE_RESCALE")
  
  input<-prepare_wham_input(asap, recruit_model = 2, model_name = asap$comments,                         
                            M=list(model="constant",re="ar1_a",est_ages=NULL),
                            age_comp = mod,
                            selectivity = make_NMFS_sel_init(asap))  
  prep_rescale(input)
  
}

NMFS_noRE_prep<-function(asap,mod="logistic-normal-miss0"){
  
  asap$comments<-paste0(asap$comments,"_NMFS_noRE_RESCALE")
  
  input<-prepare_wham_input(asap, recruit_model = 2, model_name = asap$comments,                         
                           age_comp = mod,
                            selectivity = make_NMFS_sel_init(asap))  
  prep_rescale(input)
  
}


