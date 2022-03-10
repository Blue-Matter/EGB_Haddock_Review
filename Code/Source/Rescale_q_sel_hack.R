
# rescaling hack by Tim Miller and Liz Brooks to allow selectivity to be estimated unbounded with a fixed q, then rescaled by max selectivity- thus preventing non-convergence at selectivity=1

gen.logit <- function(x, low, upp) return(log((x-low)/(upp-x)))  # (tim)

inv.logit <- function(p, low, upp) return( (exp(p)*upp + low)/(1+exp(p))  ) # (liz)

rescale.out <- function(mod) {    # (liz)
  
  n.ind <- mod$input$data$n_indices
  n.yrs <- mod$input$data$n_years_model
  n.ages <- mod$input$data$n_ages
  ind.blks <- c( (mod$input$data$n_selblocks-n.ind+1) : mod$input$data$n_selblocks)
  max_sel <- unlist(lapply(mod$rep$selAA[ind.blks], max))
  mod$orig.selAA <- mod$rep$selAA
  mod$orig.q <- mod$rep$q
  mod$rep$q_yr_rescale <- matrix(NA, nrow=n.yrs, ncol=n.ind )
  for (i in 1:n.ind) {
    
    #get max selectivity at age across years, use to scale original (fixed) q (per index)
    tmp.sel <- matrix(unlist(mod$rep$selAA[ind.blks[i]] ), nrow=n.yrs, ncol=n.ages, byrow=F )
    max.sel <- apply(tmp.sel, 1, max)
    mod$rep$q_yr_rescale[,i] <- mod$orig.q[i]*max.sel
    
    #scale selectivity at age to have max of 1 in each year (per index)
    tmp <- as.matrix(as.data.frame(mod$rep$selAA[ind.blks[i]] ))
    tmp2 <- t(apply(tmp, 1, function (x) x/max(x)   ))
    mod$rep$selAA[ind.blks[i]] <- list(tmp2)
  } # end i-loop over n.ind
  return(mod)
  
}

rescale.samp <- function(mod,output) {    # (liz)
  
  n.ind <- mod$input$data$n_indices
  n.yrs <- mod$input$data$n_years_model
  n.ages <- mod$input$data$n_ages
  ind.blks <- c( (mod$input$data$n_selblocks-n.ind+1) : mod$input$data$n_selblocks)
  max_sel <- unlist(lapply(mod$rep$selAA[ind.blks], max))
  mod$orig.selAA <- mod$rep$selAA
  mod$orig.q <- mod$rep$q
  mod$rep$q_yr_rescale <- matrix(NA, nrow=n.yrs, ncol=n.ind )
  nsim<-length(output)
  
  for(ss in 1:nsim){
  for (i in 1:n.ind) {
    
    #get max selectivity at age across years, use to scale original (fixed) q (per index)
    tmp.sel <- matrix(unlist(output[[ss]]$selAA[ind.blks[i]] ), nrow=n.yrs, ncol=n.ages, byrow=F )
     
    #scale selectivity at age to have max of 1 in each year (per index)
    tmp <- as.matrix(as.data.frame(tmp.sel))
    tmp2 <- t(apply(tmp, 1, function (x) x/max(x)   ))
    output[[ss]]$selAA[ind.blks[i]] <- list(tmp2)
  } # end i-loop over n.ind
  }  
  
  return(output)
  
}






conv_info <- function(mod,rseed=5302991)  {   # (liz)
  
  check <- capture.output(check_convergence(mod) )
  ok_sdrep = (if(mod$na_sdrep==FALSE & !is.na(mod$na_sdrep)) 1 else 0)
  conv <- ifelse(mod$opt$convergence == 0, T, F)
  pdHess <- as.logical(ok_sdrep)
  check_convergence <- list(check=check, convergence=conv, max.grad=max(abs(mod$final_gradient)), 
                            max.grad.par=names(mod$par)[which.max(mod$final_gradient)], 
                            max.grad.par.number=which(mod$final_gradient==max(mod$final_gradient)) ,
                            pdHess=pdHess , rseed=rseed,  run.time=9999   )
  
  return(check_convergence)
}  # end function conv_info
# end functions  ====