# inferenecce   
stan_inference = function(dataMats, MCMCout=1000, 
                          MCMCburnin=1000, MCMCthin=1, nChains=10, 
                          max_logLik=TRUE, nCores=1){
  
  nCtrl = nrow( dataMats$ctrl )
  nPat = nrow( dataMats$pts )
  nCrl = length( unique( dataMats$indexCtrl) )
  parameterVals = dataMats$parameterVals
  
  nSyn = 1e3
  xSyn = seq(min(c(dataMats$ctrl[,1], dataMats$pts[,1]))-2, max(c(dataMats$ctrl[,1], dataMats$pts[,1]))+2, length.out=nSyn)
  param_list = list(
    D=2, 
    K=2,
    M=nCrl+1,
    ctrl_mat = dataMats$ctrl,
    pat_mat = dataMats$pts,
    nCtrl = nCtrl,
    nCrl = nCrl,
    nPat = nPat,
    nPts=1,
    ctrlIndex = dataMats$indexCtrl,
    nSyn = nSyn,
    xSyn = xSyn,
    mean_mu_m = 1.0,
    prec_mu_m = 1/0.1^2, 
    mean_mu_c = 0.0,
    prec_mu_c = 1/0.2^2, 
    shape_tau_m = 27.56372, 
    rate_tau_m = 1.660233, 
    shape_tau_c = 1.25, 
    rate_tau_c = 0.25,
    shape_tau = 41.97618,
    rate_tau = 2.048809,
    alpha_pi = 0.0,
    beta_pi = 1.0,
    slope_lb = 0.1,
    tau_def=1e-4
  )
  
  if (!is.null(parameterVals) && is.list(parameterVals)) {
    for (param in names(parameterVals) ) {
      if (param %in% names(param_list) ) {
        param_list[[param]] = parameterVals[[param]]
      } else {
        message(paste("The parameter `", param, "` is not part of the model."))
      }
    }
  }
  
  init_list = list()
  for(i in 1:nChains){
    gen_list = list("mu_m"=rnorm(1, param_list[["mean_mu_m"]], 1/sqrt(param_list[["prec_mu_m"]]) ), 
                    "mu_c"=rnorm(1, param_list[["mean_mu_c"]], 1/sqrt(param_list[["prec_mu_c"]])),
                    "tau_m"= rgamma(1, param_list[["shape_tau_m"]], param_list[["rate_tau_m"]]),
                    "tau_c"=rgamma(1, param_list[["shape_tau_c"]], param_list[["rate_tau_c"]]),
                    "tau_norm"=rgamma(1, param_list[["shape_tau"]], param_list[["rate_tau"]]),
                    "probdiff"=runif(1,param_list[["alpha_pi"]],  param_list[["beta_pi"]] ) )
    m_init = list()
    c_init = list()
    for(j in 1:(nCrl+1)){
      m_init[[paste0("m[",j,"]")]] = rnorm(1, gen_list[["mu_m"]], 1/gen_list[["tau_m"]])
      c_init[[paste0("c[",j,"]")]] = rnorm(1, gen_list[["mu_c"]], 1/gen_list[["tau_c"]])
    }
    init_list[[i]] = c(gen_list, m_init, c_init)
  }
  
  output = stan("./bhlmm_stan.stan", data=param_list, chains=nChains, 
                init=init_list, cores=nCores,
                iter=(MCMCout+MCMCburnin)*MCMCthin, warmup=MCMCburnin, thin=MCMCthin)

  
  outmat = as.matrix(output)
  outcols = colnames(outmat)
  
  if( max_logLik ){
    logLik_mat = matrix(outmat[,"log_lik"], nrow=MCMCout, ncol=nChains, byrow=FALSE)
    ll_mean = colMeans(logLik_mat)
    maxChain = which.max(ll_mean)
    chain_ind = ((maxChain - 1)*MCMCout + 1):(maxChain*MCMCout)
    outmat = outmat[chain_ind, ]
    nChains = 1
  }
  
  classifs_mat = outmat[, grepl("classif", outcols)]
  post = outmat[, !( grepl("classif", outcols)|grepl("lp__", outcols)|
                       grepl("probvec", outcols)|grepl("dens", outcols)|
                       grepl("yPred", outcols)|grepl("_tmp", outcols)|
                       grepl("_prior", outcols) ) ]
  postpred_mat = outmat[, grepl("yPred", outcols) & !grepl("_prior", outcols)]
  
  prior = outmat[, !( grepl("classif", outcols)|grepl("lp__", outcols)|
                        grepl("probvec", outcols)|grepl("dens", outcols)|
                        grepl("yPred", outcols)|grepl("_tmp", outcols) ) & grepl("_prior", outcols) ]
  colnames(prior) = gsub("_prior", "", colnames(prior))
  
  priorpred_mat = outmat[, grepl("yPred", outcols) & grepl("_prior", outcols)]
  
  if( nChains == 1 ){
    postpred = apply(postpred_mat, 2, quantile, probs=c(0.025, 0.5, 0.975))
    postpred = cbind(xSyn, t(postpred))
    colnames(postpred) = c("mitochan", "lwrNorm", "medNorm", "uprNorm")
    
    priorpred = apply(priorpred_mat, 2, quantile, probs=c(0.025, 0.5, 0.975))
    priorpred = cbind(xSyn, t(priorpred))
    colnames(priorpred) = c("mitochan", "lwrNorm", "medNorm", "uprNorm")
  } else {
    postpred = NULL
    priorpred = NULL
    
    for( i in 1:nChains ){
      chain_ind = ((i - 1)*MCMCout + 1):(i*MCMCout)
      
      postpred_chain = apply(postpred_mat[chain_ind, ], 2, quantile, probs=c(0.025, 0.5, 0.975))
      priorpred_chain = apply(priorpred_mat[chain_ind, ], 2, quantile, probs=c(0.025, 0.5, 0.975))
      postpred = rbind(postpred, t(postpred_chain))
      priorpred = rbind(priorpred, t(priorpred_chain))
    }
    postpred = cbind(rep(xSyn, nChains), postpred)
    priorpred = cbind(rep(xSyn, nChains), priorpred)
    
    colnames(postpred) = c("mitochan", "lwrNorm", "medNorm", "uprNorm")
    colnames(priorpred) = c("mitochan", "lwrNorm", "medNorm", "uprNorm")
  }
  # colnames(priorpred) = gsub("_prior", "", colnames(priorpred))
  
  return( list(POST=post, POSTPRED=postpred, CLASSIF=classifs_mat, 
               PRIOR=prior, PRIORPRED=priorpred ) )
}
