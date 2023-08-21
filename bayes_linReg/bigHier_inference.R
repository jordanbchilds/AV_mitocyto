bigHier_inference = function (dataMats,
          parameterVals = NULL,
          MCMCout = 1000,
          MCMCburnin = 1000,
          MCMCthin = 1)
{
  modelstring = "
  model {
    for(i in 1:nCtrl){
      Yctrl[i] ~ dnorm(m[indexCtrl[i]]*Xctrl[i]+c[indexCtrl[i]], tau_norm)
    }
    for(j in 1:nPat){
      Ypat[j] ~ dnorm(m[indexPat[j]]*Xpat[j]+c[indexPat[j]], tau_hat[j])
      tau_hat[j] = ifelse(class[j]==1, tau_def, tau_norm)
      class[j] ~ dbern(probdiff[indexPat[j] - nCrl])
    }
    
    for( k in 1:nPts){
        for(l in 1:nSyn){
          Ysyn_norm[l, k] ~ dnorm(m[k + nCrl]*Xsyn[l]+c[k + nCrl], tau_norm)
          Ysyn_def[l, k] ~ dnorm(m[k + nCrl]*Xsyn[l]+c[k + nCrl], tau_def)
        }
    }
    
    for(ii in 1:(nCrl+nPts)){
      m[ii] ~ dnorm(mu_m, tau_m)
      c[ii] ~ dnorm(mu_c, tau_c)
    }
    
    for( jj in 1:nPts ){
        probdiff[jj] ~ dbeta(alpha_pi, beta_pi)
    }
    
    m_pred ~ dnorm(mu_m, tau_m)
    c_pred ~ dnorm(mu_c, tau_m)
    
    mu_m ~ dnorm( mean_mu_m, prec_mu_m )
    mu_c ~ dnorm( mean_mu_c, prec_mu_c )
    tau_m ~ dgamma( shape_tau_m, rate_tau_m )
    tau_c ~ dgamma( shape_tau_c, rate_tau_c )
    tau_norm ~ dgamma(shape_tau, rate_tau)

  }
  "
  if (is.null(names(dataMats)))
    names(dataMats) = c("ctrl", "pts", "indexCtrl", "indexPat")
  
  ctrl_mat = dataMats$ctrl
  pat_mat = dataMats$pts
  nCtrl = nrow(ctrl_mat)
  nPat = nrow(pat_mat)
  indexCtrl = dataMats$indexCtrl
  indexPat = dataMats$indexPat
  nCrl = length(unique(indexCtrl))
  nPts = length(unique(indexPat))
  nSyn = 1000
  Xsyn = seq(min(0, min(c(ctrl_mat[, 1], pat_mat[, 1]))), 
             max(c(ctrl_mat[,1], pat_mat[, 1])) * 1.5, 
             length.out = nSyn)
  data_list = list(
    Yctrl = ctrl_mat[, 2],
    Xctrl = ctrl_mat[, 1],
    Ypat = pat_mat[, 2],
    Xpat = pat_mat[, 1],
    nCtrl = nCtrl,
    nPts = nPts, 
    nCrl = nCrl,
    nPat = nPat,
    indexCtrl = indexCtrl,
    indexPat = indexPat,
    nSyn = nSyn,
    Xsyn = Xsyn,
    mean_mu_m = 1,
    prec_mu_m = 1 / 0.1 ^ 2,
    mean_mu_c = 0,
    prec_mu_c = 1 / 0.2 ^ 2,
    shape_tau_m = 27.56372,
    rate_tau_m = 1.660233,
    shape_tau_c = 1.25,
    rate_tau_c = 0.25,
    shape_tau = 41.97618,
    rate_tau = 2.048809,
    alpha_pi = 1,
    beta_pi = 1,
    tau_def = 1e-04
  )
  if (!is.null(parameterVals) && is.list(parameterVals)) {
    for (param in names(parameterVals)) {
      if (param %in% names(data_list)) {
        data_list[[param]] = parameterVals[[param]]
      }
      else {
        message(paste("The parameter `", param, "` is not part of the model."))
      }
    }
  }
  data_priorpred = data_list
  data_priorpred$Yctrl = NULL
  data_priorpred$nCtrl = 0
  data_priorpred$Ypat = NULL
  data_priorpred$nPat = 0
  model_pat = rjags::jags.model(textConnection(modelstring),
                                data = data_list,
                                n.chains = 1)
  model_pat_priorpred = rjags::jags.model(textConnection(modelstring),
                                          data = data_priorpred)
  output_post = rjags::coda.samples(
    model = model_pat,
    n.iter = MCMCout *
      MCMCthin,
    thin = MCMCthin,
    variable.names = c(
      "m",
      "c",
      "m_pred",
      "c_pred",
      "mu_m",
      "tau_m",
      "mu_c",
      "tau_c",
      "tau_norm",
      "Ysyn_norm",
      "Ysyn_def",
      "class",
      "probdiff"
    )
  )
  output_prior = rjags::coda.samples(
    model = model_pat_priorpred,
    n.iter = MCMCout,
    thin = 1,
    variable.names = c(
      "m",
      "c",
      "m_pred",
      "c_pred",
      "mu_m",
      "tau_m",
      "mu_c",
      "tau_c",
      "tau_norm",
      "Ysyn_norm",
      "Ysyn_def",
      "probdiff"
    )
  )
  post_out = as.data.frame(output_post[[1]])
  prior_out = as.data.frame(output_prior[[1]])
  summ_pat = summary(output_post)
  classifs_mat = post_out[, grepl("class", colnames(post_out))]
  
  post_names = colnames(post_out)
  post = post_out[, c(
    paste0("m[", 1:(nCrl + nPts), "]"),
    paste0("c[",1:(nCrl + nPts), "]"),
    "m_pred",
    "c_pred",
    "mu_m",
    "tau_m",
    "mu_c",
    "tau_c",
    "tau_norm",
    paste0("probdiff[", 1:nPts, "]" )
  )]
  
  prior_names = colnames(prior_out)
  prior = prior_out[, c(
    paste0("m[", 1:(nCrl + nPts), "]"),
    paste0("c[",1:(nCrl + nPts), "]"),
    "m_pred",
    "c_pred",
    "mu_m",
    "tau_m",
    "mu_c",
    "tau_c",
    "tau_norm",
    paste0("probdiff[", 1:nPts, "]" )
  )]
  
  postpred_norm = apply( post_out[, grep("Ysyn_norm", colnames(post_out))], 2, quantile, probs=c(0.025, 0.5, 0.975))
  postpred_def = apply( post_out[, grep("Ysyn_def", colnames(post_out))], 2, quantile, probs=c(0.025, 0.5, 0.975))
  priorpred_norm = apply( prior_out[, grep("Ysyn_norm", colnames(prior_out))], 2, quantile, probs=c(0.025, 0.5, 0.975))
  priorpred_def = apply( prior_out[, grep("Ysyn_def", colnames(prior_out))], 2, quantile, probs=c(0.025, 0.5, 0.975))
  
  postpred = list()
  priorpred = list()
  for( i in 1:nPts ){
    postpred[[i]] = cbind(Xsyn, t(postpred_norm[, grep(paste0(",",i,"]"), colnames(postpred_norm))]), t(postpred_def[, grep(paste0(",",i,"]"), colnames(postpred_norm))]))
    colnames(postpred[[i]]) = c("mitochan", "lwr_norm", "med_norm", "upr_norm", "lwr_def", "med_def", "upr_def")
    priorpred[[i]] =  cbind(Xsyn, t(priorpred_norm[, grep(paste0(",",i,"]"), colnames(priorpred_norm))]), t(priorpred_def[, grep(paste0(",",i,"]"), colnames(priorpred_norm))]))
    colnames(priorpred[[i]]) = c("mitochan", "lwr_norm", "med_norm", "upr_norm", "lwr_def", "med_def", "upr_def")
  }
  
  classifs = list()
  for( i in seq_along(unique(indexPat))) {
    classifs[[i]] = classifs_mat[ ,indexPat==unique(indexPat)[i] ]
  }
  
  out_list = list(
    POST = post,
    POSTPRED = postpred,
    PRIOR = prior,
    PRIORPRED = priorpred,
    CLASSIF = classifs
  )
  return(out_list)
}
