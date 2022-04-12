library(rjags)
library(MASS)
library(parallel)

modelstring = "
 model {
  # Likelihood of data given model parameters
  for (i in 1:N){
   class[i] ~ dbern(probdiff)
   Yobs[i] ~ dnorm(Y[i],tau_hat[i])
   Y[i] = m*Xobs[i] + c
   tau_hat[i] = ifelse(class[i]==0, tau, 0.01)
  }
  for (j in 1:Nsyn){
   Ysyn[j]~dnorm(Ys[j], tau)
   Ys[j] <- m*Xsyn[j] + c
  }
  # Specify prior beliefs about parameters
  m~dnorm(mu_m,tau_m)
  c~dnorm(mu_c,tau_c)
  tau~dgamma(shape_tau,rate_tau)
  probdiff~dbeta(alpha,beta)  
 }
"

dir.create(file.path("./Output"), showWarnings = FALSE)
dir.create(file.path("./Output/linReg_classif"), showWarnings = FALSE)

getData_mats = function(fulldat="Data_prepped.csv", mitochan="VDAC", chan, 
                        pts=NULL, ctrl_only=FALSE){
  data_raw = read.csv(fulldat, header=TRUE)
  
  pts_raw = unique(data_raw$patient_id)
  
  data = data_raw
  data[data[,"patient_id"] %in% pts_raw[grep("C0", pts_raw)],"patient_id"] = "control"
  pts_all = unique(data_raw$patient_id)
  
  ctrl_data = data[data$patient_id=="control", ]
  Xctrl = log(ctrl_data[[mitochan]])
  Yctrl = log(ctrl_data[[chan]])
  XY_ctrl = cbind( Xctrl, Yctrl )
  
  if(!ctrl_only){
    if(is.null(pts)){
      Ypat_all = matrix(NA, nrow=1, ncol=2)
      Npats = vector("numeric")
      for(pat in pts_all[grepl("P", pts_all)]){
        pat_data = data[data$patient_id==pat,]
        Xpat = log(pat_data[[mitochan]])
        Ypat = log(pat_data[[chan]])
        XY_pat = cbind(Xpat, Ypat)
        Ypat_all = rbind(Ypat_all, XY_pat)
        Npats[pat] = nrow(XY_pat)
      }
    } else {
      Ypat_all = matrix(NA, nrow=1, ncol=2)
      Npats = vector("numeric")
      for(pat in pts){
        pat_data = data[data$patient_id==pat,]
        Xpat = log(pat_data[[mitochan]])
        Ypat = log(pat_data[[chan]])
        XY_pat = cbind(Xpat, Ypat)
        Ypat_all = rbind(Ypat_all, XY_pat)
        Npats[pat] = nrow(XY_pat)
      }
    }
  }
  
  if(ctrl_only) return(XY_ctrl)
  return(list(ctrl=XY_ctrl, pat=Ypat_all[-1,], Npats=Npats))
}

colQuantiles = function(x, probs=0.5){
  quants = matrix(NA, nrow=ncol(x), ncol=length(probs))
  for(i in 1:ncol(x)){
    quants[i,] = quantile(x[,i], probs)
  }
  colnames(quants) = probs
  return(quants)
}

# inferenecce 
inference = function(input){
  with(c(input),{
    data_mats = getData_mats(chan=chan, pts=pat)
    Yctrl = data_mats$ctrl
    Ypat = data_mats$pat
    Nctrl = nrow(Yctrl)
    Npat = nrow(Ypat)
  
    # prior parameters for control data
    c_est = 0
    tau_c = 1/2^2
    m_est = 0
    tau_m = 1/2^2
    tau_err_shape = 3
    tau_err_rate = 1
    tau_err_mode = (tau_err_shape-1)/tau_err_rate
    N_syn = 1000
    
    ## control inference
    data_ctrl = list( Xobs=Yctrl[,1], Yobs=Yctrl[,2], N=Nctrl, Nsyn=N_syn,
                      Xsyn=seq(0.9*min(Yctrl[,1]), 1.1*max(Yctrl[,1]),length.out=N_syn),
                      mu_m=m_est, tau_m=tau_m, mu_c=c_est, tau_c=tau_c,
                      shape_tau=tau_err_shape, rate_tau=tau_err_rate,
                      alpha=1.0, beta=1.0)
    
    data_ctrl_priorpred = data_ctrl
    data_ctrl_priorpred$Yobs = NULL
    data_ctrl_priorpred$N = 0
    
    model_ctrl = jags.model(textConnection(modelstring), data=data_ctrl, n.chains=n.chains)
    model_ctrl_priorpred = jags.model(textConnection(modelstring), data=data_ctrl_priorpred)
    
    update(model_ctrl,n.iter=MCMCBurnin)
    output_ctrl = coda.samples(model=model_ctrl, n.iter=MCMCOut*MCMCThin,thin=MCMCThin,
                               variable.names=c("m","c","tau","Ysyn","class","probdiff"))
    output_ctrl_priorpred = coda.samples(model=model_ctrl_priorpred, n.iter=MCMCOut,thin=1,
                                         variable.names=c("m","c","tau","Ysyn","probdiff"))
    
    posterior_ctrl = as.data.frame(output_ctrl[[1]])
    prior_ctrl = as.data.frame(output_ctrl_priorpred[[1]])

    summ_ctrl = summary(output_ctrl)
    classifs_ctrl = summ_ctrl$statistics[grepl("class",rownames(summ_ctrl$statistics)),"Mean"]
    
    ### patient inference
    # pateint priors
    flex = 0.01
    c_est = mean(posterior_ctrl$c)
    tau_c = flex/(sd(posterior_ctrl$c)^2)
    m_est = mean(posterior_ctrl$m)
    tau_m = flex/(sd(posterior_ctrl$m)^2)
    delta = 1.5*as.numeric(quantile(posterior_ctrl$tau,0.5)) # Choose this value so that Tiago's replication dataset never predicts over-expression of CI or CIV
    tau_err_mean = mean(posterior_ctrl$tau) + delta # Precision tau = (1/sd)^2
    tau_err_sd = sd(posterior_ctrl$tau) # Deviation from prior tau should require a lot of contradictory data
    tau_err_shape = (tau_err_mean^2)/(tau_err_sd^2)
    tau_err_rate = tau_err_mean/(tau_err_sd^2)
    N_syn = 1000
    
    data_pat = list( Xobs=Ypat[,1], Yobs=Ypat[,2], N=Npat, Nsyn=N_syn,
                     Xsyn=seq(0.9*min(Ypat[,1]), 1.1*max(Ypat[,1]), length.out=N_syn),
                     mu_m=m_est, tau_m=tau_m, 
                     mu_c=c_est, tau_c=tau_c,
                     shape_tau=tau_err_shape, rate_tau=tau_err_rate, 
                     alpha=12.5, beta=12.5)
    
    data_pat_priorpred = data_pat
    data_pat_priorpred$Yobs = NULL
    data_pat_priorpred$N = 0
    
    model_pat = jags.model(textConnection(modelstring), data=data_pat, n.chains=n.chains)
    model_pat_priorpred = jags.model(textConnection(modelstring), data=data_pat_priorpred)
    update(model_pat, n.iter=MCMCBurnin)
    # converge_pat = coda.samples(model=model_pat,variable.names=c("m","c","tau_par","class","probdiff"),n.iter=MCMCUpdates_Report,thin=MCMCUpdates_Thin)
    output_pat_post = coda.samples(model=model_pat, n.iter=MCMCOut*MCMCThin, thin=MCMCThin,
                              variable.names=c("m","c","tau","Ysyn","class","probdiff"))
    output_pat_prior = coda.samples(model=model_pat_priorpred,n.iter=MCMCOut, thin=1,
                                        variable.names=c("m","c","tau","Ysyn","probdiff"))
    
    posterior_pat = as.data.frame(output_pat_post[[1]])
    prior_pat = as.data.frame(output_pat_prior[[1]])
    
    summ_pat = summary(output_pat_post)
    classifs_pat = summ_pat$statistics[grepl("class",rownames(summ_pat$statistics)),"Mean"]
    
    posterior_ctrl_names = colnames(posterior_ctrl)
    post_ctrl = posterior_ctrl[,!(grepl("class", posterior_ctrl_names)|grepl("Ysyn", posterior_ctrl_names))]
    
    postpred_ctrl = colQuantiles( posterior_ctrl[,grepl("Ysyn", posterior_ctrl_names)], probs=c(0.025, 0.5, 0.975) )
    
    prior_ctrl_names = colnames(prior_ctrl)
    priorpred_ctrl = colQuantiles(prior_ctrl[, grepl("Ysyn", prior_ctrl_names)], probs=c(0.025,0.5,0.975))
    prior_control = prior_ctrl[,!grepl("Ysyn", prior_ctrl_names)]
    
    posterior_pat_names = colnames(posterior_pat)
    post_pat = posterior_pat[,!(grepl("class", posterior_pat_names)|grepl("Ysyn",posterior_pat_names))]
    postpred_pat = colQuantiles(posterior_pat[,grepl("Ysyn", posterior_pat_names)], probs=c(0.025,0.5,0.975))
    
    prior_pat_names = colnames(prior_pat)
    priorpred_pat = colQuantiles(prior_pat[,grepl("Ysyn", prior_pat_names)], probs=c(0.025,0.5,0.975))
    prior_patient = prior_pat[,!grepl("Ysyn",prior_pat_names)]
    
    ctrl_inference = list(post=post_ctrl, postpred=postpred_ctrl, 
                          prior=prior_control, priorpred=priorpred_ctrl,
                          classif=classifs_ctrl)
    pat_inference = list(post=posterior_pat, postpred=postpred_pat,
                         prior=prior_patient, priorpred=priorpred_pat,
                         classif=classifs_pat)
    
    return( list(ctrl=ctrl_inference, pat=pat_inference) )
  })
}

fp_data = "Data_prepped.csv"
mc_raw = read.csv(fp_data, header=TRUE)
pts_raw = unique(mc_raw$patient_id)
mc_data = mc_raw
mc_data[mc_raw[,"patient_id"] %in% pts_raw[grep("C0", pts_raw)],"patient_id"] = "control"
pts_all = unique(mc_data$patient_id)
pts  = pts_all[pts_all!="control"]
mitochan = "VDAC"
channels = c("MTCO1", "NDUFB8", "CYB")

inputs = list()
{
  input0 = list()
  input0$MCMCOut = 10000
  input0$MCMCBurnin = 2000
  input0$MCMCThin = 1
  input0$n.chains = 1
  i = 0
  for(chan in channels){
    for(pat in pts){
      i = i + 1
      inputs[[i]] = input0
      inputs[[i]]$chan = chan
      inputs[[i]]$pat = pat
    } # pts
  } # chans
}

ncores = 12
cl  = makeCluster(ncores) 
{
  clusterExport(cl, c("getData_mats", "modelstring", "colQuantiles"))
  clusterEvalQ(cl, {
    library("rjags")
  })
  linreg_output = parLapply(cl, inputs, inference)
}
stopCluster(cl)

# save output
for(index in 1:length(linreg_output)){
  output = linreg_output[[index]]
  chan = inputs[[index]]$chan
  pat = gsub("_",".",inputs[[index]]$pat)
  
  outroot_pat = paste(chan, "CONTROL", sep="_")
  outroot_ctrl = paste(chan, pat, sep="_")
  
  output = linreg_output[[index]]
  
  write.table(output[["pat"]]$post, paste0("./Output/linReg_classif/", outroot_pat,"_POST.txt"), 
              row.names=FALSE, col.names=TRUE)
  write.table(output[["pat"]]$postpred, paste0("./Output/linReg_classif/", outroot_pat,"_POSTPRED.txt"), 
              row.names=FALSE, col.names=TRUE)
  write.table(output[["pat"]]$classif, paste0("./Output/linReg_classif/",outroot_pat,"_CLASS.txt"), 
              row.names=FALSE, col.names=TRUE)
  write.table(output[["pat"]]$prior, paste0("./Output/linReg_classif/", outroot_pat,"_PRIOR.txt"), 
              row.names=FALSE, col.names=TRUE)
  write.table(output[["pat"]]$priorpred, paste0("./Output/linReg_classif/", outroot_pat,"_PRIORPRED.txt"), 
              row.names=FALSE, col.names=TRUE)
  
  write.table(output[["ctrl"]]$post, paste0("./Output/linReg_classif/", outroot_ctrl,"_POST.txt"), 
              row.names=FALSE, col.names=TRUE)
  write.table(output[["ctrl"]]$postpred, paste0("./Output/linReg_classif/", outroot_ctrl,"_POSTPRED.txt"), 
              row.names=FALSE, col.names=TRUE)
  write.table(output[["ctrl"]]$classif, paste0("./Output/linReg_classif/", outroot_ctrl,"_CLASS.txt"), 
              row.names=FALSE, col.names=TRUE)
  write.table(output[["ctrl"]]$prior, paste0("./Output/linReg_classif/", outroot_ctrl,"_PRIOR.txt"), 
              row.names=FALSE, col.names=TRUE)
  write.table(output[["ctrl"]]$priorpred, paste0("./Output/linReg_classif/", outroot_ctrl,"_PRIORPRED.txt"), 
              row.names=FALSE, col.names=TRUE)
}






