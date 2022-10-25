library(rjags)
library(MASS)
library(parallel)
source("helper_functions.R", local=TRUE)

folder = "linReg_classifier"

dir.create(file.path("Output"), showWarnings = FALSE)
dir.create(file.path("Output", folder), showWarnings = FALSE)

modelstring = "
 model {
  # Likelihood of data given model parameters
  for (i in 1:N){
   class[i] ~ dbern(probdiff)
   Yobs[i] ~ dnorm(Y[i], tau_hat[i])
   Y[i] = m*Xobs[i] + c
   tau_hat[i] = ifelse(class[i]==1, tau_ctrl, tau_notctrl)
  }
  for(j in 1:Nsyn){
   Ysyn[j] ~ dnorm(Ys[j], tau_ctrl)
   Ys[j] <- m*Xsyn[j] + c
  }
  # Specify prior beliefs about parameters
  m ~ dnorm(mu_m, tau_m)
  c ~ dnorm(mu_c, tau_c)
  tau_ctrl ~ dgamma(shape_tau, rate_tau)
  notctrl_var ~ dgamma(shape_gamma, rate_gamma)
  tau_notctrl = 1/notctrl_var
  p ~ dbeta(alpha, beta)
  probdiff = ifelse(ctrl_ind==1, 1, p)
 }
"
# inferenecce 
inference = function(input){
  with(c(input),{
    data_mats = getData_mats(chan=chan, pts=pat)
    ctrl_mat = data_mats$ctrl
    pat_mat = data_mats$pts
    Nctrl = nrow(ctrl_mat)
    Npat = nrow(pat_mat)
    
    # prior parameters for control data
    c_est = 0
    tau_c = 1/2^2
    m_est = 0
    tau_m = 1/2^2
    shape_tau = 10
    rate_tau = 1
    shape_gamma = 20 
    rate_gamma = 10 
    N_syn = 1e4
    
    ## control inference
    data_ctrl = list( Xobs=ctrl_mat[,1], Yobs=ctrl_mat[,2], N=Nctrl, Nsyn=N_syn,
                      Xsyn=seq(0, max(c(ctrl_mat[,1], pat_mat[,1]))*1.5, length.out=N_syn),
                      mu_m=m_est, tau_m=tau_m, mu_c=c_est, tau_c=tau_c,
                      shape_tau=shape_tau, rate_tau=rate_tau,
                      shape_gamma=shape_gamma, rate_gamma=rate_gamma,
                      alpha=1.0, beta=1.0,
                      ctrl_ind=1)
    
    data_ctrl_priorpred = data_ctrl
    data_ctrl_priorpred$Yobs = NULL
    data_ctrl_priorpred$N = 0
    
    model_ctrl = jags.model(textConnection(modelstring), data=data_ctrl, n.chains=n.chains)
    model_ctrl_priorpred = jags.model(textConnection(modelstring), data=data_ctrl_priorpred)
    
    update(model_ctrl,n.iter=MCMCBurnin)
    output_ctrl = coda.samples(model=model_ctrl, n.iter=MCMCOut*MCMCThin,thin=MCMCThin,
                               variable.names=c("m","c","tau_ctrl","tau_notctrl","Ysyn","class","probdiff"))
    output_ctrl_prior = coda.samples(model=model_ctrl_priorpred, n.iter=MCMCOut,thin=1,
                                     variable.names=c("m","c","tau_ctrl","tau_notctrl","Ysyn","probdiff"))
    
    posterior_ctrl = as.data.frame(output_ctrl[[1]])
    prior_ctrl = as.data.frame(output_ctrl_prior[[1]])
    summ_ctrl = summary(output_ctrl)
    classifs_ctrl = summ_ctrl$statistics[grepl("class",rownames(summ_ctrl$statistics)),"Mean"]
    
    ### patient inference
    # pateint priors
    flex = 0.01
    c_est = mean(posterior_ctrl$c)
    tau_c = flex/(sd(posterior_ctrl$c)^2)
    m_est = mean(posterior_ctrl$m)
    tau_m = flex/(sd(posterior_ctrl$m)^2)
    # delta = 1.5*as.numeric(quantile(posterior_ctrl$tau_ctrl,0.5)) # Choose this value so that Tiago's replication dataset never predicts over-expression of CI or CIV
    delta = 0
    tau_mean = mean(posterior_ctrl$tau_ctrl) + delta # Precision tau = (1/sd)^2
    tau_sd = sd(posterior_ctrl$tau_ctrl) # Deviation from prior tau should require a lot of contradictory data
    tau_shape = (tau_mean^2)/(tau_sd^2)
    tau_rate = tau_mean/(tau_sd^2)
    
    gamma_mode = 20 
    gamma_sd = 10 
    rate_gamma = (gamma_mode+sqrt(gamma_mode^2+4*gamma_mode^2))/(2*gamma_sd^2)
    shape_gamma = 1+gamma_mode*rate_gamma

    data_pat = list(Xobs=pat_mat[,1], Yobs=pat_mat[,2], N=Npat, Nsyn=N_syn,
                    Xsyn=seq(0, max(c(ctrl_mat[,1], pat_mat[,1]))*1.5, length.out=N_syn),
                    mu_m=m_est, tau_m=tau_m, 
                    mu_c=c_est, tau_c=tau_c,
                    shape_tau=tau_shape, rate_tau=tau_rate, 
                    shape_gamma=shape_gamma, rate_gamma=rate_gamma, 
                    alpha=1, beta=1,
                    ctrl_ind=0)
    
    data_pat_priorpred = data_pat
    data_pat_priorpred$Yobs = NULL
    data_pat_priorpred$N = 0
    
    model_pat = jags.model(textConnection(modelstring), data=data_pat, n.chains=n.chains)
    model_pat_priorpred = jags.model(textConnection(modelstring), data=data_pat_priorpred)
    update(model_pat, n.iter=MCMCBurnin)
    #converge_pat=coda.samples(model=model_pat,variable.names=c("m","c","tau_par","class","probdiff"),n.iter=MCMCUpdates_Report,thin=MCMCUpdates_Thin)
    output_pat_post = coda.samples(model=model_pat, n.iter=MCMCOut*MCMCThin, thin=MCMCThin,
                            variable.names=c("m","c","tau_ctrl","tau_notctrl","Ysyn","class","probdiff"))
    output_pat_prior = coda.samples(model=model_pat_priorpred,n.iter=MCMCOut, thin=1,
                                      variable.names=c("m","c","tau_ctrl","tau_notctrl","Ysyn","probdiff"))
    
    posterior_pat = as.data.frame(output_pat_post[[1]])
    prior_pat = as.data.frame(output_pat_prior[[1]])
    
    summ_pat = summary(output_pat_post)
    classifs_pat = summ_pat$statistics[grepl("class",rownames(summ_pat$statistics)),"Mean"]
    
    posterior_ctrl_names = colnames(posterior_ctrl)
    post_ctrl = posterior_ctrl[,c("m", "c", "tau_ctrl", "tau_notctrl", "probdiff")]
    postpred_ctrl = colQuantiles( posterior_ctrl[,grepl("Ysyn", posterior_ctrl_names)], probs=c(0.025, 0.5, 0.975) )
    
    prior_ctrl_names = colnames(prior_ctrl)
    priorpred_ctrl = colQuantiles(prior_ctrl[, grepl("Ysyn", prior_ctrl_names)], probs=c(0.025,0.5,0.975))
    prior_control = prior_ctrl[,c("m", "c", "tau_ctrl", "tau_notctrl", "probdiff")]
    
    posterior_pat_names = colnames(posterior_pat)
    post_pat = posterior_pat[,c("m", "c", "tau_ctrl", "tau_notctrl", "probdiff")]
    postpred_pat = colQuantiles(posterior_pat[,grepl("Ysyn", posterior_pat_names)], probs=c(0.025,0.5,0.975))
    
    prior_pat_names = colnames(prior_pat)
    priorpred_pat = colQuantiles(prior_pat[,grepl("Ysyn", prior_pat_names)], probs=c(0.025,0.5,0.975))
    prior_patient = prior_pat[,c("m", "c", "tau_ctrl", "tau_notctrl", "probdiff")]
    
    ctrl_list = list(post=post_ctrl, postpred=postpred_ctrl, 
                     prior=prior_control, priorpred=priorpred_ctrl,
                     classif=classifs_ctrl)
    pat_list = list(post=post_pat, postpred=postpred_pat,
                    prior=prior_patient, priorpred=priorpred_pat,
                    classif=classifs_pat)
    
    return( list(ctrl=ctrl_list, pat=pat_list) )
  })
}

dat = read.csv("../Data_prepped.csv", header=TRUE)

cord = c("NDUFB8", "CYB", "MTCO1")
chlabs = c("CI", "CII","CIV")
names(chlabs) = cord
mitochan = "VDAC1"

ptsAll = unique(dat$patient_id)
pts = ptsAll[grepl("P", ptsAll)]

inputs = list()
{
  input0 = list()
  input0$MCMCOut = 2000
  input0$MCMCBurnin = 1000
  input0$MCMCThin = 50
  input0$n.chains = 1
  
  for(chan in cord){
    for(pat in pts){
      outroot = paste(chan, pat, sep="_")
      inputs[[outroot]] = input0
      inputs[[outroot]]$chan = chan
      inputs[[outroot]]$pat = pat
    } # pts
  } # chans
}

ncores = detectCores() - 1
cl  = makeCluster(ncores)
{
  clusterExport(cl, c("modelstring"))
  clusterEvalQ(cl, {
    library("rjags")
    source("helper_functions.R", local=TRUE)
  })
  linreg_output = parLapply(cl, inputs, inference)
}
stopCluster(cl)

for(outroot in names(linreg_output)){
  output_saver(linreg_output[[outroot]], folder, outroot)
}




