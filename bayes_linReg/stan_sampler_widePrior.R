
library(analysis2Dmito)
library("parallel")
library("dplyr")
library("readr")
library("tidyr")
library("data.table")
library("rstan")

folder = "stan_sampler_wide2"

dir.create(file.path("Output"), showWarnings = FALSE)
dir.create(file.path("Output", folder), showWarnings = FALSE)

mitochan = "VDAC"
channels = c("NDUFB8", "CYB", "MTCO1")
nChan = length(channels)

raw_data = read.csv("../Data_prepped.csv", header=TRUE)

raw_data = raw_data[,c("ID", "patient_id", mitochan, channels)]
colnames(raw_data) = c("fibreID", "sampleID", mitochan, channels)

data_lng = pivot_longer(raw_data, cols=c(mitochan, channels), names_to="channel")
data_lng = as.data.frame(data_lng)

data = data_lng
data$value = log(data$value)

sbj = unique(data$sampleID)
ctrlID = grep("C", sbj, value=TRUE)
pts = sbj[!(sbj %in% ctrlID)]

priorInflate = 10
{
  grad = matrix(NA, nrow=length(channels), ncol=length(ctrlID))
  colnames(grad) = ctrlID
  rownames(grad) = channels
  
  inter = grad
  prec = grad
  
  for( chan in channels ){
    for( crl in ctrlID ){
      xCtrl = data[data$channel==mitochan & data$sampleID==crl, "value"]
      yCtrl = data[data$channel==chan & data$sampleID==crl, "value"]
      dd = data.frame(mitochan=xCtrl, chan=yCtrl)
      xSyn = data.frame(mitochan=seq(min(data$value)-2, max(data$value)+2, length.out=1e3))
      
      mod = lm(chan ~ mitochan, data=dd)
      grad[chan, crl] = mod$coefficients[2]
      inter[chan, crl] = mod$coefficients[1]
      prec[chan, crl] = 1 / summary(mod)$sigma^2
      
    }
  }
  
  grad_mean = apply(grad, 1, mean)
  inter_mean = apply(inter, 1, mean)
  prec_mean = apply(prec, 1, mean)
  
  grad_var = apply(grad, 1, var)
  inter_var = apply(inter, 1, var)
  prec_var = apply(prec, 1, var)
  
  tau_c_vars = rep(50, nChan)^2 *priorInflate
  names(tau_c_vars) = channels
  
  tau_m_vars = rep(50, nChan)^2 *priorInflate
  names(tau_m_vars) = channels
  
  tau_vars = rep(10, nChan)^2 *priorInflate
  names(tau_vars) = channels
  
  grad_prec =  rep(50, nChan)
  inter_prec = rep(50, nChan)
  names(grad_prec) = channels
  names(inter_prec) = channels
}

for (chan in channels) {
  mean_mu_m = grad_mean[chan]
  prec_mu_m = 1/0.25^2 * (1/priorInflate)
  
  mean_mu_c = inter_mean[chan]
  prec_mu_c = 1/0.25^2 * (1/priorInflate)
  
  tau_m_mode = grad_prec[chan]
  tau_m_var = tau_m_vars[chan]
  rate_tau_m = 0.5*(tau_m_mode + sqrt(tau_m_mode^2+4*tau_m_var)) / tau_m_var
  shape_tau_m = 1 + tau_m_mode*rate_tau_m
  
  tau_c_mode = inter_prec[chan]
  tau_c_var = tau_c_vars[chan]
  rate_tau_c = 0.5*(tau_c_mode + sqrt(tau_c_mode^2+4*tau_c_var)) / tau_c_var
  shape_tau_c = 1 + tau_c_mode*rate_tau_c
  
  tau_mode = prec_mean[chan]
  tau_var = tau_vars[chan]
  rate_tau = 0.5*(tau_mode + sqrt(tau_mode^2+4*tau_var)) / tau_var
  shape_tau = 1 + tau_mode*rate_tau
  
  tau_def = 0.0001
  
  paramVals = list(
    mean_mu_m=mean_mu_m, 
    prec_mu_m=prec_mu_m, 
    mean_mu_c=mean_mu_c, 
    prec_mu_c=prec_mu_c,
    shape_tau_m=shape_tau_m, 
    rate_tau_m=rate_tau_m, 
    shape_tau_c=shape_tau_c, 
    rate_tau_c=rate_tau_c, 
    shape_tau = shape_tau,
    rate_tau = rate_tau, 
    alpha_pi = 1.0,
    beta_pi = 1.0,
    tau_def=tau_def
  )
  
  data_list = list()
  for (i in seq_along(pts)) {
    rt = paste(chan, pts[i], sep="__")
    data_list[[rt]] = getData_mats(
      data,
      pts = pts[i],
      channels = c(mitochan, chan),
      ctrlID = ctrlID,
      getIndex = TRUE
    )
  }
  
  # ncores = 8
  ncores = detectCores() - 4
  cl  = makeCluster(ncores)
  {
    output = parLapply(
      cl,
      data_list,
      stan_inference,
      warmup=50000,
      iter=52000,
      parameterVals = paramVals
    )
  }
  stopCluster(cl)

  for( rt in names(output) ){
    list_saver(output[[rt]], file.path("Output", folder, rt), rootSep="_")
  }
}


