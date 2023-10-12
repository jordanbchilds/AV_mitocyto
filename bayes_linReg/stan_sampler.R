# install.packages("devtools")
# library("devtools")
# install_github("jordanbchilds/analysis2Dmito", force=TRUE)

library(analysis2Dmito)
library("parallel")
library("dplyr")
library("readr")
library("tidyr")
library("data.table")
library("rstan")
source("stan_sampler_function.R", local=TRUE)

folder = "stan_sampler_multiChain"

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
    xSyn = data.frame(mitochan=seq(min(data$value)-1, max(data$value)+1, length.out=1e3))
    
    mod = lm(chan ~ mitochan, data=dd)
    pred = predict.lm(mod, newdata=xSyn, interval="prediction")
    for( pat in pts ){
      # xPat = data[data$channel==mitochan & data$sampleID==pat, "value"]
      # yPat = data[data$channel==chan & data$sampleID==pat, "value"]
      # 
      # plot(xCtrl,yCtrl, pch=20, col=alphaBlack(1.0),
      #      xlab=mitochan, ylab=chan,
      #      xlim=range(data$value), ylim=range(data$value))
      # points(xPat, yPat, pch=20, col=alphaBlue(0.7))
      # lines(xSyn$mitochan, pred[,"fit"], lty="solid", col="pink", lwd=3)
      # lines(xSyn$mitochan, pred[,"lwr"], lty="dashed", col="pink", lwd=3)
      # lines(xSyn$mitochan, pred[,"upr"], lty="dashed", col="pink", lwd=3)
    }
      grad[chan, crl] = mod$coefficients[2]
      inter[chan, crl] = mod$coefficients[1]
      prec[chan, crl] = 1 / summary(mod)$sigma^2
  }
}

grad_mean = apply(grad, 1, mean)
inter_mean = apply(inter, 1, mean)
prec_mean = apply(prec, 1, mean)

# grad_var = apply(grad, 1, var)
# inter_var = apply(inter, 1, var)
# prec_var = apply(prec, 1, var)

tau_c_vars = rep(50, nChan)^2 
names(tau_c_vars) = channels

tau_m_vars = rep(50, nChan)^2
names(tau_m_vars) = channels

tau_vars = rep(10, nChan)^2 
names(tau_vars) = channels

grad_var =  1 / rep(50, nChan)
inter_var = 1 / rep(50, nChan)
names(grad_var) = channels
names(inter_var) = channels

data_list = list()
for (chan in channels) {
  mean_mu_m = grad_mean[chan]
  prec_mu_m = 1/0.25^2
  
  mean_mu_c = inter_mean[chan]
  prec_mu_c = 1/0.25^2
  
  tau_m_mode = 1 / grad_var[chan]
  tau_m_var = tau_m_vars[chan] 
  rate_tau_m = 0.5*(tau_m_mode + sqrt(tau_m_mode^2+4*tau_m_var)) / tau_m_var
  shape_tau_m = 1 + tau_m_mode*rate_tau_m
  
  tau_c_mode = 1 / inter_var[chan]
  tau_c_var = tau_c_vars[chan] 
  rate_tau_c = 0.5*(tau_c_mode + sqrt(tau_c_mode^2+4*tau_c_var)) / tau_c_var
  shape_tau_c = 1 + tau_c_mode*rate_tau_c
  
  tau_mode = prec_mean[chan]
  tau_var = tau_vars[chan] 
  rate_tau = 0.5*(tau_mode + sqrt(tau_mode^2+4*tau_var)) / tau_var
  shape_tau = 1 + tau_mode*rate_tau
  
  tau_def = 0.01
  
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
    alpha_pi = 0.0,
    beta_pi = 0.5,
    tau_def=tau_def,
    slope_lb = 0.1
  )
  
  for (i in seq_along(pts)) {
    data_list[[paste(chan, pts[i], sep = "__")]] = getData_mats(
      data,
      pts = pts[i],
      channels = c(mitochan, chan),
      ctrlID = ctrlID,
      getIndex = TRUE
    )
    data_list[[paste(chan, pts[i], sep="__")]]$parameterVals = paramVals
  }
}

ncores = detectCores() - 1
cl  = makeCluster(ncores)
{
  clusterExport(cl, c("stan"))
  output = parLapply(
    cl,
    data_list,
    stan_inference,
    MCMCburnin = 5000,
    MCMCout = 4000,
    MCMCthin=1,
    nChains=5,
    max_logLik=FALSE,
    nCores=1
  )
}
stopCluster(cl)


for( rt in names(output) ){
  list_saver(output[[rt]], file.path("Output", folder, rt))
}




