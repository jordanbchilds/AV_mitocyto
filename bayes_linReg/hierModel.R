
install.packages("devtools")
library("devtools")
install_github("jordanbchilds/analysis2Dmito")
library(analysis2Dmito)

library("parallel")
library("dplyr")
library("readr")
library("tidyr")
library("data.table")


folder = "analysis2Dmito_loglog"

dir.create(file.path("Output"), showWarnings = FALSE)
dir.create(file.path("Output", folder), showWarnings = FALSE)

raw_data = read.csv("../Data_prepped.csv", header=TRUE)

sbj = unique(raw_data$patient_id)
ctrlID = grep("C", sbj, value=TRUE)
pts = sbj[!(sbj %in% ctrlID)]

mitochan = "VDAC"
channels = c("NDUFB8", "CYB", "MTCO1")

raw_data = raw_data[,c("ID", "patient_id", mitochan, channels)]
colnames(raw_data) = c("fibreID", "sampleID", mitochan, channels)

data = pivot_longer(raw_data, cols=c(mitochan, channels), names_to="channel")
data = as.data.frame(data)

data$value = log(log(data$value))

grad = matrix(NA, nrow=length(channels), ncol=length(ctrlID))
colnames(grad) = ctrlID
rownames(grad) = channels

inter = grad
prec = grad

for( chan in channels ){
  for( crl in ctrlID ){
    x = data[data$channel==mitochan & data$sampleID==crl, "value"]
    y = data[data$channel==chan & data$sampleID==crl, "value"]
    

    dd = data.frame(mitochan=x, chan=y)
    
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

data_list = list()
for(chan in channels){
  for(pat in pts){
    data_list[[paste(chan, pat, sep = "__")]] = getData_mats(
      data,
      pts = pat,
      channels = c(mitochan, chan),
      ctrlID = ctrlID,
      getIndex = TRUE
    )
  }
}

mean_mu_m = mean(grad)
prec_mu_m = 1/0.05^2
curve(dnorm(x, mean_mu_m, 1/sqrt(prec_mu_m)), from=0.9, to=1.4, main="mu_m Prior")

mean_mu_c = mean(inter)
prec_mu_c = 1/0.05^2
curve(dnorm(x, mean_mu_c, 1/sqrt(prec_mu_c)), from=0, to=1, main="mu_c Prior")

tau_m_mode = 1 / var(c(grad))
tau_m_var = 1
rate_tau_m = 0.5*(tau_m_mode + sqrt(tau_m_mode^2+4*tau_m_var)) / tau_m_var
shape_tau_m = 1 + tau_m_mode*rate_tau_m
curve(dgamma(x, shape_tau_m, rate_tau_m), from=0, to=100, main="tau_m Prior")

tau_c_mode = 1 / var(c(inter))
tau_c_var = 1
rate_tau_c = 0.5*(tau_c_mode + sqrt(tau_c_mode^2+4*tau_c_var)) / tau_c_var
shape_tau_c = 1 + tau_c_mode*rate_tau_c
curve(dgamma(x, shape_tau_c, rate_tau_c), from=0, to=100, main="tau_c Prior")

tau_mode = mean(prec)
tau_var = 1
rate_tau = 0.5*(tau_mode + sqrt(tau_mode^2+4*tau_var)) / tau_var
shape_tau = 1 + tau_mode*rate_tau
curve(dgamma(x, shape_tau, rate_tau), from=0, to=100, main="tau Prior")

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
  tau_def=tau_def
)

ncores = 6
cl  = makeCluster(ncores)
{
  #clusterExport(cl, c("colQuantiles"))
  output = parLapply(
    cl,
    data_list,
    inference,
    MCMCthin = 10,
    parameterVals = paramVals
  )
}
stopCluster(cl)

for( root in names(output)){
  list_saver(output[[root]], file.path("Output", folder, root))
}







