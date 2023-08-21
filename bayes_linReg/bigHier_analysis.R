
install.packages("devtools")
library("devtools")
install_github("jordanbchilds/analysis2Dmito")
library(analysis2Dmito)

library("parallel")
library("dplyr")
library("readr")
library("tidyr")
library("data.table")


folder = "analysis2Dmito_bigHier"

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

data$value = log(data$value)

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
    prec[chan, crl] = 1/summary(mod)$sigma
  }
}

grad_mean = apply(grad, 1, mean)
inter_mean = apply(inter, 1, mean)
prec_mean = apply(prec, 1, mean)

grad_var = apply(grad, 1, var)
inter_var = apply(inter, 1, var)
prec_var = apply(prec, 1, var)

data_list = list()
for (chan in channels) {
  data_list[[paste(chan)]] = getData_mats(
    data,
    pts = pts,
    channels = c(mitochan, chan),
    ctrlID = ctrlID,
    getIndex = TRUE
  )
}

tau_mode = 50
tau_var = 5
rate_tau = 0.5*(tau_mode + sqrt(tau_mode^2+4*tau_var)) / tau_var
shape_tau = 1 + tau_mode*rate_tau
tau_def = 0.01


ncores = length(channels)
cl  = makeCluster(ncores)
{
  output = parLapply(
    cl,
    data_list,
    inference,
    MCMCthin = 100,
    parameterVals = list(shape_tau = shape_tau, rate_tau = rate_tau, tau_def=tau_def)
  )
}
stopCluster(cl)

for( root in names(output) ){
  list_saver(output[[root]], file.path("Output", folder, root))
}
