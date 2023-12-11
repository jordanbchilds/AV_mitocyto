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

folder = "frequentist_linReg"

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


model_est = list() 
for (chan in channels) {
  tmp = matrix(NA, nrow=length(ctrlID), ncol=5)
  colnames(tmp) = c("m", "mse", "c", "cse", "tau")
  rownames(tmp) = ctrlID
  for( crl in ctrlID ){
    x = data[ data$sampleID==crl & data$channel==mitochan, "value"]
    y = data[ data$sampleID==crl & data$channel==chan, "value"]
    
    mod = lm(y ~ x)
    
    tmp[ crl, "m"] = mod$coefficients["x"]
    tmp[ crl, "mse"] = summary(mod)$coefficients["x", "Std. Error"]
    tmp[ crl, "c"] = mod$coefficients["(Intercept)"]
    tmp[ crl, "cse"] = summary(mod)$coefficients["(Intercept)", "Std. Error"]
    
    tmp[ crl, "tau"] = 1 / summary(mod)$sigma ^2
    model_est[[chan]] = tmp
  }
}

op = par(mfrow=c(1,1))
{
  for( chan in channels ){
    plot( NA, xlim=c(-1.0,2.4), ylim=c(0,5),
          xaxt=NULL, yaxt="n", xlab="", ylab="Control ID", main=chan)
    axis(side=2, at=1:4, labels=ctrlID)
    for(i in 1:length(ctrlID)){
      lines( model_est[[chan]][i, "m"]+c(-2*model_est[[chan]][i,"mse"], 2*model_est[[chan]][i,"mse"]), c(i,i), col="red", lwd=4)
      lines( model_est[[chan]][i, "c"]+c(-2*model_est[[chan]][i,"cse"], 2*model_est[[chan]][i,"cse"]), c(i,i), col="blue", lwd=4)
    }
    legend("topleft", lty=1, legend=c("m", "c"), col=c("red", "blue"), lwd=4)
  }
}
par(op)

chan = "NDUFB8"
chan_range = model_est[[chan]][,"m"] + 2*cbind(-model_est[[chan]][,"mse"], model_est[[chan]][,"mse"])
chan_range = model_est[[chan]][,"c"] + 2*cbind(-model_est[[chan]][,"cse"], model_est[[chan]][,"cse"])


for (chan in channels) {
  
  data_ctrl = getData_mats(data, channels=c(mitochan, chan), 
                           ctrlID = ctrlID, 
                           ctrl_only = TRUE,
                           getIndex = FALSE)
  
  df_ctrl = data.frame(mitochan=data_ctrl[,1], chan=data_ctrl[,2])
  
  mod = lm(chan~ mitochan, data=df_ctrl)
  
  df_syn = data.frame(mitochan=seq(min(data$value)-1, max(data$value)+1, length.out=1e3))
  
  model_pred = predict.lm(mod, newdata=df_syn, interval="prediction")
  model_pred = cbind( df_syn$mitochan, model_pred)
  colnames(model_pred) = c("mitochan", "fit", "lwr", "upr")
  
  write.table(model_pred, file=file.path("Output", folder, paste0(chan, "_POSTPRED.txt")),
              sep="\t", row.names=FALSE, col.names=TRUE)
  
  for(pat in pts ){
    data_pat = getData_mats(data, channels=c(mitochan, chan), 
                            ctrlID=ctrlID, 
                            ctrl_only = FALSE,
                            pts = pat, 
                            getIndex=FALSE)
    data_pat = data_pat$pts
    df_pat = data.frame(mitochan=data_pat[,1])
    pat_pred = predict.lm(mod, newdata=df_pat, interval="prediction")
    
    classif = ( data_pat[,2] > pat_pred[,"upr"] ) | ( data_pat[,2] < pat_pred[,"lwr"] )
    classif = as.numeric( classif )
    
    write.table(classif, file=file.path("Output", folder, paste0(chan,"__", pat, "_CLASSIF.txt")),
                sep="\t", row.names=FALSE, col.names=FALSE)
  }
}






