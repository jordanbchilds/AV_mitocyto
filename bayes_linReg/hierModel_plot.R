library("devtools")
install_github("jordanbchilds/analysis2Dmito", force=TRUE)
library(analysis2Dmito)

library("parallel")

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


folder = "analysis2Dmito_loglog"

dir.create(file.path("PDF"), showWarnings = FALSE)
dir.create(file.path("PDF", folder), showWarnings = FALSE)

# ------ MCMC plot 
pdf(file.path("PDF", folder, "MCMCplot.pdf"), width=13, height=8)
{
    for(chan in channels){
      for(pat in pts) {
        outroot = paste(chan, pat, sep="__")
        fn_post = file.path("Output", folder, paste0(outroot, "_POST.txt"))
        if(file.exists(fn_post)){
          post = as.data.frame(fread(fn_post))
          fn_prior = file.path("Output", folder, paste0(outroot, "_PRIOR.txt"))
          prior = as.data.frame(fread(fn_prior))
          analysis2Dmito::MCMCplot(post, prior, nRow=3, lag=100)
          title(main=paste(chan, pat), outer=TRUE, line=-1)
        }
      }
    }
}
dev.off()

# ------ postPlot
pdf(file.path("PDF", folder, "postPlot.pdf"), width=13, height=8)
{
    for(chan in channels){
      for( pat in pts ){
        
        outroot = paste(chan, pat, sep="__")
        fn_post = file.path("Output", folder, paste0(outroot, "_POST.txt"))
        if(file.exists(fn_post)){
          post = as.data.frame(fread(fn_post))
          fn_prior = file.path("Output", folder, paste0(outroot, "_PRIOR.txt"))
          prior = as.data.frame(fread(fn_prior))
          fn_postpred = file.path("Output", folder, paste0(outroot, "_POSTPRED.txt"))
          postpred = as.data.frame(fread(fn_postpred))
          
          fn_class = file.path("Output", folder, paste0(outroot, "_CLASSIF.txt"))
          class = apply(as.matrix(fread(fn_class)), 2, mean)
          
          
          dataMats = analysis2Dmito::getData_mats(data=data, 
                                                  channels=c(mitochan, chan), 
                                                  ctrlID=ctrlID,
                                                  pts=pat)
          
          op = par(mfrow=c(3,3), mar=c(4,4,3,3), cex.main=2, cex.lab=1.5, cex.axis=1.5)
          analysis2Dmito::postPlot(post=post, prior=prior,
                                   postpred=postpred, 
                                   classifs=class,
                                   dataMats = dataMats,
                                   var.names=c("mu_m", "tau_m", "tau_norm", "mu_c", "tau_c", "probdiff", "m", "c"),
                                   mitoPlot_xlab=paste0("log(", mitochan, ")"),
                                   mitoPlot_ylab=paste0("log(", chan, ")"))
          title(main=paste(chan, pat), outer=TRUE, line=-2)
          par(op)
        }
      }
    }
}
dev.off()

# ------ Classif plot
pdf(file.path("PDF", folder, "rawExpression_classif.pdf"), width=8, height=5)
{
  op = par(mfrow=c(1,1), mar=c(6,6,3,3), cex.main=2, cex.axis=1.5, cex.lab=2)
    for(chan in channels){
      for( pat in pts ){
        outroot = paste(chan, pat, sep="__")
        fn_class = file.path("Output", folder, paste0(outroot, "_CLASSIF.txt"))
        if(file.exists(fn_class)){
          class = apply(as.matrix(fread(fn_class)), 2, mean)
          
          fn_postpred = file.path("Output", folder, paste0(outroot, "_POSTPRED.txt"))
          postpred = as.data.frame(fread(fn_postpred))
          
          rawDat = data
          logDat = data
          
          rawDat$value = exp(exp(rawDat$value))
          logDat$value = exp(logDat$value)
          
          dataMats = analysis2Dmito::getData_mats(data=rawDat, 
                                                  channels=c(mitochan, chan), 
                                                  ctrlID=ctrlID,
                                                  pts=pat)
          
          analysis2Dmito::classif_plot(dataMats=dataMats,
                                       postpred=NULL, 
                                       classifs=class,
                                       xlab="VDAC1", 
                                       ylab="NDUFB8")
        }
      }
    }
  par(op)
}
dev.off()
