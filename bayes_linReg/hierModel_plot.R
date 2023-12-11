# library("devtools")
# install_github("jordanbchilds/analysis2Dmito", force=TRUE)

library("analysis2Dmito")
library("tidyr")
library("data.table")

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

folder = "stan_sampler"

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
          post = post[,colnames(post)!="log_lik"]
          
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
                                   postpred=postpred[1:1000,], 
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
pdf(file.path("PDF", folder, "classif.pdf"), width=13, height=8)
{
  op = par(mfrow=c(1,1), mar=c(6,6,6,3), cex.main=2, cex.axis=1.5, cex.lab=2)
    for(chan in channels){
      for( pat in pts ){
        outroot = paste(chan, pat, sep="__")
        fn_class = file.path("Output", folder, paste0(outroot, "_CLASSIF.txt"))
        
        if(file.exists(fn_class)){
          class = apply(as.matrix(fread(fn_class)), 2, mean)
          
          fn_post = file.path("Output", folder, paste0(outroot, "_POST.txt"))
          post = as.data.frame( fread(fn_post))
          pi_est = round(mean( post[,"probdiff"] ), 3)
          
          fn_postpred = file.path("Output", folder, paste0(outroot, "_POSTPRED.txt"))
          postpred = as.data.frame(fread(fn_postpred))
          
          dataMats = analysis2Dmito::getData_mats(data=data, 
                                                  channels=c(mitochan, chan), 
                                                  ctrlID=ctrlID,
                                                  pts=pat)
          
          analysis2Dmito::classif_plot(dataMats=dataMats,
                                       postpred=NULL, 
                                       classifs=class,
                                       xlab="VDAC1", 
                                       ylab="NDUFB8")
          title(main=bquote( atop(.(pat), 
                                  "nPat:"~.(nrow(dataMats$pts))*",   E("*pi*"|X)="*.(pi_est))))
        }
      }
    }
  par(op)
}
dev.off()

# ------ Prob deficient density
pdf(file.path("PDF", folder, "probDef.pdf"), width=8, height=5)
{
  op = par(mfrow=c(1,1), mar=c(6,6,3,3), cex.main=2, cex.axis=1.5, cex.lab=2)
  for(chan in channels){
    for( pat in pts ){
      outroot = paste(chan, pat, sep="__")
      fn_class = file.path("Output", folder, paste0(outroot, "_CLASSIF.txt"))
      if(file.exists(fn_class)){
        class_mat = as.matrix(fread(fn_class))
        class = apply(class_mat, 2, mean)
        
        dd = density(class)
        domain = dd$x > 0.0 & dd$x<1.0
        plot( dd$x[domain] ,dd$y[ domain ], type='l', 
              main=paste(chan, pat), ylab="Density", xlab="Probability Deficient")
      }
    }
  }
  par(op)
}
dev.off()






