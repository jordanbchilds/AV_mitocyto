# library("devtools")
# install_github("jordanbchilds/analysis2Dmito", force=TRUE)

install.packages("../../analysis2Dmito", repo=NULL, type="source")
library("analysis2Dmito")

library("stringr")
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

dat = pivot_longer(raw_data, cols=c(mitochan, channels), names_to="Channel", values_to="Value")
dat = as.data.frame(dat)
dat$Value = log(dat$Value)

dir.create(file.path("PDF"), showWarnings = FALSE)

bhlmm_predict = function(post, xSyn){
  nSbj = sum( grepl("m\\[", colnames(post)) )
  nSyn = length(xSyn)
  
  nSample = nrow(post)
  
  m = post[,paste0("m[", nSbj, "]")]
  c = post[,paste0("c[", nSbj, "]")]
  tau = post[,"tau_norm"]
  
  predQnt = matrix(NA, nrow=nSyn, ncol=3)
  
  for(i in 1:nSyn){
    pred = rnorm(nSample, mean= m*xSyn[i] + c, sd=1/sqrt(tau) )
    predQnt[i,] = quantile(pred, probs=c(0.025, 0.5, 0.975))
  }
  pred = cbind(xSyn, predQnt)
  colnames(pred) = c("mitochan", "lwrNorm", "medNorm", "uprNorm")
  return( pred )
}

# ------ MCMC plot 
pdf(file.path("PDF", "MCMCplot.pdf"), width=13, height=8)
{
    for(chan in channels){
      for(pat in pts) {
        for(ch in 1:5){
          outroot = paste0(chan, "_", pat, "_chain_", str_pad(ch, width=2, pad="0"))
          fn_post = file.path("Output", paste0(outroot, "_POST.csv"))
          if(file.exists(fn_post)){
            post = as.data.frame(fread(fn_post))
            
            fn_prior = file.path("Output", paste0(outroot, "_PRIOR.csv"))
            prior = as.data.frame(fread(fn_prior))
            analysis2Dmito::MCMCplot(post, prior, nRow=3, lag=100)
            title(main=outroot, outer=TRUE, line=-1)
          }
        }
      }
    }
}
dev.off()

# ------ postPlot
pdf(file.path("PDF", "postPlot.pdf"), width=13, height=8)
{
    for(chan in channels){
      for( pat in pts ){
        for(ch in 1:5){
          outroot = paste0(chan, "_", pat, "_chain_", str_pad(ch, width=2, pad="0"))
          fn_post = file.path("Output", paste0(outroot, "_POST.csv"))
          
          if(file.exists(fn_post)){
            post = as.data.frame(fread(fn_post))
            fn_prior = file.path("Output", paste0(outroot, "_PRIOR.csv"))
            prior = as.data.frame(fread(fn_prior))
            
            # fn_postpred = file.path("Output", paste0(outroot, "_POSTPRED.csv"))
            # postpred = as.data.frame(fread(fn_postpred))
            # 
            fn_class = file.path("Output", paste0(outroot, "_CLASSIF.csv"))
            class = apply(as.matrix(fread(fn_class)), 2, mean)
            
            dataMats = analysis2Dmito::getData_mats(data=dat, 
                                                    channels=c(mitochan, chan), 
                                                    ctrlID=ctrlID,
                                                    pts=pat)
            
            nSyn = 1000
            xSyn_range = range( c(dataMats$ctrl[,1], dataMats$pts[,1]) )
            xSyn = seq(xSyn_range[1], xSyn_range[2], length.out=nSyn)
            postpred = bhlmm_predict(post, xSyn)
            
            op = par(mfrow=c(3,3), mar=c(4,4,1,1), cex.main=2, cex.lab=1.5, cex.axis=1.5)
            analysis2Dmito::postPlot(post=post, 
                                     prior=prior,
                                     postpred=postpred, 
                                     classifs=class,
                                     dataMats = dataMats,
                                     var.names=c("mu_m", "tau_m", "tau_norm", "mu_c", "tau_c", "probdiff", "m", "c"),
                                     mitochan=paste0("log(", mitochan, ")"),
                                     chan=paste0("log(", chan, ")"))
            title(main=outroot, line=-1, outer=TRUE)
            par(op)
          }
        }
      }
    }
}
dev.off()

# ------ Classif plot
pdf(file.path("PDF", "classif.pdf"), width=13, height=8)
{
  op = par(mfrow=c(1,1), mar=c(6,6,2,2), cex.main=2, cex.axis=1.5, cex.lab=2)
    for(chan in channels){
      for( pat in pts ){
        for( ch in 1:10 ){
          outroot = paste0(chan, "_", pat, "_chain_", str_pad(ch, width=2, pad="0"))
          fn_class = file.path("Output", paste0(outroot, "_CLASSIF.csv"))
          
          if(file.exists(fn_class)){
            class = apply(as.matrix(fread(fn_class)), 2, mean)
            
            fn_post = file.path("Output", paste0(outroot, "_POST.csv"))
            post = as.data.frame( fread(fn_post))
            pi_est = round(mean( post[,"probdiff"] ), 3)
            
            fn_postpred = file.path("Output", folder, paste0(outroot, "_POSTPRED.txt"))
            postpred = as.data.frame(fread(fn_postpred))
            
            dataMats = analysis2Dmito::getData_mats(data=dat, 
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
