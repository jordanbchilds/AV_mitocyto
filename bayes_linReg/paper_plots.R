
library("analysis2Dmito")
# install.packages(c("data.table", "dplyr", "readr", "tidyr", "plyr", "devtools"))
library("data.table")
library("dplyr")
library("readr")
library("tidyr")
library("devtools")

# --- getting data

folder = "stan_sampler"

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
ctrlIDs = grep("C", sbj, value=TRUE)
pts = sbj[!(sbj %in% ctrlIDs)]

# --- log non-log data
png("./PDF/paper/controlData.png", width=9, height=6, units="in", res=300, pointsize=11)
{
  colVec = c(alphaCol(0.2, "lightpink"), 
             alphaCol(0.2, "purple3"), 
             alphaCol(0.2, "orange"), 
             alphaCol(0.2, "darkgreen")
  )
  data_exp = data
  data_exp$value = exp(data$value)
  op = par(mfrow=c(2,3), mar=c(5,5,1,1))
  for( chan in channels ){
    xlm = range( data_exp[data_exp$sampleID%in% ctrlIDs & data_exp$channel==mitochan, "value"] ) + c(0,1)
    ylm = range( data_exp[data_exp$sampleID%in% ctrlIDs & data_exp$channel!=mitochan, "value"] ) + c(0,1)
    
    plot(NA, xlim=xlm, ylim=ylm, 
         xlab=mitochan, ylab=chan)
    for( i in seq_along(ctrlIDs) ){
      points( data_exp[data_exp$sampleID==ctrlIDs[i] & data_exp$channel==mitochan, "value"], 
              data_exp[data_exp$sampleID==ctrlIDs[i] & data_exp$channel==chan, "value"],
              pch=20, col=colVec[i])
    }
    if(chan == channels[1])
      legend("topleft", 
             legend=ctrlIDs, col=c("lightpink", "purple3", "orange", "darkgreen"), 
             pch=20, 
             cex=1.0)
  }
  
  for( chan in channels ){
    xlm = range( data[data$sampleID%in% ctrlIDs & data$channel==mitochan, "value"] ) + c(0,0.5)
    ylm = range( data[data$sampleID%in% ctrlIDs & data$channel!=mitochan, "value"] ) + c(0,0.5)
    
    plot(NA, xlim=xlm, ylim=ylm, 
         xlab=paste0("log(", mitochan,")"), ylab=paste0("log(", chan, ")") )
    for( i in seq_along(ctrlIDs) ){
      points( data[data$sampleID==ctrlIDs[i] & data$channel==mitochan, "value"], 
              data[data$sampleID==ctrlIDs[i] & data$channel==chan, "value"],
              pch=20, col=colVec[i])
    }
    if(chan == channels[1]) 
      legend("topleft", 
             legend=ctrlIDs, col=c("lightpink", "purple3", "orange", "darkgreen"), 
             pch=20, 
             cex=1.0, 
             )
  }
  par(op)
}
dev.off()

# --- patient data example
pat = "P09"

png("./PDF/paper/pat_example.png", width=9, height=3, units="in", res=300, pointsize=11)
{
  op = par(mfrow=c(1,3), mar=c(5,5,2,1))
  for( chan in channels ){
    xlm = range( data[data$sampleID%in% c(ctrlIDs, pat) & data$channel==mitochan, "value"] ) + c(0,0.5)
    ylm = range( data[data$sampleID%in% c(ctrlIDs, pat) & data$channel!=mitochan, "value"] ) + c(0,0.5)
    
    plot(NA, xlim=xlm, ylim=ylm, 
         xlab=paste0("log(", mitochan,")"), ylab=paste0("log(", chan, ")") )
    points( data[data$sampleID%in% ctrlIDs & data$channel==mitochan, "value"],
            data[data$sampleID%in% ctrlIDs & data$channel==chan, "value"],
            col=alphaBlack(0.05), bg=alphaBlack(0.05), 
            pch=21, cex=0.5)
    points( data[data$sampleID==pat & data$channel==mitochan, "value"],
            data[data$sampleID==pat & data$channel==chan, "value"],
            pch=21, col=alphaCol(0.1, "hotpink"), bg=alphaCol(0.1, "hotpink"))
    
    if(chan == channels[1])
      legend("topleft", 
             col=c(alphaBlack(0.7), alphaCol(0.7, "hotpink")), 
             pch=20, 
             legend=c("Control", "patient"))
  }
  par(op)
}
dev.off()

# --- frequentist classification
png("./PDF/paper/freq_example.png", width=9, height=3, units="in", res=300, pointsize=11)
{
  colVec = c(alphaBlue(0.5), alphaRed(0.5))
  
  op = par(mfrow=c(1,3), mar=c(5,5,2,1))
  for( chan in channels ){
    freq_class = as.data.frame( fread( paste0("./Output/frequentist_linReg/", chan, "__", pat, "_CLASSIF.txt") ) )
    freq_postpred = as.data.frame( fread( paste0("./Output/frequentist_linReg/", chan, "_POSTPRED.txt" )))
    
    xlm = range( data[data$sampleID%in% c(ctrlIDs, pat) & data$channel==mitochan, "value"] ) + c(0,0.5)
    ylm = range( data[data$sampleID%in% c(ctrlIDs, pat) & data$channel!=mitochan, "value"] ) + c(0,0.5)
    
    pi_est = round(colMeans(freq_class), 3)
    
    plot(NA, xlim=xlm, ylim=ylm, 
         xlab=paste0("log(", mitochan,")"), ylab=paste0("log(", chan, ")"),
         main=bquote(pi*"="*.(pi_est)))
    points( data[data$sampleID%in% ctrlIDs & data$channel==mitochan, "value"],
            data[data$sampleID%in% ctrlIDs & data$channel==chan, "value"],
            pch=21, col=alphaBlack(0.05), bg=alphaBlack(0.05), cex=0.5)
    
    points( data[data$sampleID==pat & data$channel==mitochan, "value"],
            data[data$sampleID==pat & data$channel==chan, "value"],
            pch=21, col=colVec[freq_class[[1]]+1], bg=colVec[freq_class[[1]]+1])
    
    lines(freq_postpred$mitochan, freq_postpred$fit, 
          col=alphaGreen(1.0), lty=1, lwd=2)
    lines(freq_postpred$mitochan, freq_postpred$lwr, 
          col=alphaGreen(1.0), lty=2, lwd=2)
    lines(freq_postpred$mitochan, freq_postpred$upr, 
          col=alphaGreen(1.0), lty=2, lwd=2)
    
    if( chan == channels[1] )
      legend("topleft", 
             col=c(alphaBlack(0.7), alphaBlue(0.7), alphaRed(0.7)),
             pch=20, 
             legend=c("Control", "Healthy pat.", "Deficient pat.")
             )
  }
  par(op)
}
dev.off()

# --- bayes classification example
png("./PDF/paper/bayes_example.png", width=9, height=3, units="in", res=300, pointsize=11)
{
  op = par(mfrow=c(1,3), mar=c(5,5,2,1))
  for( chan in channels ){
    bayes_class_mat = as.matrix( fread( paste0("./Output/stan_sampler/", chan, "__", pat, "_CLASSIF.txt") ) )
    bayes_class = colMeans( bayes_class_mat )
    
    bayes_postpred = as.data.frame( fread( paste0("./Output/stan_sampler/", chan, "__", pat, "_POSTPRED.txt" )))
    
    xlm = range( data[data$sampleID%in% c(ctrlIDs, pat) & data$channel==mitochan, "value"] ) + c(0,0.5)
    ylm = range( data[data$sampleID%in% c(ctrlIDs, pat) & data$channel!=mitochan, "value"] ) + c(0,0.5)
    
    pi_est = round( colMeans( fread( paste0("./Output/stan_sampler/", chan, "__", pat, "_POST.txt" ))[, "probdiff"] ), 3)
    
    plot(NA, xlim=xlm, ylim=ylm, 
         xlab=paste0("log(", mitochan,")"), ylab=paste0("log(", chan, ")"),
         main=bquote("E("*pi*"|Y)="*.(pi_est)))
    points( data[data$sampleID%in% ctrlIDs & data$channel==mitochan, "value"],
            data[data$sampleID%in% ctrlIDs & data$channel==chan, "value"],
            pch=20, col=alphaBlack(0.05), cex=0.5)
    
    points( data[data$sampleID==pat & data$channel==mitochan, "value"],
            data[data$sampleID==pat & data$channel==chan, "value"],
            pch=20, col=classcols(bayes_class, alphaLevel=0.5))
    
    lines(bayes_postpred$mitochan, bayes_postpred$medNorm, 
          col=alphaGreen(1.0), lty=1, lwd=2)
    lines(bayes_postpred$mitochan, bayes_postpred$lwrNorm, 
          col=alphaGreen(1.0), lty=2, lwd=2)
    lines(bayes_postpred$mitochan, bayes_postpred$uprNorm, 
          col=alphaGreen(1.0), lty=2, lwd=2)
    
    if( chan == channels[1] )
      legend("topleft", 
             col=c(alphaBlack(0.7), alphaBlue(0.7), alphaRed(0.7)),
             pch=20, 
             legend=c("Control", "Healthy pat.", "Deficient pat.")
      )
  }
  par(op)
}
dev.off()

# --- prior-posteriors all n one
png("./PDF/paper/prior_post.png", width=9, height=3*4, units="in", res=300, pointsize=11)
{
  pat = "P09"
  chan = "NDUFB8"
  
  op = par(mfrow=c(4,2), mar=c(5,5,1,1))
  
  post = as.data.frame( fread( paste0("./Output/stan_sampler/", chan, "__", pat, "_POST.txt") ) )
  prior = as.data.frame( fread( paste0("./Output/stan_sampler/", chan, "__", pat, "_PRIOR.txt") ) )
  
  for(para in c("probdiff", "tau_norm")){
    dpost = density(post[,para])
    dprior = density(prior[,para])
    plot(NA, 
         xlim=range(c(dpost$x, dprior$x)), ylim=range(c(dpost$y, dprior$y)),
         ylab="Density", xlab=xlb[[para]]
    )
    if(para == "probdiff") curve(dunif(x), col=alphaPink(0.9), add=TRUE)
    else lines(dprior, col=alphaPink(0.9) )
    lines(dpost, col=alphaGreen(0.9) )
  }
  
  xlb = list(mu_m=bquote(mu[m]), mu_c=bquote(mu[c]), tau_m=bquote(tau[m]), tau_c=bquote(tau[c]))
  for(para in c("mu_m", "mu_c", "tau_m", "tau_c")){
    dpost = density(post[,para])
    dprior = density(prior[,para])
    plot(NA, 
         xlim=range(c(dpost$x, dprior$x)), ylim=range(c(dpost$y, dprior$y)),
         ylab="Density", xlab=xlb[[para]]
    )
    lines(dprior, col=alphaPink(0.9))
    lines(dpost, col=alphaGreen(0.9))
  }
    xlb = list(m="Slope", c="Intercept")
    for(para in c("m", "c")){
      xlm = NULL
      ylm = NULL
      
      dpost_list = list()
      for(i in 1:(length(ctrlIDs)+1) ){
        dpost_list[[i]] = density( post[, paste0(para, "[", i, "]")] )
        xlm = range( c(xlm, range(dpost_list[[i]]$x) ))
        ylm = range( c(ylm, range(dpost_list[[i]]$y) ))
      }
      
      plot(NA, 
           xlim=xlm, ylim=ylm,
           ylab="Density", xlab=xlb[[para]]
      )
      lines(density(prior[,paste0(para, "_pred")]), lty=3, col=alphaPink(1.0))
      lines(density(post[,paste0(para, "_pred")]), lty=3, col=alphaGreen(1.0))
      
      for(i in 1:length(ctrlIDs) ){
        lines( dpost_list[[i]],  lty=1, col=alphaGreen(0.4))
      }
      lines( dpost_list[[length(ctrlIDs)+1]], lty=1, col=alphaGreen(1.0), )
    }
  par(op)
}
dev.off()

# --- different gamma
png("./PDF/paper/postComp_gamma.png", width=9, height=3, units="in", res=300, pointsize=11)
{
  op = par(mfrow=c(1,3), mar=c(5,5,1,1))
  chan = "NDUFB8"
  pat = "P09"
  root = paste0(chan, "_", pat)
  
  tau_defs = c("000001", "00001", "0001", "001", "01", "10")
  
  post_list = list(m = list(), 
                   c = list(), 
                   tau = list())
  
  m_ylm = NULL; c_ylm = NULL; tau_ylm = NULL;
  m_xlm = NULL; c_xlm = NULL; tau_xlm = NULL;
  
  for(i in 1:length(tau_defs)){
    tau_str = as.character(tau_defs[i])
    fn_post = file.path("Output", 
                        paste0("stan_sampler_gamma", gsub("\\.", "", tau_str)),
                        paste0(root, "_POST.txt"))
    post = as.data.frame( fread(fn_post) )
    post_list$m[[i]] = density( post[,"m[5]"] )
    post_list$c[[i]] = density( post[,"c[5]"] )
    post_list$tau[[i]] = density( post[,"tau_norm"] )
    
    m_ylm = range(c(m_ylm, range(post_list$m[[i]]$y)))
    c_ylm = range(c(c_ylm, range(post_list$c[[i]]$y)))
    tau_ylm = range(c(tau_ylm, range(post_list$tau[[i]]$y)))
    
    m_xlm = range(c(m_xlm, range(post_list$m[[i]]$x)))
    c_xlm = range(c(c_xlm, range(post_list$c[[i]]$x)))
    tau_xlm = range(c(tau_xlm, range(post_list$tau[[i]]$x)))
  }
  
  plot(NA, xlim=m_xlm, ylim=m_ylm, xlab="Slope", ylab="Density")
  for(i in 1:length(tau_defs)){
    lines( post_list$m[[i]], col=i)
  }
  legend("topleft", 
         lty=1, col=1:length(tau_defs), 
         legend=tau_defs)
  
  plot(NA, xlim=c_xlm, ylim=c_ylm, xlab="Intercept", ylab="Density")
  for(i in 1:length(tau_defs)){
    lines( post_list$c[[i]], col=i)
  }
  
  plot(NA, xlim=tau_xlm, ylim=tau_ylm, xlab="tau", ylab="Density")
  for(i in 1:length(tau_defs)){
    lines( post_list$tau[[i]], col=i)
  }
}
dev.off()



# --- freq, bayes, manual comparison of pi

mdat = as.data.frame( fread( "./dat_manClass_prepped.txt") )

pi_list = list()
piFreq_list = list()
piMan_list = list()
for(pat in pts){
  for( chan in channels ){
    root = paste0(chan, "__", pat)
    fn_post = file.path("Output", folder, paste0(root, "_POST.txt"))
    post = as.matrix( fread( fn_post ))
    pi_list[[root]] = post[, "probdiff"]
    fn_freqClass = file.path("Output", "frequentist_linReg", paste0(root, "_CLASSIF.txt"))
    freqClass = as.matrix( fread(fn_freqClass) )
    piFreq_list[[root]] = mean( freqClass )
    manClass = mdat[ mdat$channel==chan & mdat$patient_id==pat, grep("Class", colnames(mdat))]
    piMan_list[[root]] = mean(manClass)
  }
}

png("./PDF/paper/proportion_comparison.png", width=9, height=6, units="in", res=300, pointsize=11)
{
  op = par(mar=c(4,4,1,9), xpd=TRUE)
  
  colVec = c(alphaCol(0.05, name="lightseagreen"), 
             alphaCol(0.05, name="deeppink3"), 
             alphaCol(0.05, name="darkorange"))
    
  myAt = c(1:44)[ 1:44 %% 4 !=0 ]
  myTicks = 4* 1:11 - 2
  plot(NA, xlim=c(1,43), ylim=c(0,1), xaxt="n", 
       ylab="Deficient Proportion", xlab="Patient")
  
  for(i in seq(0,10, by=2)) rect( xleft=i*4, xright=(i+1)*4, ybottom=par("usr")[3], ytop=par("usr")[4], border=NA, col=alphaBlack(0.1))
    stripchart(pi_list, vertical=TRUE, pch=20, cex=0.5, 
               method="jitter", main="", ylim=c(0,1), xlim=c(1,43),
               xlab="Patient", ylab="Deficient Proportion", 
               col=rep(colVec, length(pi_list)), 
               at=myAt, xaxt="n", add=TRUE)
    axis(side=1, at=myTicks, labels=pts)
    points(myAt, unlist(piFreq_list), pch=23, 
           bg=rep(c(alphaCol(0.7, name="lightseagreen"), 
                    alphaCol(0.7, name="deeppink3"), 
                    alphaCol(0.7, name="darkorange")), length(pi_list)), 
           col="black",
           cex=0.8)
    points(myAt, unlist(piMan_list), pch=25, 
           bg=rep(c(alphaCol(0.7, name="lightseagreen"), 
                    alphaCol(0.7, name="deeppink3"), 
                    alphaCol(0.7, name="darkorange")), length(pi_list)),
           col="black",
           cex=0.8)
    legend("topleft", 
           inset=c(1.025, 0), 
           legend=c("Bayesian est.", "Frequentist est.", "Mannual est."), 
           pch=c(21, 23, 25), 
           col=c(alphaBlack(1.0), alphaBlack(1.0), alphaBlack(1.0)))
    legend("topleft", 
           inset=c(1.025, 0.17), 
           legend=channels, 
           pch=20, 
           col=c(alphaCol(0.7, name="lightseagreen"), 
                 alphaCol(0.7, name="deeppink3"), 
                 alphaCol(0.7, name="darkorange")))
    par(op)
}
dev.off()

# --- better comparison of classifs

mdat = as.data.frame( fread( file.path("../dat_with_class_prepped.txt") ) )

bayes_diff = list()
freq_diff = list()
for( chan in channels ){
  for( pat in pts ){
    root = paste0(chan, "__", pat)
    fn_post = file.path("Output", folder, paste0(root, "_POST.txt"))
    probdiff_post = as.data.frame( fread( fn_post ))[,"probdiff"]
    
    man_est = mean( mdat[ mdat$patient_id==pat & mdat$channel==chan, "jbcClass"] )
    
    fn_freq = file.path("Output", "frequentist_linReg", paste0(root, "_CLASSIF.txt"))
    freq_class = as.data.frame( fread( fn_freq ) )
    
    bayes_diff[[root]] = probdiff_post - man_est
    freq_diff[[root]] = mean( freq_class[,1] )
  }
}

stripchart(bayes_diff, vertical=TRUE, pch=20, cex=0.2, col=alphaBlack(0.01), method="jitter",
           ylim=c(-0.05, 1.0))
stripchart(freq_diff, vertical=TRUE, pch=23, cex=0.4, col=alphaRed(1.0), 
           method="jitter", add=TRUE)
abline( h=0, col="red", lty=3)

png("./PDF/paper/proportion_comparison_better.png", width=9, height=6, units="in", res=300, pointsize=11)
{
  op = par(mar=c(4,4,1,9), xpd=TRUE)
  
  colVec = c(alphaCol(0.05, name="lightseagreen"), 
             alphaCol(0.05, name="deeppink3"), 
             alphaCol(0.05, name="darkorange"))
  
  myAt = c(1:44)[ 1:44 %% 4 !=0 ]
  myTicks = 4* 1:11 - 2
  plot(NA, xlim=c(1,43), ylim=c(-0.1,1), xaxt="n", 
       ylab="Difference from manual classification", xlab="Patient")
  
  for(i in seq(0,10, by=2)) rect( xleft=i*4, xright=(i+1)*4, ybottom=par("usr")[3], ytop=par("usr")[4], border=NA, col=alphaBlack(0.1))
  stripchart(bayes_diff, vertical=TRUE, pch=20, cex=0.5, 
             method="jitter", main="", ylim=c(0,1), xlim=c(1,43),
             xlab="Patient", ylab="Deficient Proportion", 
             col=rep(colVec, length(bayes_diff)), 
             at=myAt, xaxt="n", add=TRUE)
  axis(side=1, at=myTicks, labels=pts)
  points(myAt, unlist(freq_diff), pch=17, cex=1.2,
         col=rep(c(alphaCol(0.9, name="lightseagreen"), 
                   alphaCol(0.9, name="deeppink3"), 
                   alphaCol(0.9, name="darkorange")), length(freq_diff)))
  lines(x=c(0,44), c(0,0), lty=3, col="black", lwd=3)
  legend("topleft", 
         inset=c(1.025, 0), 
         legend=channels, 
         title="Bayesian",
         pch=20, 
         col=c(alphaCol(0.7, name="lightseagreen"), 
               alphaCol(0.7, name="deeppink3"), 
               alphaCol(0.7, name="darkorange")))
  legend("topleft", 
         inset=c(1.025, 0.21), 
         legend=channels, 
         title = "Frequentist",?layo
         pch=17, 
         col=c(alphaCol(0.7, name="lightseagreen"), 
               alphaCol(0.7, name="deeppink3"), 
               alphaCol(0.7, name="darkorange")) )
  par(op)
}
dev.off()






