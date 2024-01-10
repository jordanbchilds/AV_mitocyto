
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
  pat = "P09"
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

# --- prior-posteriors all in one
png("./PDF/paper/prior_post.png", width=9, height=3*4, units="in", res=300, pointsize=11)
{
  pat = "P09"
  chan = "NDUFB8"
  
  op = par(mfrow=c(4,2), mar=c(5,5,1,1))
  
  post = as.data.frame( fread( paste0("./Output/stan_sampler/", chan, "__", pat, "_POST.txt") ) )
  prior = as.data.frame( fread( paste0("./Output/stan_sampler/", chan, "__", pat, "_PRIOR.txt") ) )
  
  xlb = list(probdiff=bquote(pi), tau_norm=bquote(tau))
  names(xlb) = c("probdiff", "tau_norm")
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


# --- better comparison of classifs
install.packages("HDInterval")
library("HDInterval")

mdat = as.data.frame( fread( file.path("../dat_with_class_prepped.txt") ) )

bayes_post = list()
bayes_diff = list()
bayes_mean = list()
freq_diff = list()
nonZero = list()
for( pat in pts ){
  for( chan in channels ){
    root = paste0(chan, "__", pat)
    
    fn_post = file.path("Output", folder, paste0(root, "_POST.txt"))
    probdiff_post = as.data.frame( fread( fn_post ))[,"probdiff"]
    
    man_est = mean( unlist(mdat[ mdat$patient_id==pat & mdat$channel==chan, "classJBC"]) )
    
    diff_vec = probdiff_post - man_est
    
    fn_freq = file.path("Output", "frequentist_linReg", paste0(root, "_CLASSIF.txt"))
    freq_class = as.data.frame( fread( fn_freq ) )
    
    bayes_post[[root]] = probdiff_post
    bayes_diff[[root]] = diff_vec
    bayes_mean[[root]] = mean(diff_vec)
    freq_diff[[root]] = mean( freq_class[,1] ) - man_est
    diff_hdi = hdi(diff_vec, credMass = 0.95)
    nonZero[[root]] = (diff_hdi["lower"]>0.0) | (diff_hdi["upper"]<0.0)
  }
}

png("./PDF/paper/proportion_comparison_better.png", width=9, height=6, units="in", res=300, pointsize=11)
{
  op = par(mar=c(4,4,1,9), xpd=TRUE)
  
  colVec = c(alphaCol(0.01, name="lightseagreen"), 
             alphaCol(0.01, name="deeppink3"), 
             alphaCol(0.01, name="darkorange"))
  
  myAt = c(1:44)[ 1:44 %% 4 !=0 ]
  myTicks = 4* 1:11 - 2
  plot(NA, xlim=c(1,43), ylim=c(-0.1,1), xaxt="n", 
       ylab="Difference from manual classification", xlab="Patient")
  
  for(i in seq(0,10, by=2)) rect( xleft=i*4, xright=(i+1)*4, ybottom=par("usr")[3], ytop=par("usr")[4], border=NA, col=alphaBlack(0.1))
  lines(x=c(0,44), c(0,0), lty=3, col="black", lwd=3)
  
  stripchart(bayes_diff, vertical=TRUE, pch=20, cex=0.5, 
             method="jitter", main="", ylim=c(0,1), xlim=c(1,43),
             xlab="Patient", ylab="Deficient Proportion", 
             col=rep(colVec, length(bayes_diff)), 
             at=myAt, xaxt="n", add=TRUE)
  axis(side=1, at=myTicks, labels=pts)
  points(myAt, unlist(freq_diff), pch=24, cex=1.2,
         col="black",
         bg=rep(c(alphaCol(0.9, name="lightseagreen"), 
                   alphaCol(0.9, name="deeppink3"), 
                   alphaCol(0.9, name="darkorange")), length(freq_diff)))
  #points(x=myAt[unlist(nonZero)], y=rep(-0.1, sum(unlist(nonZero))), pch=8)
  legend(x=45, y=0.99, 
         # inset=c(0, 0), 
         legend=channels, 
         # title="Bayesian\ndistributions",
         #title.adj=0.1, 
         box.col=alphaBlack(0.0),
         text.width = 10,
         pch=20, 
         col=c(alphaCol(0.7, name="lightseagreen"), 
               alphaCol(0.7, name="deeppink3"), 
               alphaCol(0.7, name="darkorange")))
  text(x=45, y=1.0, labels=c("Bayesian distribution"), pos=4)
  rect(xleft=45, xright=55.7, 
       ybottom=0.82, ytop=1.025, 
       lwd=1, col=NA, border="black")
  
  legend(x=45, y=0.725, 
         inset=c(0.0, 0.0), 
         legend=channels, 
         # title = "Frequentist point\nestimate",
         #title.adj=0.1, 
         box.col=alphaBlack(0.0), 
         text.width = 10,
         pch=24, 
         col="black",
         pt.bg=c(alphaCol(0.7, name="lightseagreen"), 
               alphaCol(0.7, name="deeppink3"), 
               alphaCol(0.7, name="darkorange")) )
  text(x=45, y=0.75, labels=c("Frequentist point\nestimate"), pos=4)
  rect(xleft=45, xright=55.7, 
       ybottom=0.56, ytop=0.81, 
       lwd=1, col=NA, border="black")
  par(op)
}
dev.off()

# --- gamma comparison - MEAN ABSOLUTE DIFFERENCE
mad = list()
gamma_labs = c("000001", "00001", "0001", "001", "01", "10")
gamma_vals = c("0.00001", "0.0001", "0.001", "0.01", "0.1", "1.0")

for( gam in gamma_labs){
  fldr = paste0("stan_sampler_gamma", gam)
  diffs = NULL
  for( pat in pts ){
    for( chan in channels ){
      mclass = unlist( mdat[ mdat$channel==chan & mdat$patient_id==pat, "classJBC"] )
      root = paste0(chan, "__", pat)
      fn_post = file.path("Output", fldr, paste0(root, "_POST.txt"))
      
      post = as.matrix( fread( fn_post) )
      diffs = c(diffs, abs(post[,"probdiff"] - mean(mclass)))
    }
  }
  mad[[gam]] = mean(diffs)
}

png("./PDF/paper/gamma_comparison.png", width=6, height=4, units="in", res=300, pointsize=11)
{
  op = par(mfrow=c(1,1), mar=c(6,4,2,2))
  plot(1:length(gamma_vals), unlist(mad), ylim=c(0,0.1),
       ylab="Mean absolute difference", xlab=bquote(gamma), 
       xaxt='n',
       pch=20)
  axis(side=1, at=1:length(gamma_vals), labels=gamma_vals)
  
  par(op)
}
dev.off()

# --- COMPARING PI TO THE MANUAL CLASSIFICATION

mdat = fread( "../dat_with_class_prepped.txt" )

pis_diff = list()
nonZero = list()
for( gam in c("000001", "00001", "0001", "001", "01")){
  fldr = paste0("stan_sampler_gamma", gam)
  pis_diff[[gam]] = list()
  for( pat in pts ){
    for( chan in channels ){
      mclass = unlist( mdat[ mdat$channel==chan & mdat$patient_id==pat, "classJBC"] )
      root = paste0(chan, "__", pat)
      fn_post = file.path("Output", fldr, paste0(root, "_POST.txt"))
      
      post = as.matrix( fread( fn_post) )
      diff_vec = post[,"probdiff"] - mean(mclass)
      pis_diff[[gam]][[root]] = diff_vec
      diff_hdi = hdi(diff_vec, ci=0.99)
      nonZero[[gam]][[root]] = 0.0<(diff_hdi["lower"]) | 0.0>(diff_hdi["upper"])
    }
  }
}

png("./PDF/paper/gamma_comparison.png", width=9, height=6, units="in", res=300, pointsize=11)
{
  op = par(mfrow=c(2,1), mar=c(6,4,6,2), xpd=TRUE)
  gam_lbs = c("0.00001", "0.0001", "0.001", "0.01", "0.1", "1.0")
  names(gam_lbs) = c("000001", "00001", "0001", "001", "01", "10")
  
  colVec = c(alphaCol(0.05, name="lightseagreen"), 
             alphaCol(0.05, name="deeppink3"), 
             alphaCol(0.05, name="darkorange"))
  
  myAt = c(1:44)[ 1:44 %% 4 !=0 ]
  myTicks = 4* 1:11 - 2
  
  opp = par(mar=c(3,4,1,6), xpd=TRUE)
  {
    plot(NA, xlim=c(1,43), ylim=c(-0.1,0.2), 
         xaxt="n", 
         ylab="", xlab="", 
         main="")
    for(i in seq(0,10, by=2)) rect( xleft=i*4, xright=(i+1)*4, ybottom=par("usr")[3], ytop=par("usr")[4], border=NA, col=alphaBlack(0.1))
    
    stripchart(pis_diff[["001"]], vertical=TRUE, pch=20, cex=0.5, 
               method="jitter", 
               col=rep(colVec, length(pis_diff)), 
               at=myAt, xaxt='n',
               add=TRUE)
    text(46, 0.05, labels=bquote(gamma*"="*.(gam_lbs["001"])),
         srt=-90)
    axis(side=1, at=myTicks, labels=pts)
    lines(c(0,44), c(0,0), lty=3, col="black", lwd=2)
    
    points(x=myAt[unlist(nonZero[["001"]])], y=rep(-0.1, sum(unlist(nonZero[["001"]]))), 
           pch=8)
    exact_xlim = par("usr")[1:2]
    exact_ylim = par("usr")[3:4]
    
    legend(exact_xlim[2] + 1.0, exact_ylim[2], 
           legend=channels, 
           title = "Channel",
           pch=20,  cex=0.7, 
           xpd=TRUE,
           col=c(alphaCol(0.7, name="lightseagreen"), 
                 alphaCol(0.7, name="deeppink3"), 
                 alphaCol(0.7, name="darkorange")) )
  }
  par(opp)
  
  opp = par(mar=c(4,4,0,6), xpd=TRUE)
  {
    plot(NA, xlim=c(1,43), ylim=c(-0.1,0.2), xaxt="n", 
         ylab="", xlab="Patient",
         main="")
    for(i in seq(0,10, by=2)) rect( xleft=i*4, xright=(i+1)*4, ybottom=par("usr")[3], ytop=par("usr")[4], border=NA, col=alphaBlack(0.1))
    
    stripchart(pis_diff[["00001"]], vertical=TRUE, pch=20, cex=0.5, 
               method="jitter",
               col=rep(colVec, length(pis_diff)), 
               at=myAt, xaxt="n", add=TRUE)
    axis(side=1, at=myTicks, labels=pts)
    text(46, 0.05, labels=bquote(gamma*"="*.(gam_lbs["00001"])),
         srt=-90)
    lines(c(0,44), c(0,0), lty=3, col="black", lwd=2)
    points(x=myAt[unlist(nonZero[["00001"]])], y=rep(-0.1, sum(unlist(nonZero[["00001"]]))), 
           pch=8)
  }
  par(opp)
  
  mtext("Difference from manual classification", side=2, line=-1, outer=TRUE, adj=0.6, padj=0.5)
  
  par(op)
}
dev.off()

# --- comparing priors against manual classifications

folders = c("stan_sampler_wide", "stan_sampler_narrow")
names(folders) = c("wide", "narrow")

mdat = fread( "../dat_with_class_prepped.txt" )

pis_diff = list()
for( name in names(folders) ){
  pis_diff[[name]] = list()
  for( pat in pts ){
    for( chan in channels ){
      mclass = unlist( mdat[ mdat$channel==chan & mdat$patient_id==pat, "classJBC"] )
      
      root = paste0(chan, "__", pat)
      
      fn_post = file.path("Output", folders[name], paste0(root, "_POST.txt"))
      
      post = as.matrix( fread( fn_post) )
      diff_vec = post[,"probdiff"] - mean(mclass)
      pis_diff[[name]][[root]] = diff_vec
    }
  }
}

png("./PDF/paper/prior_comparison_with_manual.png", width=9, height=6, units="in", res=300, pointsize=11)
{
  op = par(mfrow=c(2,1), mar=c(6,4,6,2), xpd=TRUE)
  
  colVec = c(alphaCol(0.05, name="lightseagreen"), 
             alphaCol(0.05, name="deeppink3"), 
             alphaCol(0.05, name="darkorange"))
  
  myAt = c(1:44)[ 1:44 %% 4 !=0 ]
  myTicks = 4* 1:11 - 2
  
  opp = par(mar=c(3,4,1,6), xpd=TRUE)
  {
    plot(NA, xlim=c(1,43), ylim=c(-0.1,0.2), 
         xaxt="n", 
         ylab="", xlab="", 
         main="")
    for(i in seq(0,10, by=2)) rect( xleft=i*4, xright=(i+1)*4, ybottom=par("usr")[3], ytop=par("usr")[4], border=NA, col=alphaBlack(0.1))
    
    stripchart(pis_diff[["wide"]], vertical=TRUE, pch=20, cex=0.5, 
               method="jitter", 
               col=rep(colVec, length(pis_diff)), 
               at=myAt, xaxt='n',
               add=TRUE)
    text(46, 0.05, labels="Wide prior", srt=-90)
    axis(side=1, at=myTicks, labels=pts)
    lines(c(0,44), c(0,0), lty=3, col="black", lwd=2)
    exact_xlim = par("usr")[1:2]
    exact_ylim = par("usr")[3:4]
    
    legend(exact_xlim[2] + 1.0, exact_ylim[2], 
           legend=channels, 
           title = "Channel",
           pch=20,  cex=0.7, 
           xpd=TRUE,
           col=c(alphaCol(0.7, name="lightseagreen"), 
                 alphaCol(0.7, name="deeppink3"), 
                 alphaCol(0.7, name="darkorange")) )
  }
  par(opp)
  
  opp = par(mar=c(4,4,0,6), xpd=TRUE)
  {
    plot(NA, xlim=c(1,43), ylim=c(-0.1,0.2), xaxt="n", 
         ylab="", xlab="Patient",
         main="")
    for(i in seq(0,10, by=2)) rect( xleft=i*4, xright=(i+1)*4, ybottom=par("usr")[3], ytop=par("usr")[4], border=NA, col=alphaBlack(0.1))
    
    stripchart(pis_diff[["narrow"]], vertical=TRUE, pch=20, cex=0.5, 
               method="jitter",
               col=rep(colVec, length(pis_diff)), 
               at=myAt, xaxt="n", add=TRUE)
    axis(side=1, at=myTicks, labels=pts)
    text(46, 0.05, labels="Narrow prior", srt=-90)
    lines(c(0,44), c(0,0), lty=3, col="black", lwd=2)
  }
  par(opp)
  
  mtext("Difference from manual classification", side=2, line=-1, outer=TRUE, adj=0.6, padj=0.5)
  
  par(op)
}
dev.off()

# --- prior comparison against OG posterior

folders = c("stan_sampler_wide", "stan_sampler_narrow")
names(folders) = c("wide", "narrow")

pis_diff = list()
pis_diff[["wide"]] = list()
pis_diff[["narrow"]] = list()
for(pat in pts) {
  for (chan in channels) {
    root = paste0(chan, "__", pat)
    
    fn_post_og = file.path("Output", "stan_sampler", paste0(root, "_POST.txt"))
    pi_post_og = as.matrix(fread(fn_post_og))[, "probdiff"]
    

    for (pp in names(folders)) {
      fn_post_two = file.path("Output", folders[pp], paste0(root, "_POST.txt"))
      pi_post_two = as.matrix(fread(fn_post_two))[, "probdiff"]
      pis_diff[[pp]][[root]] = pi_post_og - pi_post_two
    }
  }
}

png("./PDF/paper/prior_comparison.png", width=9, height=6, units="in", res=300, pointsize=11)
{
  op = par(mfrow=c(2,1), mar=c(6,4,6,2), xpd=TRUE)
  
  colVec = c(alphaCol(0.01, name="lightseagreen"), 
             alphaCol(0.01, name="deeppink3"), 
             alphaCol(0.01, name="darkorange"))
  
  myAt = c(1:44)[ 1:44 %% 4 !=0 ]
  myTicks = 4* 1:11 - 2
  
  opp = par(mar=c(3,4,1,6), xpd=TRUE)
  {
    plot(NA, xlim=c(1,43), ylim=c(-0.2,0.2), 
         xaxt="n", 
         ylab="", xlab="", 
         main="")
    for(i in seq(0,10, by=2)) rect( xleft=i*4, xright=(i+1)*4, ybottom=par("usr")[3], ytop=par("usr")[4], border=NA, col=alphaBlack(0.1))
    
    stripchart(pis_diff[["wide"]], vertical=TRUE, pch=20, cex=0.5, 
               method="jitter", 
               col=rep(colVec, length(pis_diff)), 
               at=myAt, xaxt='n',
               add=TRUE)
    text(46, -0.01, labels="Wide prior", srt=-90)
    axis(side=1, at=myTicks, labels=pts)
    lines(c(0,44), c(0,0), lty=3, col="black", lwd=2)
    exact_xlim = par("usr")[1:2]
    exact_ylim = par("usr")[3:4]
    
    legend(exact_xlim[2] + 1.0, exact_ylim[2], 
           legend=channels, 
           title = "Channel",
           pch=20,  cex=0.7, 
           xpd=TRUE,
           col=c(alphaCol(0.7, name="lightseagreen"), 
                 alphaCol(0.7, name="deeppink3"), 
                 alphaCol(0.7, name="darkorange")) )
  }
  par(opp)
  
  opp = par(mar=c(4,4,0,6), xpd=TRUE)
  {
    plot(NA, xlim=c(1,43), ylim=c(-0.2,0.2), xaxt="n", 
         ylab="", xlab="Patient",
         main="")
    for(i in seq(0,10, by=2)) rect( xleft=i*4, xright=(i+1)*4, ybottom=par("usr")[3], ytop=par("usr")[4], border=NA, col=alphaBlack(0.1))
    
    stripchart(pis_diff[["narrow"]], vertical=TRUE, pch=20, cex=0.5, 
               method="jitter",
               col=rep(colVec, length(pis_diff)), 
               at=myAt, xaxt="n", add=TRUE)
    axis(side=1, at=myTicks, labels=pts)
    text(46, -0.01, labels="Narrow prior", srt=-90)
    lines(c(0,44), c(0,0), lty=3, col="black", lwd=2)
  }
  par(opp)
  
  mtext("Difference distribution", side=2, line=-1, outer=TRUE, adj=0.6, padj=0.5)
  
  par(op)
}
dev.off()

# --- ALL POSTERIOR PREDICTIVIONS - supplementary

png("./PDF/paper/allPost_classifs_1.png", width=9, height=3*4, units="in", res=300, pointsize=11)
{
  op = par(mfrow=c(4,3), mar=c(4,4,2,1), xpd=FALSE)
  
  for( pat in pts[1:4] ){
    for( chan in channels ){
      root = paste0(chan, "__", pat)
      fn_post = file.path("Output", "stan_sampler", paste0(root, "_POST.txt"))
      fn_postpred = file.path("Output", "stan_sampler", paste0(root, "_POSTPRED.txt"))
      fn_classif = file.path("Output", "stan_sampler", paste0(root, "_CLASSIF.txt"))
      
      class_mats = as.matrix( fread(fn_classif) )
      class = colMeans( class_mats ) 
      
      postpred = as.matrix( fread(fn_postpred) )
      post = as.matrix( fread(fn_post) )
      pi_mean = round(mean( post[,"probdiff"] ), 3)
      
      dataMats = getData_mats(data, channels=c(mitochan, chan), 
                              pts=pat, ctrlID=ctrlIDs)
      
      classif_plot(dataMats, 
                   classifs = class,
                   postpred = postpred, 
                   xlab=paste0("log(", mitochan, ")"), 
                   ylab=paste0("log(", chan, ")"),
                   main=bquote("E("*pi*"|Y)="*.(pi_mean)))
      if( chan == "NDUFB8" ){
        text(par("usr")[1], par("usr")[4], labels=pat, adj=c(-0.2,1.3), cex=2)
      }
    }
  }
  par(op)
}
dev.off()

png("./PDF/paper/allPost_classifs_2.png", width=9, height=3*4, units="in", res=300, pointsize=11)
{
  op = par(mfrow=c(4,3), mar=c(4,4,2,1), xpd=FALSE)
  
  for( pat in pts[5:8] ){
    for( chan in channels ){
      root = paste0(chan, "__", pat)
      fn_post = file.path("Output", "stan_sampler", paste0(root, "_POST.txt"))
      fn_postpred = file.path("Output", "stan_sampler", paste0(root, "_POSTPRED.txt"))
      fn_classif = file.path("Output", "stan_sampler", paste0(root, "_CLASSIF.txt"))
      
      class_mats = as.matrix( fread(fn_classif) )
      class = colMeans( class_mats ) 
      
      postpred = as.matrix( fread(fn_postpred) )
      post = as.matrix( fread(fn_post) )
      pi_mean = round(mean( post[,"probdiff"] ), 3)
      
      dataMats = getData_mats(data, channels=c(mitochan, chan), 
                              pts=pat, ctrlID=ctrlIDs)
      
      classif_plot(dataMats, 
                   classifs = class,
                   postpred = postpred, 
                   xlab=paste0("log(", mitochan, ")"), 
                   ylab=paste0("log(", chan, ")"),
                   main=bquote("E("*pi*"|Y)="*.(pi_mean)))
      if( chan == "NDUFB8" ){
        text(par("usr")[1], par("usr")[4], labels=pat, adj=c(-0.2,1.3), cex=2)
      }
    }
  }
  par(op)
}
dev.off()

png("./PDF/paper/allPost_classifs_3.png", width=9, height=3*3, units="in", res=300, pointsize=11)
{
  op = par(mfrow=c(3,3), mar=c(4,4,2,1), xpd=FALSE)
  
  for( pat in pts[9:11] ){
    for( chan in channels ){
      root = paste0(chan, "__", pat)
      fn_post = file.path("Output", "stan_sampler", paste0(root, "_POST.txt"))
      fn_postpred = file.path("Output", "stan_sampler", paste0(root, "_POSTPRED.txt"))
      fn_classif = file.path("Output", "stan_sampler", paste0(root, "_CLASSIF.txt"))
      
      class_mats = as.matrix( fread(fn_classif) )
      class = colMeans( class_mats ) 
      
      postpred = as.matrix( fread(fn_postpred) )
      post = as.matrix( fread(fn_post) )
      pi_mean = round(mean( post[,"probdiff"] ), 3)
      
      dataMats = getData_mats(data, channels=c(mitochan, chan), 
                              pts=pat, ctrlID=ctrlIDs)
      
      classif_plot(dataMats, 
                   classifs = class,
                   postpred = postpred, 
                   xlab=paste0("log(", mitochan, ")"), 
                   ylab=paste0("log(", chan, ")"),
                   main=bquote("E("*pi*"|Y)="*.(pi_mean)))
      if( chan == "NDUFB8" ){
        text(par("usr")[1], par("usr")[4], labels=pat, adj=c(-0.2,1.3), cex=2)
      }
    }
  }
  par(op)
}
dev.off()

# --- pi posterior - supplementary
bayes_post = list()
for( pat in pts ){
  for( chan in channels ){
    root = paste0(chan, "__", pat)
    fn_post = file.path("Output", folder, paste0(root, "_POST.txt"))
    probdiff_post = as.data.frame( fread( fn_post ))[,"probdiff"]
    bayes_post[[root]] = probdiff_post
  }
}

png("./PDF/paper/piPost.png", width=9, height=6, units="in", res=300, pointsize=11)
{
  op = par(mar=c(4,4,1,1), xpd=FALSE)
  
  colVec = c(alphaCol(0.01, name="lightseagreen"), 
             alphaCol(0.01, name="deeppink3"), 
             alphaCol(0.01, name="darkorange"))
  
  myAt = c(1:44)[ 1:44 %% 4 !=0 ]
  myTicks = 4* 1:11 - 2
  plot(NA, xlim=c(1,43), ylim=c(0.0,1), xaxt="n", 
       ylab="Deficiency proportion", xlab="Patient")
  
  for(i in seq(0,10, by=2)) rect( xleft=i*4, xright=(i+1)*4, ybottom=par("usr")[3], ytop=par("usr")[4], border=NA, col=alphaBlack(0.1))
  
  stripchart(bayes_post, vertical=TRUE, pch=20, cex=0.5, 
             method="jitter", main="", ylim=c(0,1), xlim=c(1,43),
             xlab="Patient", ylab="Deficient Proportion", 
             col=rep(colVec, length(bayes_post)), 
             at=myAt, xaxt="n", add=TRUE)
  axis(side=1, at=myTicks, labels=pts)
  legend("topleft", 
         bg = "transparent",
         legend=channels, 
         box.col=alphaBlack(0.0),
         text.width = 10,
         pch=20, 
         col=c(alphaCol(0.7, name="lightseagreen"), 
               alphaCol(0.7, name="deeppink3"), 
               alphaCol(0.7, name="darkorange")))
  text(x=45, y=1.0, labels=c("Bayesian distribution"), pos=4)
  rect(xleft=45, xright=55.7, 
       ybottom=0.82, ytop=1.025, 
       lwd=1, col=NA, border="black")
  par(op)
}
dev.off()





