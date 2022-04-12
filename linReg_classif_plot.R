
dir.create(file.path("./PDF"), showWarnings = FALSE)
dir.create(file.path("./PDF/linReg_classif"), showWarnings = FALSE)

getData_mats = function(fulldat="Data_prepped.csv", mitochan="VDAC", chan, pat){
  data_raw = read.csv(fulldat, header=TRUE)
  
  pts_raw = unique(data_raw$patient_id)
  
  data = data_raw
  data[data[,"patient_id"] %in% pts_raw[grep("C0", pts_raw)],"patient_id"] = "control"
  
  ctrl_data = data[data$patient_id=="control", ]
  Xctrl = log(ctrl_data[[mitochan]])
  Yctrl = log(ctrl_data[[chan]])
  XY_ctrl = cbind( Xctrl, Yctrl )
  
  pat_data = data[data$patient_id==pat,]
  Xpat = log(pat_data[[mitochan]])
  Ypat = log(pat_data[[chan]])
  XY_pat = cbind(Xpat, Ypat)
  
  return(list(Yctrl=XY_ctrl, Ypat=XY_pat))
}

fp_data = "Data_prepped.csv"
mc_raw = read.csv(fp_data, header=TRUE)
pts_raw = unique(mc_raw$patient_id)
mc_data = mc_raw
mc_data[mc_raw[,"patient_id"] %in% pts_raw[grep("C0", pts_raw)],"patient_id"] = "control"
pts_all = unique(mc_data$patient_id)
pts  = pts_all[pts_all!="control"]
mitochan = "VDAC"
channels = c("MTCO1")


######################
### PRIORS AND POSTERIORS
######################
post = list()
classif = list()
prior = list()
{
  post$ctrl = list()
  post$pat = list()
  classif$pat = list()
  classif$ctrl = list()
  prior$ctrl = list()
  prior$pat = list()
  
  for(chan in channels){
    outroot_ctrl = paste(chan, "CONTROL", sep="_")
    
    fp_ctrlprior = paste0("./Output/linReg_classif/", outroot_ctrl, "_PRIOR.txt" )
    fp_ctrlpost = paste0("./Output/linReg_classif/", outroot_ctrl, "_POST.txt" )      
    fp_ctrlclassif = paste0("./Output/linReg_classif/", outroot_ctrl, "_CLASSIF.txt" )
    if(file.exists(fp_ctrlprior)) prior$ctrl[[outroot_ctrl]] = read.table(fp_ctrlprior, header=TRUE, stringsAsFactors = FALSE)
    if(file.exists(fp_ctrlpost)) post$ctrl[[outroot_ctrl]] = read.table(fp_ctrlpost, header=TRUE, stringsAsFactors = FALSE)
    if(file.exists(fp_ctrlclassif)) classif$ctrl[[outroot_ctrl]] = read.table(fp_ctrlclassif, header=TRUE, stringsAsFactors = FALSE)
      
    for(pat in pts){
      outroot_pat = paste(chan, pat, sep="_")
      
      fp_prior = paste0("./Output/linReg_classif/", outroot_pat, "_PRIOR.txt" )
      fp_post = paste0("./Output/linReg_classif/", outroot_ctrl, "_POST.txt" )
      fp_classif = paste0("./Output/linReg_classif/", outroot_ctrl, "_CLASSIF.txt" )
      if(file.exists(fp_post)) post[[outroot_pat]] = read.table(fp_post, header=TRUE, stringsAsFactors=FALSE)
      if(file.exists(fp_classif)) classif[[outroot_pat]] = read.table(fp_classif, header=TRUE, stringsAsFactors=FALSE)
      if(file.exists(fp_prior)) prior[[outroot_pat]] = read.table(fp_prior, header=TRUE, stringsAsFactors = FALSE)
    } #  patients
  } # channels
}

##
## COLOURS
##
myBlack = function(alpha) rgb(0,0,0, alpha)
myDarkGrey = function(alpha) rgb(169/255,169/255,159/255, alpha)
myGrey = function(alpha) rgb(66/255,66/255,66/255, alpha)
myBlue = function(alpha) rgb(0,0,128/255, alpha)
myRed = function(alpha) rgb(1,0,0, alpha)
myGreen = function(alpha) rgb(0,100/255,0, alpha)
myYellow = function(alpha) rgb(225/255,200/255,50/255, alpha)

cramp = colorRamp(c(myRed(0.2),myBlue(0.2)), alpha=TRUE)
# rgb(...) specifies a colour using standard RGB, where 1 is the maxColorValue
# 0.25 determines how transparent the colour is, 1 being opaque 
# cramp is a function which generates colours on a scale between two specifies colours

classcols = function(classif){
  # A function using the cramp specified colours to generate rgb colour names
  # input: a number in [0,1]
  # output: rgb colour name
  rgbvals = cramp(classif)/255.0
  return(rgb(rgbvals[,1],rgbvals[,2],rgbvals[,3], alpha=rgbvals[,4]))
}

##
## plotting functions
##
{
  classif_plot = function(ctrl_data, pat_data, classifs_pat, chan, mitochan, pat){
    
    xrange = range(c(ctrl_data[,1], pat_data[,1]))
    yrange = range(c(ctrl_data[,2], pat_data[,2]))
    
    plot(ctrl_data, pch=20, col=myDarkGrey(0.2), xlim=xrange, ylim=yrange, 
         xlab=paste(mitochan), ylab=paste(chan), 
         main=paste(pat), cex.lab=1.2, cex.main=1.4)
    points(pat_data, pch=20, col=classcols(classifs_pat))
  }
  
  MCMCplot = function(post, prior, title){
    col.names = colnames(post)
    n.chains = length(post)
    
    par(mfrow=c(3,3), mar = c(5.5,5.5,4,4))
    for(param in col.names){
      plot(ts(post[[param]]), xlab="Iteration", ylab=paste(param), 
           main="", cex.lab=1.2)
      acf(post[[param]], xlab="lag index", ylab="ACF", main="",
          cex.lab=1.2)
      plot(density(post[[param]]), lwd=2, col="blue", xlab=paste(param), ylab="Density",
           main="")
      if(param %in% colnames(prior)) lines(density(prior[[param]]), lwd=2, col="green")
      
      title(main=title, line=-2, outer=TRUE, cex.main=1.6)
    }
  }
  
  rowQuantiles = function(x, probs=0.5){
    quants = matrix(NA, nrow=nrow(x), ncol=length(probs))
    for(i in 1:nrow(x)){
      quants[i,] = quantile(x[i,], probs)
    }
    quants
  }
  priorpost_comp = function(ctrl_data, pat_data, prior_ctrl, post_ctrl,
                            prior_pat, post_pat, classif_ctrl, classif_pat, 
                            chan, mitochan="VDAC1", title){
    
    xlims = range(c(ctrl_data[,1], pat_data[,1]))
    ylims = range(c(ctrl_data[,2], pat_data[,2]))
    Xsyn = seq(0.9*min(ctrl_data[,1]), 1.1*max(ctrl_data[,1]), length.out=1000)
    
    par(mfrow=c(2,2))
    
    plot(ctrl_data, pch=20, col=myGrey(0.2), xlim=xlims, ylim=ylims,
         xlab=paste0("log(",mitochan,")"), ylab=paste0("log(",chan,")"), 
         main="Control Prior")
    priorpred_ctrl_qnts = rowQuantiles(prior_ctrl[,grepl("Ysyn", colnames(prior_ctrl))])
    lines(Xsyn, priorpred_ctrl_qnts[,c(1,3)], lty=2, col="green", lwd=2)
    lines(Xsyn, priorpred_ctrl_qnts[,2], lty=1, col="green", lwd=3)
    
    plot(ctrl_data, pch=20, col=classcols(classif_ctrl), xlim=xlims, ylim=ylims,
         xlab=paste0("log(",mitochan,")"), ylab=paste0("log(",chan,")"), 
         main="Control Prior")
    postpred_ctrl_qnts = rowQuantiles(post_ctrl[,grepl("Ysyn", colnames(post_ctrl))])
    lines(Xsyn, postpred_ctrl_qnts[,c(1,3)], lty=2, col="blue", lwd=2)
    lines(Xsyn, postpred_ctrl_qnts[,2], lty=1, col="blue", lwd=3)
    
    plot(ctrl_data, pch=20, col=myGrey(0.1), xlim=xlims, ylim=ylims,
         xlab=paste0("log(",mitochan,")"), ylab=paste0("log(",chan,")"), 
         main="Control Prior")
    points(pat_data, pch=20, col=myYellow(0.2))
    priorpred_pat_qnts = rowQuantiles(prior_pat[,grepl("Ysyn", colnames(prior_pat))])
    lines(Xsyn, priorpred_pat_qnts[,c(1,3)], lty=2, col="green", lwd=2)
    lines(Xsyn, priorpred_pat_qnts[,2], lty=1, col="green", lwd=3)
    
    plot(ctrl_data, pch=20, col=myGrey(0.1), xlim=xlims, ylim=ylims,
         xlab=paste0("log(",mitochan,")"), ylab=paste0("log(",chan,")"), 
         main="Control Prior")
    points(pat_data, pch=20, col=classcols(classif_pat))
    postpred_pat_qnts = rowQuantiles(post_pat[,grepl("Ysyn", colnames(post_pat))])
    lines(Xsyn, postpred_pat_qnts[,c(1,3)], lty=2, col="blue", lwd=2)
    lines(Xsyn, postpred_pat_qnts[,2], lty=1, col="blue", lwd=3)
  }
  
  pi_post = function(post_pat, chan, title){
    # do something
  }
}


##
## THE GOODS
##
pdf("./PDF/PCA_naiveBayes_GMM/MCMC.pdf", width=13, height=8)
{
  for(chanpat in names(chanpat_post)){
    MCMCplot(chanpat_post[[chanpat]], prior, chanpat)
  }
}
dev.off()


