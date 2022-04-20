library(MASS)

dir.create(file.path("./PDF"), showWarnings = FALSE)
dir.create(file.path("./PDF/linReg_classif"), showWarnings = FALSE)

getData_mats = function(fulldat="Data_prepped.csv", mitochan="VDAC", chan, 
                        pts=NULL, ctrl_only=FALSE){
  data_raw = read.csv(fulldat, header=TRUE)
  
  pts_raw = unique(data_raw$patient_id)
  
  data = data_raw
  data[data[,"patient_id"] %in% pts_raw[grep("C0", pts_raw)],"patient_id"] = "control"
  pts_all = unique(data_raw$patient_id)
  
  ctrl_data = data[data$patient_id=="control", ]
  Xctrl = log(ctrl_data[[mitochan]])
  Yctrl = log(ctrl_data[[chan]])
  XY_ctrl = cbind( Xctrl, Yctrl )
  
  if(!ctrl_only){
    if(is.null(pts)){
      Ypat_all = matrix(NA, nrow=1, ncol=2)
      Npats = vector("numeric")
      for(pat in pts_all[grepl("P", pts_all)]){
        pat_data = data[data$patient_id==pat,]
        Xpat = log(pat_data[[mitochan]])
        Ypat = log(pat_data[[chan]])
        XY_pat = cbind(Xpat, Ypat)
        Ypat_all = rbind(Ypat_all, XY_pat)
        Npats[pat] = nrow(XY_pat)
      }
    } else {
      Ypat_all = matrix(NA, nrow=1, ncol=2)
      Npats = vector("numeric")
      for(pat in pts){
        pat_data = data[data$patient_id==pat,]
        Xpat = log(pat_data[[mitochan]])
        Ypat = log(pat_data[[chan]])
        XY_pat = cbind(Xpat, Ypat)
        Ypat_all = rbind(Ypat_all, XY_pat)
        Npats[pat] = nrow(XY_pat)
      }
    }
  }
  
  if(ctrl_only) return(XY_ctrl)
  return(list(ctrl=XY_ctrl, pat=Ypat_all[-1,], Npats=Npats))
}

fp_data = "Data_prepped.csv"
mc_raw = read.csv(fp_data, header=TRUE)
pts_raw = unique(mc_raw$patient_id)
mc_data = mc_raw
mc_data[mc_raw[,"patient_id"] %in% pts_raw[grep("C0", pts_raw)],"patient_id"] = "control"
pts_all = unique(mc_data$patient_id)
pts  = pts_all[pts_all!="control"]
mitochan = "VDAC"
channels = c("MTCO1", "NDUFB8", "CYB")


######################
### PRIORS AND POSTERIORS
######################

output_reader = function(fulldats="", out_type, channels=channels){
  out = list()
  for(fulldat in fulldats){
    ctrl_pats = c("CONTROL", pts)
    for(chan in channels){
      for(pat in ctrl_pats){
        outroot = paste(froot, chan, pat, sep="_")
        fp = paste0("./Output/linReg_classif/", outroot, "_", out_type, ".txt")
        if(file.exists(fp)) out[[outroot]] = read.table(fp, header=TRUE, stringsAsFactors=FALSE)
      }
    }
  }
  return(out)
}

post = output_reader(out_type="POST")
postpred = output_reader(out_type="POSTPRED")
classif = output_reader(out_type="CLASS")
prior = output_reader(out_type="PRIOR")
priorpred = output_reader(out_type="PRIORPRED")


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
myPink = function(alpha) rgb(255/255,62/255,150/255, alpha)

cramp = colorRamp(c(myBlue(0.2),myRed(0.2)), alpha=TRUE)
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
    
    par(mfrow=c(2,3), mar = c(5.5,5.5,4,4))
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
  priorpost_comp = function(ctrl_data, pat_data=NULL, priorpred, postpred,
                            classif=NULL, 
                            chan, mitochan="VDAC1", title){
    
    xlims = range(c(ctrl_data[,1], pat_data[,1]))
    ylims = range(c(ctrl_data[,2], pat_data[,2]))
    Xsyn = seq(min(ctrl_data[,1])-1, max(ctrl_data[,1])+1, length.out=1000)
    
    par(mfrow=c(1,2))
    plot(ctrl_data, pch=20, cex=0.7, col=myGrey(0.1), xlim=xlims, ylim=ylims,
         xlab=paste0("log(",mitochan,")"), ylab=paste0("log(",chan,")"), 
         main="Prior Predictive")
    if(!is.null(pat_data)) points(pat_data, pch=20, cex=1.2, col=myYellow(0.2))
    lines(Xsyn, priorpred[,1], lty=2, col=myGreen(0.6), lwd=3)
    lines(Xsyn, priorpred[,2], lty=1, col=myGreen(0.6), lwd=4)
    lines(Xsyn, priorpred[,3], lty=2, col=myGreen(0.6), lwd=3)
    
    plot(ctrl_data, pch=20, col=myGrey(0.1), xlim=xlims, ylim=ylims,
         xlab=paste0("log(",mitochan,")"), ylab=paste0("log(",chan,")"), 
         main="Posterior Predictive")
    if(!is.null(pat_data)) points(pat_data, pch=20, cex=1.2, col=classcols(classif))
    lines(Xsyn, postpred[,1], lty=2, col=myPink(0.6), lwd=3)
    lines(Xsyn, postpred[,2], lty=1, col=myPink(0.6), lwd=4)
    lines(Xsyn, postpred[,3], lty=2, col=myPink(0.6), lwd=3)
    
  }
  
  pi_post = function(pipost_list, title){
    # do something
  }
}


##
## THE GOODS
##
pdf("./PDF/linReg_classif/MCMC.pdf", width=13, height=8)
{
  for(chan in channels){
    outroot_ctrl = paste(chan, "CONTROL", sep="_")
    MCMCplot(post[[outroot_ctrl]], prior[[outroot_ctrl]], 
             title=paste(chan, "CONTROL"))
    for(pat in pts){
      outroot_pat = paste(chan, pat, sep="_")
      MCMCplot(post[[outroot_pat]], prior[[outroot_pat]], 
               title=paste(chan, pat))
    }
  }
}
dev.off()

pdf("./PDF/linReg_classif/model_post.pdf", width=13, height=8)
{
  for(chan in channels){
    outroot_ctrl = paste(chan, "CONTROL", sep="_")
    ctrl_data = getData_mats(chan=chan, ctrl_only=TRUE)
    priorpost_comp(ctrl_data=ctrl_data, 
                   priorpred=priorpred[[outroot_ctrl]], postpred=postpred[[outroot_ctrl]],
                   chan=chan, mitochan="VDAC1", title=paste(chan, 'CONTROL'))
    
    for(pat in pts){
      outroot_pat = paste(chan, pat, sep="_")
      pat_data = getData_mats(chan=chan, pts=pat)$pat
      
      priorpost_comp(ctrl_data=ctrl_data, pat_data=pat_data,
                     classif=classif[[outroot_pat]][,1],
                     priorpred=priorpred[[outroot_pat]], postpred=postpred[[outroot_pat]],
                     chan=chan, mitochan="VDAC1", title=paste(chan, pat))
    }
  }
}
dev.off()

pipost_list = list()
piprior_list = list()
for(chan in channels){
  temp_post = list()
  temp_prior = list()
  for(pat in pts){
    temp_post[[pat]] = 1 - post[[paste(chan, pat, sep="_")]][,"probdiff"]
    temp_prior[[pat]] = 1 - prior[[paste(chan, pat, sep="_")]][,"probdiff"]
  }
  pipost_list[[chan]] = temp_post
  piprior_list[[chan]] = temp_prior
}

pdf("./PDF/linReg_classif/pi_post.pdf", width=13, height=8)
{
  for(chan in channels){
    stripchart(pipost_list[[chan]], vertical=TRUE, pch=20, cex=0.7,
               col=myPink(0.02), method="jitter", jitter=0.1, ylim=c(0,1), 
               main=paste(chan), ylab="Like Control Proportion", 
               group.names=pts)
    # stripchart(pipost_list[[chan]], vertical=TRUE, pch=20, cex=0.5, col=myPink(0.035), 
    #            method="jitter", jitter=0.2, add=TRUE, at=c(1:length(pts)))
  }
}
dev.off()




















