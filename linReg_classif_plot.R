
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
channels = c("MTCO1")

######################
### PRIORS AND POSTERIORS
######################
post = list()
postpred = list()
prior = list()
priorpred = list()
classif = list()
{
  post$ctrl = list()
  post$pat = list()
  postpred$ctrl = list()
  postpred$pat = list()
  classif$pat = list()
  classif$ctrl = list()
  prior$ctrl = list()
  prior$pat = list()
  priorpred$ctrl = list()
  priorpred$pat = list()
  
  for(chan in channels){
    outroot_ctrl = paste(chan, "CONTROL", sep="_")
    
    fpc_prior = paste0("./Output/linReg_classif/", outroot_ctrl, "_PRIOR.txt" )
    fpc_post = paste0("./Output/linReg_classif/", outroot_ctrl, "_POST.txt" )      
    fpc_classif = paste0("./Output/linReg_classif/", outroot_ctrl, "_CLASSIF.txt" )
    fpc_priorpred = paste0("./Output/linReg_classif/", outroot_ctrl, "_PRIORPRED.txt")
    fpc_postpred = paste0("./Output/linReg_classif/", outroot_ctrl, "_POSTPRED.txt")
    # if(file.exists(fpc_prior)) prior$ctrl[[outroot_ctrl]] = read.table(fpc_prior, header=TRUE, stringsAsFactors = FALSE)
    # if(file.exists(fpc_post)) post$ctrl[[outroot_ctrl]] = read.table(fpc_post, header=TRUE, stringsAsFactors = FALSE)
    # if(file.exists(fpc_classif)) classif$ctrl[[outroot_ctrl]] = read.table(fpc_classif, header=TRUE, stringsAsFactors = FALSE)
    if(file.exists(fpc_priorpred)) priorpred$ctrl[[outroot_ctrl]] = read.table(fpc_priorpred, header=TRUE, stringsAsFactor=FALSE)
    if(file.exists(fpc_postpred)) postpred$ctrl[[outroot_ctrl]] = read.table(fpc_postpred, header=TRUE, stringsAsFactor=FALSE)
    # 
    for(pat in pts){
      outroot_pat = paste(chan, pat, sep="_")
      fpp_prior = paste0("./Output/linReg_classif/", outroot_pat, "_PRIOR.txt" )
      fpp_post = paste0("./Output/linReg_classif/", outroot_pat, "_POST.txt" )
      fpp_classif = paste0("./Output/linReg_classif/", outroot_pat, "_CLASSIF.txt" )
      fpp_priorpred = paste0("./Output/linReg_classif/", outroot_pat, "_PRIORPRED.txt")
      fpp_postpred = paste0("./Output/linReg_classif/", outroot_pat, "_POSTPRED.txt")
      # if(file.exists(fpp_post)) post$pat[[outroot_pat]] = read.table(fpp_post, header=TRUE, stringsAsFactors=FALSE)
      # if(file.exists(fpp_classif)) classif$pat[[outroot_pat]] = read.table(fpp_classif, header=TRUE, stringsAsFactors=FALSE)
      # if(file.exists(fpp_prior)) prior$pat[[outroot_pat]] = read.table(fpp_prior, header=TRUE, stringsAsFactors = FALSE)
      if(file.exists(fpp_priorpred)) priorpred$pat[[outroot_ctrl]] = read.table(fpp_priorpred, header=TRUE, stringsAsFactor=FALSE)
      if(file.exists(fpp_postpred)) postpred$pat[[outroot_ctrl]] = read.table(fpp_postpred, header=TRUE, stringsAsFactor=FALSE)

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
  priorpost_comp = function(data, priorpred, postpred, pat_data = NULL, 
                            classif, lims, Xsyn, 
                            chan, mitochan="VDAC1", title){
     
    xlims = lims[[1]]
    ylims = lims[[2]]
    
    par(mfrow=c(1,2))
    plot(data, pch=20, col=myGrey(0.1), xlim=xlims, ylim=ylims,
         xlab=paste0("log(",mitochan,")"), ylab=paste0("log(",chan,")"), 
         main="Control Prior")
    if(!is.null(pat_data)) points(pat_data, pch=20, col=myYellow(0.2))
    priorpred_qnts = rowQuantiles(priorpred)
    lines(Xsyn, priorpred_qnts[,c(1,3)], lty=2, col="green", lwd=2)
    lines(Xsyn, priorpred_qnts[,2], lty=1, col="green", lwd=3)
    
    plot(ctrl_data, pch=20, col=myGrey(0.1), xlim=xlims, ylim=ylims,
         xlab=paste0("log(",mitochan,")"), ylab=paste0("log(",chan,")"), 
         main="Control Prior")
    if(!is.null(pat_data)) points(pat_data, pch=20, col=classcols(classif))
    postpred_qnts = rowQuantiles(post_pat)
    lines(Xsyn, postpred_pat_qnts[,c(1,3)], lty=2, col="blue", lwd=2)
    lines(Xsyn, postpred_pat_qnts[,2], lty=1, col="blue", lwd=3)
    
    title(main=paste(title), line=-1, outer=TRUE, cex.main=1.6)
  }
  
  pi_post = function(post_pat, chan, title){
    # do something
  }
}


##
## THE GOODS
##
pdf("./PDF/linReg_classif/MCMC.pdf", width=13, height=8)
{
  
  for(chan in channels){
    froot_ctrl = paste(chan, "CONTROL", sep="_")
    MCMCplot(post$ctrl[[froot_ctrl]], prior$ctrl[[froot_ctrl]], froot_ctrl)
    for(pat in pts){
      froot_pat = paste(chan, pat, sep="_")
      MCMCplot(post$pat[[chanpat]], prior$pat[[chanpat]], froot_pat)
    }
  }
}
dev.off()


pdf("./PDF/linReg_classif/model_post.pdf", width=13, height=8)
{
  for(chan in c("MTCO1")){
    froot = paste(chan, "CONTROL", sep="_")
    data = getData_mats(chan=chan)
    
    lims = list(c(0,6), c(0,6))
    Xsyn = seq(0.9*min(data$ctrl[,1]), 1.1*max(data$ctrl[,1]), length.out=1000)

    priorpost_comp(data=ctrl_data$Yctrl, priorpred=priorpred$ctrl[[froot]], 
                   postpred=postpred$ctrl[[froot]], classif=classif$ctrl[[froot]], 
                   lims=lims, Xsyn=Xsyn, 
                   chan=chan, title=paste(chan, "CONTROL"))
    for(pat in pts){
      froot = paste(chan, pat, sep="_")
      pat_data = getData_mats(chan=chan, pat=pat)
      priorpost_comp(data=ctrl_data$Yctrl, pat_data=pat_data$Ypat, 
                     post=postpred, prior=priorpred, 
                     classif=classif$pat[[froot]], 
                     title=paste(chan, pat))
    }
  }
}
dev.off()




