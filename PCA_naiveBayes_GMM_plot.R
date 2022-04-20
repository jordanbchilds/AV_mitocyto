library("MASS")
source("../BootStrapping/parseData.R", local = TRUE)

dir.create(file.path("./PDF"), showWarnings = FALSE)
dir.create(file.path("./PDF/PCA_naiveBayes_GMM"), showWarnings = FALSE)

getData_mats = function(fulldat="Data_prepped.csv", mitochan="VDAC", chan, 
                        pts=NULL, ctrl_only=FALSE, PCA=TRUE){
  data_raw = read.csv(fulldat, header=TRUE)
  
  pts_raw = unique(data_raw$patient_id)
  
  data = data_raw
  data[data[,"patient_id"] %in% pts_raw[grep("C0", pts_raw)],"patient_id"] = "control"
  pts_all = unique(data_raw$patient_id)
  
  ctrl_data = data[data$patient_id=="control", ]
  Xctrl = log(ctrl_data[[mitochan]])
  Yctrl = log(ctrl_data[[chan]])
  ctrl_mat = cbind( Xctrl, Yctrl )
  
  if(!ctrl_only){
    if(is.null(pts)){
      pat_mat = matrix(NA, nrow=1, ncol=2)
      Npats = vector("numeric")
      for(pat in pts_all[grepl("P", pts_all)]){
        pat_data = data[data$patient_id==pat,]
        Xpat = log(pat_data[[mitochan]])
        Ypat = log(pat_data[[chan]])
        XY_pat = cbind(Xpat, Ypat)
        pat_mat = rbind(pat_mat, XY_pat)
        Npats[pat] = nrow(XY_pat)
      }
    } else {
      pat_mat = matrix(NA, nrow=1, ncol=2)
      Npats = vector("numeric")
      for(pat in pts){
        pat_data = data[data$patient_id==pat,]
        Xpat = log(pat_data[[mitochan]])
        Ypat = log(pat_data[[chan]])
        XY_pat = cbind(Xpat, Ypat)
        pat_mat = rbind(pat_mat, XY_pat)
        Npats[pat] = nrow(XY_pat)
      }
    }
  }
  
  pat_mat = pat_mat[-1,]
  
  if(!PCA){
    if(ctrl_only) return(ctrl_mat)
    else return(list(ctrl=ctrl_mat, pat=pat_mat, Npats=Npats))
  } else {
    # PCA
    pca = prcomp(ctrl_mat, center=FALSE, scale=FALSE)
    pca_mean = colMeans(pca$x)
    ctrlmat_cen = sweep(pca$x, 2, pca_mean)
    if(ctrl_only) return(ctrlmat_cen)
    else {
      patmat_cen = sweep(pat_mat%*%pca$rotation, 2, pca_mean)  
      return(list(ctrl=ctrlmat_cen, pat=patmat_cen, Npats=Npats))
    }
  }
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

output_reader = function(out_type, chans=channels){
  out = list()
  for(chan in chans){
    for(pat in pts){
      outroot = paste(chan, pat, sep="_")
      fp = paste0("./OUTPUT/PCA_naiveBayes_GMM/_", outroot, "_", out_type, ".txt")
      if(file.exists(fp)) out[[outroot]] = read.table(fp, header=TRUE, stringsAsFactors=FALSE)
    }
  }
  return(out)
}

post = output_reader(out_type="POST")
postpred = output_reader(out_type="POSTPRED")
classif = output_reader(out_type="CLASS")
prior = read.table("./OUTPUT/PCA_naiveBayes_GMM/allData_PRIOR.txt", header=TRUE, stringsAsFactors=FALSE)
priorpred = read.table("./OUTPUT/PCA_naiveBayes_GMM/allData_PRIORPRED.txt", header=TRUE, stringsAsFactors=FALSE)

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

# plotting functions
{
percentiles = function(xdat, ydat, probs=c(0.975, 0.5, 0.025)){
  dens = kde2d(xdat, ydat, n=200); ## estimate the z counts
  dx = diff(dens$x[1:2])
  dy = diff(dens$y[1:2])
  sz = sort(dens$z)
  c1 = cumsum(sz) * dx * dy
  levs = sapply(probs, function(x) {
    approx(c1, sz, xout = 1 - x)$y
  })
  return( list(dens=dens, levels=levs, probs=probs))
}

priorpost = function(data, prior, posterior, classifs, ctrl=NULL,
                     title){
  # output: plots the prior and posterior regression lines and data
  op = par(mfrow=c(1,2), mar = c(5.5,5.5,3,3))
  if( is.null(ctrl) ){
    x.lim = range(data[,1]) + c(-1,1)
    y.lim = range(data[,2]) + c(-1,1)
    
    plot(data[,1], data[,2], col=myDarkGrey, pch=20, cex.axis=1.5,
         xlab=paste("log(",mitochan,")"), ylab=paste("log(",chan,")"), main='Prior Predictive',
         xlim=x.lim, ylim=y.lim)
    contours = percentiles(prior[,"predOne[1]"], prior[,"predOne[2]"])
    contour( contours$dens, levels=contours$levels, labels=contours$probs, add=TRUE, lwd=2 )
    
    plot(data[,1], data[,2], col=myDarkGrey, pch=20, cex.axis=1.5,
         xlab=paste("log(",mitochan,")"), ylab=paste("log(",chan,")"), main='Posterior Predictive',
         xlim=x.lim, ylim=y.lim)
    contours = percentiles(posterior[,"predOne[1]"], posterior[,"predOne[2]"])
    contour( contours$dens, levels=contours$levels, labels=contours$probs, add=TRUE, lwd=2 )
    
    title(main=title, line = -1, outer = TRUE)
  } else {
    x.lim = range( data[,1], ctrl[,1] ) + c(-1,1)
    y.lim = range( data[,2], ctrl[,2] ) + c(-1,1)
    
    plot(ctrl[,1], ctrl[,2], col=myDarkGrey, pch=20, cex.axis=1.5,
         xlab=paste("log(",mitochan,")"), ylab=paste("log(",chan,")"), main='Prior Predictive',
         xlim=x.lim, ylim=y.lim)
    points( data[,1], data[,2], col=myYellow, pch=20)
    contours = percentiles(prior[,"predOne[1]"], prior[,"predOne[2]"])
    contour( contours$dens, levels=contours$levels, labels=contours$probs, add=TRUE, lwd=2 )
    
    plot(ctrl[,1], ctrl[,2], col=myDarkGrey, pch=20, cex.axis=1.5,
         xlab=paste("log(",mitochan,")"), ylab=paste("log(",chan,")"), main='Posterior Predictive',
         xlim=x.lim, ylim=y.lim)
    
    points( data[,1], data[,2], col=classcols(classifs), pch=20)
    contours = percentiles(posterior[,"predOne[1]"], posterior[,"predOne[2]"])
    contour( contours$dens, levels=contours$levels, labels=contours$probs, add=TRUE, lwd=2 )
    title(main=title, line = -1, outer = TRUE)
  }
  
  par(op)
} 

priorpost_marginals = function(prior, posterior, title){
  op = par(mfrow=c(1,2), mar = c(5.5,5.5,3,3))
  par(mfrow=c(2,2))
  ## mu_1
  # prior
  contour( kde2d(prior[,'mu[1,1]'], prior[,'mu[2,1]'], n=100), cex.lab=2, cex.axis=1.5,
           xlab=expression(mu[11]), ylab=expression(mu[12]), nlevels=5,
           main=expression(mu[1]~'Prior Density') )
  # posterior
  contour( kde2d(posterior[,'mu[1,1]'], posterior[,'mu[2,1]'], n=100), cex.lab=2, cex.axis=1.5,
           xlab=expression(mu[11]), ylab=expression(mu[12]), nlevels=5,
           main=expression(mu[1]~'Posterior Density') )
  ## mu_2
  # prior
  contour( kde2d(prior[,'mu[1,2]'], prior[,'mu[2,2]'], n=100), cex.lab=2, cex.axis=1.5,
           xlab=expression(mu[21]), ylab=expression(mu[22]), nlevels=5,
           main=expression(mu[2]~'Prior Density') )
  # posterior
  contour( kde2d(posterior[,'mu[1,2]'], posterior[,'mu[2,2]'], n=100), cex.lab=2, cex.axis=1.5,
           xlab=expression(mu[21]), ylab=expression(mu[22]), nlevels=5,
           main=expression(mu[2]~'Posterior Density') )
  title(main=title, line = -1, outer = TRUE)
  
  
  par(mfrow=c(1,3))
  ## tau_1
  # prior
  plot( density(prior[,"tau[1,1,1]"]), cex.lab=2, cex.axis=1.5, lwd=2, 
        xlab=expression(tau[11]), ylab="", main=expression(T[1]~"Prior Density"))
  lines( density(posterior[,"tau[1,1,1]"]), col=myGreen, lwd=2)
  
  plot( density(prior[,"tau[1,2,1]"]), cex.lab=2, cex.axis=1.5, lwd=2, 
        xlab=expression(tau[12]), ylab="", main=expression(T[1]~"Prior Density"))
  lines( density(posterior[,"tau[1,1,1]"]), lwd=2, col=myGreen)
  
  plot( density(prior[,"tau[2,2,1]"]), cex.lab=2, cex.axis=1.5, lwd=2, 
        xlab=expression(tau[22]), ylab="", main=expression(T[1]~"Prior Density"))
  lines( density(posterior[,"tau[2,2,1]"]), lwd=2, col=myGreen)
  
  ## tau_1
  # posterior
  plot( density(prior[,"tau[1,1,2]"]), cex.lab=2, cex.axis=1.5, lwd=2, 
        xlab=expression(tau[11]), ylab="", main=expression(T[2]~"Prior Density"))
  lines( density(posterior[,"tau[1,1,2]"]), col=myGreen, lwd=2)
  plot( density(prior[,"tau[1,2,2]"]), cex.lab=2, cex.axis=1.5, lwd=2, 
        xlab=expression(tau[12]), ylab="", main=expression(T[2]~"Prior Density"))
  lines( density(posterior[,"tau[1,2,2]"]), col=myGreen, lwd=2)
  plot( density(prior[,"tau[2,2,2]"]), cex.lab=2, cex.axis=1.5, lwd=2, 
        xlab=expression(tau[22]), ylab="", main=expression(T[2]~"Prior Density"))
  lines( density(posterior[,"tau[2,2,2]"]), col=myGreen, lwd=2)
  title(main=title, line = -1, outer = TRUE)
  
  par(mfrow=c(1,2))
  plot( density(posterior[,"probdiff_pat"]), cex.lab=2, cex.axis=1.5, xlim=c(0,1),          
        xlab="probdiff_pat", ylab="", lwd=2, col="red", main="probdiff_pat Density")
  lines( density(prior[,"probdiff_pat"]), lwd=2, col="green")
  title(main=title, line = -1, outer = TRUE)
  
  par(op)
}

component_densities = function(ctrl_data, pat_data, 
                               post_one, post_two, prior_one, prior_two, 
                               classifs, title, chan, mitochan ){
  
  xrange = range(c(ctrl_data[,1], pat_data[,1])) + c(-1,1)
  yrange = range(c(ctrl_data[,2], pat_data[,2])) + c(-1,1)
  
  par(mfrow=c(2,2))
  plot(ctrl_data, pch=20, col=myDarkGrey(0.6),
       xlab=paste("log(",mitochan,")"), ylab=paste("log(",chan,")"),
       main="Component One", xlim=xrange, ylim=yrange)
  points(pat_data, pch=20, col=myYellow(0.6))
  prior_one_dens = percentiles(prior_one[,1], prior_one[,2])
  contour(prior_one_dens$dens, levels=prior_one_dens$levels, labels=prior_one_dens$probs,
          col=myBlue(0.4), lwd=2, add=TRUE)
  
  plot(ctrl_data, pch=20, col=myDarkGrey(0.6),
       xlab=paste("log(",mitochan,")"), ylab=paste("log(",chan,")"),
       main="Component Two", xlim=xrange, ylim=yrange)
  points( pat_data, pch=20, col=myYellow(0.6))
  prior_two_dens = percentiles(prior_two[,1], prior_two[,2])
  contour(prior_two_dens$dens, levels=prior_two_dens$levels, labels=prior_two_dens$probs, 
          col=myRed(0.4), lwd=2, add=TRUE)
  
  plot(ctrl_data, pch=20, col=myDarkGrey(0.6),
       xlab=paste("log(",mitochan,")"), ylab=paste("log(",chan,")"),
       main="Component One", xlim=xrange, ylim=yrange)
  points( pat_data, pch=20, col=classcols(classifs))
  post_one_dens = percentiles(post_one[,1], post_one[,2])
  contour(post_one_dens$dens, levels=post_one_dens$levels, labels=post_one_dens$probs,
          col=myBlue(0.4), lwd=2, add=TRUE)
  
  plot(ctrl_data, pch=20, col=myDarkGrey(0.6),
       xlab=paste("log(",mitochan,")"), ylab=paste("log(",chan,")"),
       main="Component Two", xlim=xrange, ylim=yrange)
  points( pat_data, pch=20, col=classcols(classifs))
  post_two_dens = percentiles(post_two[,1], post_two[,2])
  contour(post_two_dens$dens, levels=post_two_dens$levels, labels=post_two_dens$probs, 
          col=myRed(0.9), lwd=2, add=TRUE)
  
  title(main=title, line = -2, outer = TRUE)
}

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
  for(param in col.names[col.names!="lp__"]){
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
}

pdf("./PDF/PCA_naiveBayes_GMM/model_post.pdf", width=12, height=10)
{
  for(outroot in names(post)){
    chan = strsplit(outroot, split="_")[[1]][1]
    pat = strsplit(outroot, split="_")[[1]][2]
    
    mitochan = "VDAC1"
    
    data_mats = getData_mats(chan=chan, pts=pat)
    ctrl_data = data_mats$ctrl
    pat_data = data_mats$pat
    
    post_cp = post[[outroot]]
    post_one = as.matrix(post_cp[,c("comp.1.1.","comp.2.1.")])
    post_two = as.matrix(post_cp[,c("comp.1.2.","comp.2.2.")])
    prior_one = as.matrix(priorpred[,c("comp.1.1.","comp.2.1.")])
    prior_two = as.matrix(priorpred[,c("comp.1.2.","comp.2.2.")])
    
    # par(mfrow=c(2,2))
    # plot(post_one, col=myBlack(0.05), main="post one")
    # plot(post_two, col=myBlack(0.05), main="post two")
    # plot(prior_one, col=myBlack(0.05), main="prior one")
    # plot(prior_two, col=myBlack(0.05), main="prior two")
    # 
    component_densities(ctrl_data, pat_data,
                        post_one, post_two, prior_one, prior_two,
                        classif[[outroot]][,1], title=paste(chan, pat), chan, mitochan)
    
  }
}
dev.off()

pdf("./PDF/PCA_naiveBayes_GMM/classifs.pdf", width=13, height=8)
{
  par(mfrow=c(1,2))
  for(outroot in names(post)){
    chan = strsplit(outroot, split="_")[[1]][1]
    pat = strsplit(outroot, split="_")[[1]][2]
    mitochan = "VDAC1"
    
    data_mats = getData_mats(chan=chan, pts=pat, PCA=FALSE)
    ctrl_data = data_mats$ctrl
    pat_data = data_mats$pat
    
    classif_pat = classif[[outroot]][,"x"]
    
    classif_plot(ctrl_data=ctrl_data, pat_data=pat_data, classifs_pat=classif_pat, 
                 chan, mitochan, pat)
  }
}
dev.off()

pdf("./PDF/PCA_naiveBayes_GMM/pipost_pat.pdf", width=13, height=8)
{
  par(mfrow=c(1,1))
  for(pat in pts){
    pipost_pat = list()
    for(chan in channels){
      pipost_pat[[chan]] = post[[paste(chan,pat, sep="_")]][,"probctrl"]
    }
    stripchart(pipost_pat, method="jitter", pch=".", col=myBlack(0.1), vertical=TRUE,
                ylab="like ctrl Proportion", xlab="Protein", 
                main=paste0(pat), cex.lab=1.4, cex.main=1.4, 
                ylim=c(0,1))
  }
}
dev.off()

pdf("./PDF/PCA_naiveBayes_GMM/pipost_chan.pdf", width=13, height=8)
{
  par(mfrow=c(1,1))
  for(chan in channels){
    pipost_chan = list()
    for(pat in pts){
      pipost_chan[[pat]] = post[[paste(chan,pat, sep="_")]][,"probctrl"]
    }
    stripchart(pipost_chan, method="jitter", pch=".", col=myBlack(0.1), vertical=TRUE,
               ylab="like ctrl Proportion", xlab="Patient", 
               main=paste0(chan), cex.lab=1.4, cex.main=1.4, 
               ylim=c(0,1))
  }
}
dev.off()

pdf("./PDF/PCA_naiveBayes_GMM/MCMC.pdf", width=13, height=8)
{
  for(outroot in names(post)){
    MCMCplot(post[[outroot]], prior, outroot)
  }
}
dev.off()









