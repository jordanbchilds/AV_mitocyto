library(data.table)
library(MASS)

##
## DATA COLLEECTION
##

getData_mats = function(chan, mitochan="VDAC", 
                        pts=NULL, ctrl_only=FALSE, 
                        data_transform=log_transform){

  dat = read.csv("../Data_prepped.csv", header=TRUE)
  dat_chan = dat[, c("ID", "patient_id", mitochan, chan)]
  Xctrl = dat_chan[grepl("C", dat_chan$patient_id), mitochan]
  Yctrl = dat_chan[grepl("C", dat_chan$patient_id), chan]
  ctrl_mat = cbind( Xctrl, Yctrl )
  pts_all = unique( dat_chan$patient_id )
  
  if(!ctrl_only){
    if(is.null(pts)){
      Ypts = list()
      for(pat in pts_all[grepl("P", pts_all)]){
        Xpat = dat_chan[dat_chan$patient_id == pat, mitochan]
        Ypat = dat_chan[dat_chan$patient_id == pat, chan]
        XY_pat = cbind(Xpat, Ypat)
        Ypts[[pat]] = XY_pat
      }
    } else {
      Ypts = list()
      for(pat in pts){
        Xpat = dat_chan[dat_chan$patient_id == pat, mitochan]
        Ypat = dat_chan[dat_chan$patient_id == pat, chan]
        XY_pat = cbind(Xpat, Ypat)
        Ypts[[pat]] = XY_pat
      }
    }
    pat_mat = list2matrix(Ypts)
  } else { pat_mat=NULL }
  
  if(!is.null(data_transform)) return( data_transform(ctrl_mat, pat_mat) )
  
  if(ctrl_only) return( ctrl_mat )
  list(ctrl=ctrl_mat, pts=pat_mat)
}

##
## COLOURS
##

pal = palette()
palette(c(pal, "darkblue", "darkorange", "deeppink"))

myBlack = function(alpha) rgb(0,0,0, alpha)
myDarkGrey = function(alpha) rgb(169/255,169/255,159/255, alpha)
myGrey = function(alpha) rgb(66/255,66/255,66/255, alpha)
myBlue = function(alpha) rgb(0,0,128/255, alpha)
myRed = function(alpha) rgb(1,0,0, alpha)

myGreen = function(alpha) rgb(0,100/255,0, alpha)
myYellow = function(alpha) rgb(225/255,200/255,50/255, alpha)
myPink = function(alpha) rgb(255/255,62/255,150/255, alpha)
myPurple = function(alpha) rgb(160/255, 32/255, 240/255, alpha)

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

colQuantiles = function(x, probs=0.5){
  quants = matrix(NA, nrow=ncol(x), ncol=length(probs))
  for(i in 1:ncol(x)){
    quants[i,] = quantile(x[,i], probs)
  }
  colnames(quants) = probs
  return(quants)
}

##
## GEN FUNCTIONS
##

vec_rep = function(x, n, byRow=TRUE){
  out = matrix(x, nrow=n, ncol=length(x), byrow=TRUE)
  if(byRow) return( out )
  else return( t(out) )
}

list2matrix = function(X, rowBind=TRUE){
  # X must list of matrices or a matrix
  if(is.matrix(X)) return(X)
  if(!is.list(X)) stop("X must be a list of matrices or a matrix")
  if(length(X)==1) return(X[[1]])
  
  n = length(X)
  out = X[[1]]
  if(rowBind){
    for(i in 2:n){ out = rbind(out, X[[i]]) }
    return(out)
  } else {
    for(i in 2:n){ out = cbind(out, X[[i]]) }
    return(out)
  }
}

log_transform = function(ctrl_mat, pat_mat=NULL){
  if( !is.null(pat_mat) ) return( list(ctrl=log(ctrl_mat), pts=log(pat_mat)) )
  log( ctrl_mat )
}

##
## READ & SAVERS AOUTPUT
##

output_saver = function(output, folder, outroot){
  split = strsplit(outroot, split="_")[[1]]
  chan = split[1]
  pat = split[2]
  
  for(ctrl_pat in names(output)){
    out_ctrlpat = ifelse(ctrl_pat=="ctrl", "CONTROL", pat)
    for(out_type in names(output[[ctrl_pat]])){
      filename = paste(chan, out_ctrlpat, toupper(out_type), sep="_")
      write.table(output[[ctrl_pat]][[out_type]], paste0(file.path("Output", folder, filename), ".txt"),
                  row.names=FALSE, col.names=TRUE)
    }
  }
}

output_reader = function(folder, chan, pat, out_type){
  outroot = paste(chan, pat, sep="_")
  
  fp = file.path("Output", folder, paste0(outroot, "_", out_type, ".txt"))
  if(file.exists(fp)) return( read.table(fp, header=TRUE, stringsAsFactors=FALSE) )
  else stop("file does not exist")
}

##
## PLOTTERS
##

  
classif_plot = function(ctrl_data, pat_data, classifs_pat, chan, mitochan, pat){
    
    xrange = range(c(ctrl_data[,1], pat_data[,1]))
    yrange = range(c(ctrl_data[,2], pat_data[,2]))
    
    plot(ctrl_data, pch=20, col=myDarkGrey(0.2), xlim=xrange, ylim=yrange, 
         xlab=paste(mitochan), ylab=paste(chan), 
         main=paste(pat), cex.lab=1.2, cex.main=1.4)
    points(pat_data, pch=20, col=classcols(classifs_pat))
  }
  
MCMCplot = function(folder, chan, pat, title="", lag=20){
  
    post = output_reader(folder, chan, pat, out_type="POST")
    prior = output_reader(folder, chan, pat,  out_type="PRIOR")
    
    col.names = colnames(post)
    n.chains = length(post)
    par(mfrow=c(2,3))
    for(param in col.names){
      post_vec = post[[param]]
      plot(ts(post_vec), xlab="Iteration", ylab=paste(param), 
           main="", cex.lab=1.2)
      if(sum(post_vec==post_vec[1])!=length(post_vec)){
        acf(post[[param]], xlab="lag index", ylab="ACF", main="",
            cex.lab=1.2, lag=lag)
      } else {
        plot(NA, type='n', xlim=c(0,lag), ylim=c(0,1), 
             xlab="lag index", ylab="ACF", main="")
      }
      plot(density(post[[param]]), lwd=2, col="blue", xlab=paste(param), ylab="Density",
           main="")
      if(param %in% colnames(prior)) lines(density(prior[[param]]), lwd=2, col="green")
      
      title(main=title, line=-1, outer=TRUE, cex.main=1.6)
    }
}

priorpost = function(ctrl_data, pat_data=NULL, priorpred, postpred,
                     classif=NULL, 
                     chan, mitochan="VDAC1", title="", xlims=NULL, ylims=NULL){

  N_syn = 1e4
  Xsyn = seq(0, max(c(ctrl_data[,1], pat_data[,1]))*1.5, length.out=N_syn) 
  
  op = par(mfrow=c(1,2))
  plot(ctrl_data, pch=20, cex=0.7, col=myGrey(0.1),
        xlab=paste0("log(",mitochan,")"), ylab=paste0("log(",chan,")"), 
       main="Prior Predictive", xlim=xlims, ylim=ylims)
  if(!is.null(pat_data)) points(pat_data, pch=20, cex=1.2, col=myYellow(0.2))
  lines(Xsyn, priorpred[,1], lty=2, col=myGreen(0.6), lwd=3)
  lines(Xsyn, priorpred[,2], lty=1, col=myGreen(0.6), lwd=4)
  lines(Xsyn, priorpred[,3], lty=2, col=myGreen(0.6), lwd=3)
  
  plot(ctrl_data, pch=20, col=myGrey(0.1),
       xlab=paste0("log(",mitochan,")"), ylab=paste0("log(",chan,")"), 
       main="Posterior Predictive", xlim=xlims, ylim=ylims)
  if(!is.null(pat_data)) points(pat_data, pch=20, cex=1.2, col=classcols(classif))
  lines(Xsyn, postpred[,1], lty=2, col=myPink(0.6), lwd=3)
  lines(Xsyn, postpred[,2], lty=1, col=myPink(0.6), lwd=4)
  lines(Xsyn, postpred[,3], lty=2, col=myPink(0.6), lwd=3)
  title(main=title, line=-2, outer=TRUE)
  par(op)
}


pipost_plotter = function(chan, folder, pts, alpha=0.05){
  npat = length(pts)
  pis = list()
  for(pat in pts){
    pis[[pat]] = 1 - output_reader(folder, chan, pat, out_type="POST")[,"probdiff"]
  }

  title = paste( chan )

  stripchart(pis, pch=20, method="jitter", vertical=TRUE, 
             col=rgb(t(col2rgb(palette()[colvector_gen(pts)]))/255, alpha=alpha), 
             group.names=pts,
             at=1:length(pts),
             main=title, ylim=c(0,1), ylab="Deficiency Proportion", 
             xlab="Patient Sample")
}

pipost_plotter_v2 = function(channels, pts, folder, alpha=0.01){
  npat = length(pts)
  pis = list()
  for(pat in pts){
    for(chan in channels){
      pis[[paste(pat, chan, sep="_")]] = 1 - output_reader(folder, chan, pat, out_type="POST")[,"probdiff"]
    }
  }
  
  stripchart(pis, pch=20, method="jitter", vertical=TRUE, 
             col=rgb(t(col2rgb(rep(palette()[1:length(channels) + 1], length(pts))))/255, alpha=alpha),
             at=1:(length(channels)*length(pts)), 
             xaxt="n", ylim=c(0,1),
             main="", ylab="Deficiency Proportion", 
             xlab="Patient" )
  abline( v=(0:length(pts)-1)*3+3.5, lwd=4, lty="dotted", col=myGrey(0.2))
  axis(1, at=(1:length(pts)-1)*3+2, labels=pts)
  legend("topright", legend=channels, pch=20, col=rep(palette()[1:length(channels) + 1]),
         cex=1.5, bty="o", bg="white", title="channels")
  
  
}













































































