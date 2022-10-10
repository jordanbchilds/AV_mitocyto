library(data.table)
library(MASS)

##
## DATA SPECIFIC FUNCTIONS 
##

getData = function(fname,cord,mitochan="VDAC1",updatechans=TRUE,correctnpc=FALSE){
  dat = fread(fname,sep="\t",stringsAsFactors=FALSE,header=TRUE)
  colnames(dat) = tolower(colnames(dat))
  
  if(updatechans) {
    dat$channel = gsub("GRIM19","NDUFA13",dat$channel)
    dat$channel = gsub("Ch1","Laminin",dat$channel)
    dat$channel = gsub("Ch2","MTCO1",dat$channel)
    dat$channel = gsub("Ch3","VDAC1",dat$channel)
    dat$channel = gsub("Ch4","SDHA",dat$channel)
    dat$ch = substring(dat$channel,regexpr("\\_[^\\_]*$", dat$channel)+1,nchar(dat$channel))
  }else{
    dat$ch = as.character(dat$channel)
  }
  #print("Available channels:")
  #print(unique(dat$ch))
  dat = dat[dat$ch%in%cord,]
  dat$fn = gsub(" NPC","",dat$filename)
  dat$fn = gsub(" OXPHOS","",dat$fn)
  dat$pch = paste(dat$fn,dat$ch)
  
  if (correctnpc){
    npc = dat[grepl("NPC",dat$filename),]
    oxphos = dat[grepl("OXPHOS",dat$filename),]
    agg = aggregate(npc$value,by=list(npc$pch),mean,data=dat)
    lu = agg$x
    names(lu) = agg$Group.1
    
    oxphos$filename = gsub(" OXPHOS","",oxphos$filename)
    oxphos$value = pmax(1.0,oxphos$value - lu[oxphos$pch])
    oxphos$filename = oxphos$fn
    oxphos$fn = NULL
    oxphos$pch = NULL
    dat = oxphos
  }
  
  dat$type = "Mean intensity"
  dat$type[grepl("LOG_",dat$channel)] = "Log mean intensity"
  dat$type[grepl("MED_",dat$channel)] = "Median intensity"
  dat$type[grepl("R_",dat$channel)] = "Ratio mean intensity (VDAC1)"
  dat$type[grepl("R_MED_",dat$channel)] = "Ratio median intensity (VDAC1)"
  dat$type[grepl("R_LOG_",dat$channel)] = "Ratio log mean intensity (VDAC1)"
  dat$type[grepl("Z_",dat$channel)] = "z-score"
  dat$outlier_diff = "NODIFF"
  dat$regression_diff = "NODIFF"
  dat$z_diff = "NODIFF"
  dat$z = 0
  
  dat$chstr = dat$ch
  transform = log
  dat_r = dat[dat$type=="Mean intensity",]
  dat_r$type = paste("r (",mitochan,")",sep="")
  dat_theta = dat[dat$type=="Mean intensity",]
  dat_theta$type = paste("theta (",mitochan,")",sep="")
  
  for(pid in unique(dat$patrep_id)){
    for(ch in unique(dat$ch)){
      dt = dat[(dat$patrep_id==pid)&(dat$type=="Mean intensity"),]
      
      isch = as.character(dt$ch)==ch
      ismito = as.character(dt$ch)==mitochan
      prot = dt[isch,]
      mito = dt[ismito,]
      
      x = mito$value
      y = prot$value
      dat_r$value[(dat_r$patrep_id==pid)&(as.character(dat_r$ch)==ch)] = sqrt(x^2+y^2)
      dat_r$channel[(dat_r$patrep_id==pid)&(as.character(dat_r$ch)==ch)] = paste("RADIUS",ch,sep="_")
      dat_theta$value[(dat_theta$patrep_id==pid)&(as.character(dat_theta$ch)==ch)] = 360*atan(y/x)/(2*pi)
      dat_theta$channel[(dat_theta$patrep_id==pid)&(as.character(dat_theta$ch)==ch)] = paste("THETA",ch,sep="_")
    }
  }
  dat=rbind(dat,dat_r,dat_theta)
  dat = data.frame(dat, stringsAsFactors=FALSE)
  return(dat)
}

getData_mats = function(fulldat, mitochan="VDAC1", chan, 
                        pts=NULL, ctrl_only=FALSE, 
                        data_transform=log_transform){
  froot = gsub(".RAW.txt","",fulldat)
  correctnpc = TRUE
  updatechans = FALSE
  
  data = getData(paste0("../",gsub(".RAW", ".RAW_ren", fulldat)), c(chan, mitochan),
                 mitochan = mitochan, updatechans = updatechans, correctnpc = correctnpc)
  
  data$fn = gsub("_.0", "", data$filename)
  data$pch = paste(data$fn, data$ch, sep="_")
  data$fn = gsub("_R1","",data$fn)
  data$fn = gsub("_R2","",data$fn)
  
  sbj = sort(unique(data$fn))
  crl = grep("C._H", sbj, value = TRUE)
  pts_all = grep("P._", sbj, value=TRUE)
  
  ctrl_data = data[(data$fn%in%crl)&(data$type=="Mean intensity"), ]
  Xctrl = ctrl_data$value[ctrl_data$channel==mitochan]
  Yctrl = ctrl_data$value[ctrl_data$channel==chan]
  ctrl_mat = cbind( Xctrl, Yctrl )
  
  if(!ctrl_only){
    if(is.null(pts)){
      Ypts = list()
      for(pat in pts_all[grepl("P", pts_all)]){
        pat_data = data[(data$fn==pat)&(data$type=="Mean intensity"),]
        Xpat = pat_data$value[pat_data$channel==mitochan]
        Ypat = pat_data$value[pat_data$channel==chan]
        XY_pat = cbind(Xpat, Ypat)
        Ypts[[pat]] = XY_pat
      }
    } else {
      Ypts = list()
      for(pat in pts){
        pat_data = data[(data$fn==pat)&(data$type=="Mean intensity"),]
        Xpat = pat_data$value[pat_data$channel==mitochan]
        Ypat = pat_data$value[pat_data$channel==chan]
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

getData_chanpats = function(fulldat){
  froot = gsub(".RAW.txt","",fulldat)
  mchans = c("Ch1","Ch2","Ch3","Ch4","Area")
  if(grepl(".TCF.",fulldat)){# attention confocal - swap channels
    chans=c("LAMA1","VDAC1","MTCO1","NDUFB8","Area") # swapped MTCO1 and VDAC1 for confocal
  }else{
    chans=c("LAMA1","MTCO1","VDAC1","NDUFB8","Area") # standard for CD7
  }
  names(chans) = mchans
  
  dat = read.delim(file.path("../../BootStrapping",fulldat),stringsAsFactors=FALSE)
  dat$Channel = chans[dat$Channel]
  write.table(dat,gsub(".RAW", ".RAW_ren", fulldat),row.names=FALSE,quote=FALSE,sep="\t")
  
  cord = c("NDUFB8","MTCO1","VDAC1")
  chlabs = c("CI","CIV","OMM")
  names(chlabs) = cord
  mitochan = "VDAC1"
  correctnpc = TRUE
  updatechans = FALSE
  dat = getData(gsub(".RAW", ".RAW_ren", fulldat), cord,
                mitochan = mitochan,updatechans = updatechans, correctnpc = correctnpc)
  
  dat$fn = gsub("_.0", "", dat$filename)
  dat$pch = paste(dat$fn,dat$ch,sep="_")
  # Merge different regions of the same section
  dat$fn = gsub("_R1","",dat$fn)
  dat$fn = gsub("_R2","",dat$fn)
  
  # grabbing and seperating ctrls and patients
  sbj = sort(unique(dat$fn))
  crl = grep("C._H", sbj, value = TRUE)
  pts = grep("P.", sbj, value = TRUE)
  
  return(list(channels=c("NDUFB8", "MTCO1"), patients=pts, controls=crl) )
}



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

##
## GENERAL FUNCTIONS FOR bayes GMM 
##

colQuantiles = function(x, probs=0.5){
  quants = matrix(NA, nrow=ncol(x), ncol=length(probs))
  for(i in 1:ncol(x)){
    quants[i,] = quantile(x[,i], probs)
  }
  colnames(quants) = probs
  return(quants)
}

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

rotate_mat = function(X, R, reverse=FALSE){
  if(!reverse) return( X%*%R )
  X%*%solve(R)
}

centre_mat = function(X, centre=NULL, reverse=FALSE){
  if(!reverse){
    if(!is.null(centre)){
      return( X - vec_rep(centre, nrow(X)))
    } 
    return( X - vec_rep(colMeans(X), nrow(X)) )
  }
  if(is.null(centre)){ stop("Reverse calculations require centre") }
  X + vec_rep(centre, nrow(X))
}

scale_mat = function(X, scale=NULL, reverse=FALSE){
  if(!reverse){
    if(!is.null(scale)){
      return( X / vec_rep(scale, nrow(X)) )
    }
    return( X / vec_rep(sqrt(diag(var(X))), nrow(X)) )
  }
  if(is.null(scale)){ stop("Reverse calculations require scale") }
  X * vec_rep(scale, nrow(X))
}

##
## DATA TRANSFORM PIPE
##

myData_transform = function(ctrl_mat, pat_mat=NULL){
  ctrl_mat = log(ctrl_mat)
  pca = prcomp(ctrl_mat, scale=FALSE, center=FALSE)
  ctrl_mat = pca$x
  # ctrl_mat = scale(ctrl_mat, center=TRUE, scale=TRUE)
  ctrl_mean = colMeans(ctrl_mat)
  ctrl_mat = centre_mat(ctrl_mat, ctrl_mean)
  ctrl_sd = sqrt(diag(var(ctrl_mat)))
  ctrl_mat = scale_mat(ctrl_mat, ctrl_sd)
  
  if(!is.null(pat_mat)){
    pat_mat = log(pat_mat)
    pat_mat = rotate_mat(pat_mat, pca$rotation)
    pat_mat = centre_mat(pat_mat, ctrl_mean)
    pat_mat = scale_mat(pat_mat, ctrl_sd)
    return(list(ctrl=ctrl_mat, pts=pat_mat))
  }
  ctrl_mat
}

back_transform = function(X, ctrl_raw=NULL, fulldat=NULL, folder=NULL, chan=NULL, pat=NULL){
  if( !is.null(fulldat) | !is.null(folder) | !is.null(chan) | !is.null(pat)){ 
    ctrl_raw = getData_mats(fulldat, chan=chan, pts=pat, ctrl_only=TRUE, data_transform=NULL)
  }
  if(!is.null(ctrl_raw)){
    ctrl_mat = log(ctrl_raw)
    pca = prcomp(ctrl_mat, scale=FALSE, center=FALSE)
    ctrl_mat = pca$x
    ctrl_mean = colMeans(ctrl_mat)
    ctrl_sd = sqrt(diag(var(ctrl_mat)))
    
    Xnew = scale_mat(X, scale=ctrl_sd, reverse=TRUE)
    Xnew = centre_mat(Xnew, centre=ctrl_mean, reverse=TRUE)
    return( rotate_mat(Xnew, R=pca$rotation, reverse=TRUE) )
  } else {
    stop("Must supply ctrl_raw or fulldat, folder, chan and pat")
  }
}

log_transform = function(ctrl_mat, pat_mat=NULL){
  if(!is.null(pat_mat)) return( list(ctrl=log(ctrl_mat), pts=log(pat_mat)) )
  log(ctrl_mat)
}

##
## SAVE & READ OUTPUT
##

output_saver = function(outroot, output, folder){
  split = strsplit(outroot, split="_")[[1]]
  froot = split[1]
  chan = split[2]
  pat = split[3]
  
  for(ctrl_pat in c("ctrl", "pat")){
    out_ctrlpat = ifelse(ctrl_pat=="ctrl", "CONTROL", pat)
    for(out_type in names(output[[ctrl_pat]])){
      filename = paste(froot, chan, out_ctrlpat, toupper(out_type), sep="_")
      write.table(output[[ctrl_pat]][[out_type]], paste0(file.path("Output", folder, filename), ".txt"),
                  row.names=FALSE, col.names=TRUE)
    }
  }
}

output_reader = function(folder, fulldat, chan, pat="CONTROL", out_type){
  outroot = paste(gsub(".RAW.txt", "", fulldat), chan, gsub("_",".", pat), sep="_")
  fp = paste0("./Output/", folder, "/", outroot, "_", out_type, ".txt")
  if(file.exists(fp)) return( read.table(fp, header=TRUE, stringsAsFactors=FALSE) )
  else stop("file does not exist")
}

##
## PLOTTING FUNCTIONS
##
 
MCMCplot_linreg = function(folder, fulldat, chan, pat="CONTROL", title="", lag=20){
    post = output_reader(folder, fulldat, chan, pat, out_type="POST")
    prior = output_reader(folder, fulldat, chan, pat, out_type="PRIOR")
    
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

colvector_gen = function(pts){
  colind = double(length(pts))
  pts_B = unique(gsub("_S.", "", pts))
  for(i in 1:length(pts_B)){
    colind[ gsub("_S.", "", pts) == pts_B[i] ] = i + 1
  }
  colind
}

pipost_plotter = function(fulldat, chan, folder, alpha=0.05){
  pts = getData_chanpats(fulldat)$patients
  pts_blocks = unique(gsub("_S.", "", pts))
  npat = length(pts)
  pis = list()
  for(pat in pts){
    pis[[pat]] = output_reader(folder, fulldat, chan, pat, out_type="POST")[,"probctrl"]
  }
  
  pat_labels = as.vector(rbind("", unique(gsub("_S.", "", gsub("P._", "", pts))), ""))
  
  title = paste0(substr(pts[1],1,2), " ", chan, "\n", gsub(".RAW.txt", "", fulldat) )
  # par(mar=c(6,4,4,8), xpd=TRUE)
  stripchart(pis, pch=20, method="jitter", vertical=TRUE, 
             col=rgb(t(col2rgb(palette()[colvector_gen(pts)]))/255, alpha=alpha), 
             group.names=rep("", length(pts)),
             main=title, ylim=c(0,1), ylab="like control proportion", 
             xlab="Patient Sample")
  if(!all(substr(pts, 4,4)==substr(pts,4,4)[1])){
    sep = sum(substr(pts, 4,4)==substr(pts,4,4)[1])+0.5
    arrows(x0=sep, x1=sep, y0=0, y1=1, lwd=3, lty=2, length=0, 
           col=myDarkGrey(1))
  }
  # legend(length(pts)+1, 1, lty=1, lwd=3, col=palette()[(1:length(pts_blocks)+2)], 
  #        legend=gsub("P._", "" , pts_blocks))
  text(1:length(pts), y=0, labels=pat_labels, pos=3, cex=1)
}

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

priorpost = function(ctrl_data, pat_data, prior=NULL, post, classifs_pat, title="",
                         mitochan, chan, pat, reverse_transform=NULL, ...){
  xlims = range(c(ctrl_data[,1], pat_data[,1]))
  ylims = range(c(ctrl_data[,2], pat_data[,2]))
  
  if(!is.null(reverse_transform)){
    prior = back_transform(prior, fulldat=fulldat, folder=folder, chan=chan, pat=pat)
    post = back_transform(post, fulldat=fulldat, folder=folder, chan=chan, pat=pat)
  }
  
  op = par(mfrow=c(1,2))
  if(!is.null(prior)){
    op = par(mfrow=c(2,2))
    plot(ctrl_data, pch=20, col=myDarkGrey(0.3), 
         xlab=paste0("log(", mitochan,")"), ylab=paste0("log(",chan,")"),
         xlim=xlims, ylim=ylims)
    points(pat_data, pch=20, col=myYellow(0.3))
    contours = percentiles(prior[,"comp.1.1."], prior[,"comp.2.1."])
    contour( contours$dens, levels=contours$levels, labels=contours$probs, add=TRUE, lwd=2 )
    
    plot(ctrl_data, pch=20, col=myDarkGrey(0.3), 
         xlab=paste0("log(", mitochan,")"), ylab=paste0("log(",chan,")"),
         xlim=xlims, ylim=ylims)
    points(pat_data, pch=20, col=myYellow(0.3))
    contours = percentiles(prior[["comp.1.2."]], prior[["comp.2.2."]])
    contour( contours$dens, levels=contours$levels, labels=contours$probs, add=TRUE, lwd=2 )
  }
  plot(ctrl_data, pch=20, col=myDarkGrey(0.3), 
       xlab=paste0("log(", mitochan,")"), ylab=paste0("log(",chan,")"),
       xlim=xlims, ylim=ylims)
  points(pat_data, pch=20, col=classcols(classifs_pat) )
  contours = percentiles(post[,"comp.1.1."], post[,"comp.2.1."])
  contour( contours$dens, levels=contours$levels, labels=contours$probs, add=TRUE, lwd=2 )
  
  plot(ctrl_data, pch=20, col=myDarkGrey(0.3), 
       xlab=paste0("log(", mitochan,")"), ylab=paste0("log(",chan,")"),
       xlim=xlims, ylim=ylims)
  points(pat_data, pch=20, col=classcols(classifs_pat) )
  contours = percentiles(post[,"comp.1.2."], post[,"comp.2.2."])
  contour( contours$dens, levels=contours$levels, labels=contours$probs, add=TRUE, lwd=2 )
  
  title(main=title, line=-2, outer=TRUE)
  
  par(op)
}

MCMCplot = function(folder, fulldat, chan, pat="CONTROL", title="", lag=20){
  post = output_reader(folder, fulldat, chan, pat, out_type="POST")
  prior = read.table(file.path("./Output", folder, "PRIOR.txt"), 
                     header=TRUE, stringsAsFactors=FALSE)

  col.names = colnames(post)
  n.chains = length(post)
  op = par(mfrow=c(2,3))
  for(param in col.names){
    post_vec = post[, param]
    plot(ts(post_vec), xlab="Iteration", ylab=paste(param), 
         main="", cex.lab=1.2)
    if(sum(post_vec==post_vec[1])!=length(post_vec)){
      acf(post_vec, xlab="lag index", ylab="ACF", main="",
          cex.lab=1.2, lag=lag)
    } else {
      plot(NA, type='n', xlim=c(0,lag), ylim=c(0,1), 
           xlab="lag index", ylab="ACF", main="")
    }
    plot(density(post_vec), lwd=2, col="blue", xlab=paste(param), ylab="Density",
         main="")
    if(param %in% colnames(prior)) lines(density(prior[,param]), lwd=2, col="green")
    
    title(main=title, line=-1, outer=TRUE, cex.main=1.6)
  }
  par(op)
}














































































