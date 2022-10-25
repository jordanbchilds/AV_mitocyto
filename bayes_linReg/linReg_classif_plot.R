library(MASS)
source("helper_functions.R", local = TRUE)

folder = "linReg_classifier"

dir.create("PDF", showWarnings=FALSE)
dir.create(file.path("./PDF",folder), showWarnings=FALSE)

cord = c("NDUFB8", "CYB", "MTCO1")
chlabs = c("CI", "CII","CIV")
names(chlabs) = cord
mitochan = "VDAC1"

dat = read.csv(file.path("..", "Data_prepped.csv"), stringsAsFactors=FALSE)

ptsAll = unique(dat$patient_id)
pts = ptsAll[grepl("P", ptsAll)]

######################
### the plots
######################

pdf(file.path("PDF", folder, "MCMC.pdf"), width=13, height=8)
{ 
  for(chan in cord){
    outroot_ctrl = paste(chan, "CONTROL", sep="_")
      
    # MCMCplot(folder, chan, lag=100, 
    #            title=paste(chan, "CONTROL"))
    for(pat in pts){
      outroot_pat = paste(chan, pat, sep="_")
      MCMCplot(folder, chan, pat, lag=100,
               title=paste(chan, pat))
    } # channel
  } # patient
}
dev.off()

pdf(file.path("PDF", folder, "model_post.pdf"), width=13, height=8)
{
  op = par(mfrow=c(1,1), mar=c(6,6,6,3), cex.main=2, cex.lab=2, cex.axis=1.5)
  for(chan in cord){
    data = getData_mats(chan=chan)
    ctrl_mat =  data$ctrl
    xlims = range(c(ctrl_mat[,1], data$pts[,1]))
    ylims = range(c(ctrl_mat[,2], data$pts[,2]))
    for(pat in pts){
      pat_mat = getData_mats(chan=chan, pts=pat)$pts
          
      priorpost(ctrl_data=ctrl_mat, pat_data=pat_mat,
                classif=output_reader(folder, chan, pat, out_type="CLASSIF")[[1]],
                priorpred=output_reader(folder, chan, pat, "PRIORPRED"), 
                postpred=output_reader(folder, chan, pat, "POSTPRED"),
                chan=chan, mitochan="VDAC1", title=paste("\n", chan, pat),
                xlims=xlims, ylims=ylims)
    } # patients
  } # channels
}
dev.off()

pdf(file.path("PDF", folder, "pi_post.pdf"), width=13, height=8)
{
  op = par(mfrow=c(1,1), mar=c(6,6,6,3), cex.main=2, cex.lab=2, cex.axis=1.5)
  for(chan in cord){
    pipost_plotter(chan, pts=pts, folder=folder) 
  }
  par(op)
}
dev.off()

pdf(file.path("PDF", folder, "piall_post.pdf"), width=13, height=8)
{
  op = par(mfrow=c(1,1), mar=c(6,6,6,3), cex.main=2, cex.lab=2, cex.axis=1.5)
  pipost_plotter_v2(cord, pts, folder)
  par(op)
}
dev.off()











