library(MASS)
source("./helper_functions.R", local = TRUE)

folder = "linReg_mod2"

dir.create("./PDF", showWarnings=FALSE)
dir.create(file.path("./PDF",folder), showWarnings=FALSE)

fulldats_raw = c(
  "IMV.E02.P01.PMT.M3243AG.QIF.7TR.RAW.txt",
  "IMV.E01.P02.PMT.M3243AG.QIF.TGo.RAW.txt",
  "IMV.E02.P01.PMT.M3243AG.QIF.TGo.RAW.txt",
  "IMV.E02.P01.PMT.M3243AG.QIF.TCF.RAW.txt"
)

######################
### the plots
######################

pdf(file.path("./PDF", folder, "MCMC.pdf"), width=13, height=8)
{ 
  for(fulldat in fulldats_raw){
    froot = gsub(".RAW.txt", "", fulldat)
    pts = gsub("_", ".", getData_chanpats(fulldat)$patients)
    
    for(chan in c("NDUFB8", "MTCO1")){
      outroot_ctrl = paste(froot, chan, "CONTROL", sep="_")
      
      MCMCplot(folder, fulldat, chan, lag=100, 
               title=paste(chan, "CONTROL"))
      for(pat in pts){
        outroot_pat = paste(froot, chan, pat, sep="_")
        MCMCplot(folder, fulldat, chan, pat, lag=100,
                 title=paste(chan, pat))
      }
    }
  }
}
dev.off()

pdf(file.path("./PDF", folder, "model_post.pdf"), width=13, height=8)
{
  for(fulldat in fulldats_raw){
    chanpts = getData_chanpats(fulldat)
    
    for(chan in chanpts$channels){
      output_reader0 = function(pat="CONTROL", out_type){
        output_reader(folder=folder, fulldat=fulldat, chan=chan, pat, out_type)
      }
    ctrl_data =  getData_mats(fulldat, chan=chan)$ctrl
    pts_data = getData_mats(fulldat, chan=chan)$pts
    xlims = range(c(ctrl_data[,1], pts_data[,1]))
    ylims = range(c(ctrl_data[,2], pts_data[,2]))
      
    # priorpost(ctrl_data=ctrl_data, 
    #           priorpred=output_reader0(out_type="PRIORPRED"),
    #           postpred=output_reader0(out_type="POSTPRED"),
    #           chan=chan,title=paste(chan, "CONTROL"),
    #           xlims=xlims, ylims=ylims)
    #   
    for(pat in chanpts$patients){
      pat_data = getData_mats(fulldat, chan=chan, pts=pat)$pts
        
      priorpost(ctrl_data=ctrl_data, pat_data=pat_data,
                classif=output_reader0(pat=pat, out_type="CLASSIF")[[1]],
                priorpred=output_reader0(pat=pat, "PRIORPRED"), 
                postpred=output_reader0(pat=pat, "POSTPRED"),
                chan=chan, mitochan="VDAC1", title=paste(fulldat, "\n", chan, pat),
                xlims=xlims, ylims=ylims)
      }
    }
  }
}
dev.off()

pdf(file.path("./PDF", folder, "pi_post.pdf"), width=13, height=8)
{
  par(mfrow=c(1,1))
  for(fulldat in fulldats_raw){
    for(chan in c("MTCO1", "NDUFB8")){
      pipost_plotter(fulldat, chan, folder=folder ) 
    }
  }
}
dev.off()














