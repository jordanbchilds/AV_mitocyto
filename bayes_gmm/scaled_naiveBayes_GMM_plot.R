library("MASS")
source("./helper_functions.R", local = TRUE)

folder = "scaled_naiveBayes_GMM"

dir.create(file.path("./PDF"), showWarnings = FALSE)
dir.create(file.path("./PDF", folder), showWarnings = FALSE)

fulldats_raw = c(
  "IMV.E02.P01.PMT.M3243AG.QIF.7TR.RAW.txt",
  "IMV.E01.P02.PMT.M3243AG.QIF.TGo.RAW.txt",
  "IMV.E02.P01.PMT.M3243AG.QIF.TGo.RAW.txt",
  "IMV.E02.P01.PMT.M3243AG.QIF.TCF.RAW.txt"
)

prior = read.table(file.path("./Output", folder, "PRIOR.txt"), 
                   header=TRUE, stringsAsFactors=FALSE)

pdf(file.path("./PDF", folder, "MCMC.pdf"), width=13, height=8)
{ 
  for(fulldat in fulldats_raw[1:2]){
    chanpts = getData_chanpats(fulldat)
    pts = gsub("_", ".", chanpts$patients)
      for(chan in chanpts$channels){
        for(pat in pts){
          MCMCplot(folder, fulldat, chan, pat, lag=100,
                   title=paste(chan, pat))
        }
      }
  }
}
dev.off()

pdf(file.path("./PDF", folder, "model_post_RAW.pdf"), width=13, height=8)
{
  for(fulldat in fulldats_raw[1:2]){
    chanpts = getData_chanpats(fulldat)
    for(chan in chanpts$channels){
      output_reader0 = function(pat="CONTROL", out_type){
        output_reader(folder=folder, fulldat=fulldat, chan=chan, pat, out_type)
      }
      ctrl_data =  getData_mats(fulldat, chan=chan, ctrl_only=TRUE, 
                                data_transform=myData_transform)

      for(pat in chanpts$patients){
        pat_data = getData_mats(fulldat, chan=chan, pts=pat, 
                                data_transform=myData_transform)$pts
        
        priorpost(ctrl_data=ctrl_data, pat_data=pat_data,
                  classif=output_reader0(pat=pat, out_type="CLASS")[[1]],
                  prior=prior, post=output_reader0(pat=pat, "POST"),
                  chan=chan, mitochan="VDAC1", 
                  title=paste(gsub(".RAW.txt", "", fulldat), "\n", chan, pat))
      }
    }
  }
}
dev.off()

pdf(file.path("./PDF", folder, "model_post.pdf"), width=13, height=8)
{
  for(fulldat in fulldats_raw[1:2]){
    chanpts = getData_chanpats(fulldat)
    for(chan in chanpts$channels){
      output_reader0 = function(pat="CONTROL", out_type){
        output_reader(folder=folder, fulldat=fulldat, chan=chan, pat, out_type)
      }
      ctrl_data =  getData_mats(fulldat, chan=chan, ctrl_only=TRUE, 
                                data_transform=myData_transform)
      
      for(pat in chanpts$patients){
        pat_data = getData_mats(fulldat, chan=chan, pts=pat, 
                                data_transform=myData_transform)$pts
        
        priorpost(ctrl_data=ctrl_data, pat_data=pat_data,
                  classif=output_reader0(pat=pat, out_type="CLASS")[[1]],
                  prior=prior, post=output_reader0(pat=pat, "POST"),
                  chan=chan, mitochan="VDAC1", title=paste(gsub(".RAW.txt", "", fulldat), "\n", chan, pat),
                  reverse_transform=back_transform, fulldat=fulldat, folder=folder, pat=pat)
      }
    }
  }
}
dev.off()

pdf(file.path("./PDF", folder, "pi_post.pdf"), width=13, height=8)
{
  par(mfrow=c(1,1))
  for(fulldat in fulldats_raw[1:2]){
    chans = getData_chanpats(fulldat)$channels
      for(chan in chans){
        pipost_plotter(fulldat, chan, folder=folder ) 
      }
    }
}
dev.off()












