# --- convergence check
library("tidyr")
library("data.table")
library("coda")
library("mcmcse")
library("stringr")

library("analysis2Dmito")

mitochan = "VDAC"
channels = c("NDUFB8", "CYB", "MTCO1")
nChan = length(channels)

raw_data = read.csv("../Data_prepped.csv", header=TRUE)

raw_data = raw_data[,c("ID", "patient_id", mitochan, channels)]
colnames(raw_data) = c("fibreID", "sampleID", mitochan, channels)

data_lng = pivot_longer(raw_data, cols=c(mitochan, channels), names_to="Channel", values_to="Value")
data_lng = as.data.frame(data_lng)
data_lng$Value = log( data_lng$Value)

sbjIDs = unique(data_lng$sampleID)
ctrlIDs = grep("C", sbjIDs, value=TRUE)
ptsIDs = sbjIDs[!(sbjIDs %in% ctrlIDs)]

nSbjs = length(sbjIDs)
nCtrls = length(ctrlIDs)
nPts = length( ptsIDs )

folder = "Output_BH"
for( chan in channels ){
  for( chain in 1:5 ){
    root = paste0(chan, "_chain0", chain, "_POST.csv")
    post = as.data.frame( fread( file.path(folder, root) ) )
    
    print(root)
    print( dim(post) )
  }
}

for( chan in channels ){
    for( chain in 1:5){
      root = paste0(chan, "_chain", str_pad(chain, width=2, pad="0"))
      
      fn_post = file.path("Output_BH", paste0(root, "_POST.csv"))
      
      post = as.data.frame( fread( fn_post ))
      print( paste(root, nrow(fn_post) ) )
      print(effectiveSize( post ))
      print(multiESS(post))
    }
}

dir.create("PDF_BH")

# ------ MCMC plot 
pdf(file.path("PDF_BH", "MCMCplot.pdf"), width=13, height=8)
{
  for(chan in channels){
      for(chain in 1:5){
        outroot = paste0(chan, "_chain", str_pad(chain, width=2, pad="0"))
        fn_post = file.path("Output_BH", paste0(outroot, "_POST.csv"))
        
        if(file.exists(fn_post)){
          post = as.data.frame(fread(fn_post))
          
          fn_prior = file.path("Output_BH", paste0(outroot, "_PRIOR.csv"))
          prior = as.data.frame(fread(fn_prior))
          analysis2Dmito::MCMCplot(post, prior, nRow=3, lag=100)
          title(main=outroot, outer=TRUE, line=-1)
        }
      }
    }
}
dev.off()

# ------ postPlot
pdf(file.path("PDF_BH", "postPlot.pdf"), width=13, height=8)
{
  for(chan in channels){
      for(chain in 1:5){
        outroot = paste0(chan, "_chain", str_pad(chain, width=2, pad="0"))
          fn_post = file.path("Output_BH", paste0(outroot, "_POST.csv"))
          
          if(file.exists(fn_post)){
            post = as.data.frame(fread(fn_post))
            fn_prior = file.path("Output_BH", paste0(outroot, "_PRIOR.csv"))
            prior = as.data.frame(fread(fn_prior))
            
            fn_postpred = file.path("Output_BH", paste0(outroot, "_POSTPRED.csv"))
            postpred = as.data.frame(fread(fn_postpred))
            
            fn_class = file.path("Output_BH", paste0(outroot, "_CLASSIF.csv"))
            class = apply(as.matrix(fread(fn_class)), 2, mean)
            
            for( pat in ptsIDs ){
              pat_dat = data_lng[data_lng$sampleID %in% c(ctrlIDs, pat), ]
              
              dataMats = analysis2Dmito::getData_mats(data=pat_dat, 
                                                      channels=c(mitochan, chan), 
                                                      ctrlID=ctrlIDs,
                                                      pts=pat)
              
              clab = c( paste0( "c[", ctrlIDs, "]"), paste0("c[", pat, "]") )
              mlab = c( paste0( "m[", ctrlIDs, "]"), paste0("m[", pat, "]") )
              
              patpost_cols = grep("m\\[", colnames(post), invert=TRUE, value=TRUE )
              patpost_cols = grep("c\\[", patpost_cols, invert=TRUE, value=TRUE )
              patpost_cols = grep("pi\\[", patpost_cols, invert=TRUE, value=TRUE)
              patpost_cols = c(patpost_cols, mlab, clab, paste0("pi[", pat, "]") )
              
              pat_post = post[, patpost_cols]
              colnames(pat_post) = c("tau_c", "tau_m", "mu_c", "mu_m", "tau_norm",
                                     "m_pred", "c_pred", 
                                     paste0("m[", 1:(nCtrls + 1), "]"), 
                                     paste0("c[", 1:(nCtrls + 1), "]"), 
                                     "probdiff")
              
              pat_prior = prior[, patpost_cols]
              colnames(pat_prior) = c("tau_c", "tau_m", "mu_c", "mu_m", "tau_norm",
                                     "m_pred", "c_pred", 
                                     paste0("m[", 1:(nCtrls + 1), "]"), 
                                     paste0("c[", 1:(nCtrls + 1), "]"), 
                                     "probdiff")
                        
              pat_postpred = postpred[, c("mitochan", grep(pat, colnames(postpred), value=TRUE) ) ]
              colnames(pat_postpred) = c("mitochan", "lwrNorm", "medNorm", "uprNorm")
              
              pat_class = class[ grepl(pat, names(class)) ]
              
              op = par(mfrow=c(3,3), mar=c(4,4,1,1), cex.main=2, cex.lab=1.5, cex.axis=1.5)
              analysis2Dmito::postPlot(post=pat_post, 
                                       prior=pat_prior,
                                       postpred=pat_postpred, 
                                       classifs=pat_class,
                                       dataMats = dataMats,
                                       var.names=c("mu_m", "tau_m", "tau_norm", "mu_c", "tau_c", "probdiff", "m", "c"),
                                       mitochan=paste0("log(", mitochan, ")"),
                                       chan=paste0("log(", chan, ")"))
              title(main=outroot, line=-1, outer=TRUE)
              par(op)
       }
      }
    }
  }
}
dev.off()



