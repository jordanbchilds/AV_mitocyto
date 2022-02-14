library(rjags)
library(parallel)

myBlack = function(alpha) rgb(0,0,0, alpha)
myGrey = function(alpha) rgb(66/255,66/255,66/255, alpha)

myBlue = function(alpha) rgb(0,0,128/255, alpha)
myRed = function(alpha) rgb(1,0,0, alpha)
myGreen = function(alpha) rgb(0,100/255,0, alpha)

dir.create(file.path("OUTPUT"), showWarnings = FALSE)
dir.create(file.path("OUTPUT/GMM_joint"), showWarnings = FALSE)

fp_data = "Data_prepped.csv"
mc_df = read.csv(fp_data, header=TRUE)

pts_df = unique(mc_df$patient_id)

mc_data = mc_df
mc_data[mc_df[,"patient_id"] %in% pts_df[grep("C0", pts_df)],"patient_id"] = "control"

pts = unique(mc_data$patient_id)

mitochan = "VDAC"
channels = c("CYB","NDUFB8","MTCO1")

# model string
{
  modelstring = "
model {
  # fit to ctrl data
  for( i in 1:Nctrl ){
    Yctrl[i,] ~ dmnorm( mu[,1], tau[,,1] )
  }
  # fit to patient data
  for( j in 1:Npat ){ 
    z[j] ~ dbern(probdiff_pat)
    class[j] = 2 - z[j]
    Ypat[j,] ~ dmnorm(mu[,class[j]], tau[,,class[j]] )
  }
  # component one prior
  tau[1:2,1:2,1] ~ dwish(U_1, n_1)
  mu[1:2,1] ~ dmnorm(mu1_mean, mu1_prec)
  # component two prior
  tau[1:2,1:2,2] ~ dwish(U_2, n_2)
  mu[1:2,2] ~ dmnorm(mu2_mean, mu2_prec)
  
  # classification
  probdiff_pat ~ dbeta(alpha_pat, beta_pat)
  
  # predictive distribution
  pred[1:2,1] ~ dmnorm(mu[,1], tau[,,1])
  pred[1:2,2] ~ dmnorm(mu[,2], tau[,,2])
  
  zpred ~ dbern(probdiff_pat)
  compred = 2 - zpred
  Ypred ~ dmnorm(mu[,compred], tau[,,compred])
}
"
}

inference = function(input){
  # input: list(raw_data, chan, pat, mitochan)
  with(c(input),{
    outroot = paste0(chan,"__",pat)
    ## CONTROL DATA
    ctrl_data = raw_data[(raw_data$patient_id=="control"), ]
    Xctrl = log(ctrl_data[,paste(mitochan)])
    Yctrl = log(ctrl_data[,paste(chan)])
    XY_ctrl = cbind( Xctrl, Yctrl )
    Nctrl = nrow(XY_ctrl)
    
    ## PATIENT DATA
    pat_data = raw_data[(raw_data$patient_id==pat), ] 
    Xpat = log(pat_data[,paste(mitochan)])
    Ypat = log(pat_data[,paste(chan)])
    XY_pat = cbind(Xpat, Ypat)
    Npat = nrow(XY_pat)
    
    mu1_mean = colMeans(XY_ctrl)
    mu2_mean = mu1_mean
    mu1_prec = solve( matrix(c(0.1, 0.169, 0.169, 0.3), ncol=2, nrow=2, byrow=TRUE) ) # correlation of ~95.0%
    mu2_prec = 0.1*diag(2) 
    n_1 = 2000
    U_1 = matrix(c(0.2,0.335,0.335,0.6), nrow=2,ncol=2)*n_1 # correlation of ~95% 
    n_2 = 20
    U_2 = matrix(c(6,2,2,6),nrow=2,ncol=2)*n_2
    alpha_pat = 1
    beta_pat = 1
    
    data = list(Yctrl=XY_ctrl, Nctrl=Nctrl, Ypat=XY_pat, Npat=Npat,
                mu1_mean=mu1_mean, mu1_prec=mu1_prec,
                mu2_mean=mu2_mean, mu2_prec=mu2_prec, 
                n_1=n_1, n_2=n_2, U_1=U_1, U_2=U_2,
                alpha_pat=alpha_pat, beta_pat=beta_pat)
    
    data_prior = data
    data_prior$Yctrl = NULL
    data_prior$Ypat = NULL
    data_prior$Nctrl = 0
    data_prior$Npat = 0 
    
    model_jags = jags.model(textConnection(modelstring), data=data, n.chains=n.chains, n.adapt=MCMCBurnin,
                            quiet=TRUE)
    output_post = coda.samples(model_jags, variable.names=c("mu","tau","z","probdiff_pat","pred","Ypred"), 
                               n.iter=MCMCOutput, thin=MCMCThin)
    
    model_jags_prior = jags.model(textConnection(modelstring), data=data_prior, n.chains=n.chains, n.adapt=MCMCBurnin,
                                  quiet=TRUE)
    output_prior = coda.samples(model_jags_prior, variable.names=c("mu","tau","probdiff_pat","pred","Ypred"), 
                                n.iter=MCMCOutput, thin=1)
    
    output_post =  output_post[[1]]
    output_prior = output_prior[[1]]
    
    zpat = colnames(output_post[,grep("z", colnames(output_post))])
    zpat.split = strsplit(zpat, split="")
    zpat.vec = double(length(zpat.split))
    for(i in seq_along(zpat.split)){
      rr = zpat.split[[i]][ !zpat.split[[i]] %in% c("z","[","]") ]
      zpat.vec[i] = as.numeric(paste(rr, collapse=""))
    }
    names(zpat.vec) = zpat
    zpat.vec = sort(zpat.vec)
    
    class_post = output_post[, names(zpat.vec)]
    classifs_post = data.frame("patient_id"=1:ncol(class_post), "classifs"=colMeans(class_post))
    
    colnames(output_post) = colnames(output_post)
    colnames(output_prior) = colnames(output_prior)
    
    output_post = as.data.frame(output_post[,-grep("z", colnames(output_post))])
    output_prior = as.data.frame(output_prior[,-grep("z",colnames(output_prior))])
    
    return( list("post"=output_post, "prior"=output_prior, 
                 "classifs_post"=classifs_post) )
  })
}

# burn-in, chain length, thinning lag
MCMCBurnin = 1000
MCMCOutput = 5000 + MCMCBurnin
MCMCThin = 1
n.chains = 1

input_basics = list("mitochan"=mitochan, "raw_data"=mc_data,
                    "MCMCBurnin"=MCMCBurnin, "MCMCOutput"=MCMCOutput, 
                    "MCMCThin"=MCMCThin, "n.chains"=n.chains)

chanpat_list = list() 
for(chan in channels){
  for(pat in pts){
    chan_pat = paste(chan, pat, sep="_")
    chanpat_list[[chan_pat]] = input_basics
    chanpat_list[[chan_pat]][["chan"]] = chan
    chanpat_list[[chan_pat]][["pat"]] = pat
    
  }
}

cl  = makeCluster(24) 
clusterExport(cl, c("modelstring"))
clusterEvalQ(cl, {
  library("rjags")
})
cellpat_output = parLapply(cl, chanpat_list, inference )
stopCluster(cl)

# save output
for(cellpat in names(cellpat_output)){
  cellpat_list = cellpat_output[[cellpat]]
  # posterior 
  write.csv(cellpat_list$post, file=paste0("./OUTPUT/GMM_joint/",cellpat,"_POST.csv"))
  write.csv(cellpat_list$classifs_post, file=paste0("./OUTPUT/GMM_joint", cellpat,"_CLASSIFS.csv."))
  
  # prior
  write.csv(cellpat_list$prior, file=paste0("./OUTPUT/GMM_joint",cellpat,"_PRIOR.csv"))
}


















