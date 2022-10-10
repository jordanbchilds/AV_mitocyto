library("parallel")
library("rstan")
source("./helper_functions.R", local=TRUE)

folder = "scaled_naiveBayes_GMM"
dir.create(file.path("./Output"), showWarnings = FALSE)
dir.create(file.path("./Output", folder), showWarnings = FALSE)

fulldats_raw = c(
  "IMV.E02.P01.PMT.M3243AG.QIF.7TR.RAW.txt",
  "IMV.E01.P02.PMT.M3243AG.QIF.TGo.RAW.txt",
  "IMV.E02.P01.PMT.M3243AG.QIF.TGo.RAW.txt",
  "IMV.E02.P01.PMT.M3243AG.QIF.TCF.RAW.txt"
)

# inferenecce 
inference = function(input){
  with(c(input),{
    data_mats = getData_mats(fulldat=fulldat, chan=chan, pts=pat, 
                             data_transform=myData_transform)
    Yctrl = data_mats$ctrl
    Ypat = data_mats$pts
    Nctrl = nrow(Yctrl)
    Npat = nrow(Ypat)
    
    output = stan("./nb_GMM_model.stan", data=prior_list, chains=1, 
                  iter=(MCMCOut*MCMCThin+MCMCBurnin), warmup=MCMCBurnin, thin=MCMCThin )
    
    outmat = as.matrix(output)
    outcols = colnames(outmat)
    
    classifs_df = outmat[, grepl("classif", outcols)]
    classifs_avg = colMeans(classifs_df)
    post = outmat[, !(grepl("classif", outcols)|grepl("lp__", outcols)|
                        grepl("probvec", outcols)|grepl("dens", outcols)|
                        grepl("z", outcols))]
    
    return( list(post=post, output=output, classifs=classifs_avg) )
  })
}

{
  p = 2
  # prior parameters
  mean_mu1 = double(p)
  mean_mu2 = double(p)
  var_mu1=  matrix(c(0.25, 0, 0, 0.25), ncol=p, nrow=p, byrow=TRUE) 
  var_mu2 = 5*diag(p)
  
  n_1 = 100
  S_1 = matrix(c(3^2, 0, 0, 3^2), nrow=p, ncol=p, byrow=TRUE)*(n_1-p-1)
  n_2 = 100
  S_2 = matrix(c(100,0,0,100), nrow=p, ncol=p, byrow=TRUE)*(n_2-p-1)
  
  alpha_pi = 1
  beta_pi = 1
}
prior_list = list(mean_mu1=mean_mu1, var_mu1=var_mu1,
                  mean_mu2=mean_mu2, var_mu2=var_mu2,
                  n_1=n_1, n_2=n_2, S_1=S_1, S_2=S_2,
                  alpha_pi=alpha_pi, beta_pi=beta_pi, 
                  D=2, K=2)

# inputs for the inference function - which data file name, chain length, thinning etc.
inputs = list()
{
  input0 = list()
  input0$MCMCOut = 2000
  input0$MCMCBurnin = 1000
  input0$MCMCThin = 1
  input0$n.chains = 1
  for(fulldat in fulldats_raw){
    froot = gsub(".RAW.txt", "", fulldat)
    chanpts = getData_chanpats(fulldat)
    for(chan in chanpts$channels){
      for(pat in chanpts$patients){
        outroot = paste(froot, chan, gsub("_",".",pat), sep="_")
        inputs[[outroot]] = input0
        inputs[[outroot]]$fulldat = fulldat
        inputs[[outroot]]$chan = chan
        inputs[[outroot]]$pat = pat
        inputs[[outroot]]$prior_list = prior_list
      } # pts
    } # chans
  } # fulldats 
}

ncores = 24
cl  = makeCluster(ncores) 
{
  clusterEvalQ(cl, {
    library("rstan")
    source("./helper_functions.R")
  })
  gmm_output = parLapply(cl, inputs, inference)
}
stopCluster(cl)

# save output
for(outroot in names(gmm_output)){
  output = gmm_output[[outroot]]
  
  write.table(output$post, file.path("./Output", folder, paste0(outroot,"_POST.txt")), row.names=FALSE, col.names=TRUE)
  write.table(output$classif, file.path("./Output", folder, paste0(outroot,"_CLASS.txt")), row.names=FALSE, col.names=TRUE)
}

### prior beliefs
output_prior = stan("./naiveBayes_GMM_prior.stan", data=prior_list, chains=1, 
                    iter=1e4, warmup=0, thin=1)

prior_all= as.matrix(output_prior)
prior = prior_all[, !grepl("_", colnames(prior_all))]

write.table(prior, file.path("./Output", folder, "PRIOR.txt"), 
            row.names=FALSE, col.names=TRUE)











