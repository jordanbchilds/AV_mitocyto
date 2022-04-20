library("rstan")
library("MASS")
library("parallel")
library("mvtnorm")

dir.create(file.path("./OUTPUT"), showWarnings = FALSE)
dir.create(file.path("./OUTPUT/PCA_naiveBayes_GMM"), showWarnings = FALSE)

rbern = function(n,p) rbinom(n,1,p)

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
      Ypat_all = matrix(NA, nrow=1, ncol=2)
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

# inferenecce 
inference = function(input){
  with(c(input),{
    data_mats = getData_mats(chan=chan, pts=pat)
    Yctrl = data_mats$ctrl
    Ypat = data_mats$pat
    Nctrl = nrow(Yctrl)
    Npat = nrow(Ypat)
    
    
    p = 2
    # prior parameters
    mean_mu1 = double(p)
    mean_mu2 = double(p)
    var_mu1=  matrix(c(0.5, 0, 0, 0.1), ncol=p, nrow=p, byrow=TRUE) 
    var_mu2 = 2^2*diag(p)
    n_1 = 15
    S_1 = matrix(c(2, 0, 0, 0.25), nrow=p, ncol=p, byrow=TRUE)*(n_1-p-1)
    n_2 = 200
    S_2 = matrix(c(10^2,0,0,10^2), nrow=p, ncol=p, byrow=TRUE)*(n_2-p-1)
    alpha_pi = 1
    beta_pi = 1
    
    data_stan = list(Yctrl=Yctrl, Ypat=Ypat, Nctrl=Nctrl, Npat=Npat, 
                     mean_mu1=mean_mu1, var_mu1=var_mu1,
                     mean_mu2=mean_mu2, var_mu2=var_mu2,
                     n_1=n_1, n_2=n_2, S_1=S_1, S_2=S_2,
                     alpha_pi=alpha_pi, beta_pi=beta_pi, 
                     D=2, K=2)
    
    output = stan("nb_GMM_model.stan", data=data_stan, chains=1, 
                  iter=(MCMCOut*MCMCThin+MCMCBurnin), warmup=MCMCBurnin, thin=MCMCThin )
    
    outmat = as.matrix(output)
    outcols = colnames(outmat)
    
    classifs_df = outmat[, grepl("classif", outcols)]
    classifs_avg = colMeans(classifs_df)
    post = outmat[, !(grepl("classif", outcols)|grepl("lp__", outcols)|
                        grepl("probvec", outcols)|grepl("dens", outcols)|
                        grepl("z", outcols))]
    
    return( list(post=post, classif=classifs_avg) )
  })
}

# inputs for the inference function - which data file name, chain length, thinning etc.
inputs = list()
{
  input0 = list()
  input0$MCMCOut = 5000
  input0$MCMCBurnin = 1000
  input0$MCMCThin = 1
  input0$n.chains = 1
  for(chan in channels){
    for(pat in pts){
      outroot = paste(chan, pat, sep="_")
      inputs[[outroot]] = input0
      inputs[[outroot]]$chan = chan
      inputs[[outroot]]$pat = pat
    } # pts
  } # chans
}

ncores = 18
cl  = makeCluster(ncores) 
{
  clusterExport(cl, c("getData_mats"))
  clusterEvalQ(cl, {
    library("rstan")
  })
  gmm_output = parLapply(cl, inputs, inference)
}
stopCluster(cl)

# save output
for(index in 1:length(gmm_output)){
  output = gmm_output[[index]]
  chan = inputs[[index]]$chan
  pat = gsub("_",".",inputs[[index]]$pat)
  fulldat = inputs[[index]]$fulldat
  
  outroot = paste(chan, pat, sep="_")
  
  write.table(output$post, paste0("./OUTPUT/PCA_naiveBayes_GMM/",outroot,"_POST.txt"), row.names=FALSE, col.names=TRUE)
  write.table(output$classif, paste0("./OUTPUT/PCA_naiveBayes_GMM/",outroot,"_CLASS.txt"), row.names=FALSE, col.names=TRUE)
}

### prior beliefs
{
  p = 2
  # prior parameters
  mean_mu1 = double(p)
  mean_mu2 = double(p)
  var_mu1=  matrix(c(0.5, 0, 0, 0.1), ncol=p, nrow=p, byrow=TRUE) 
  var_mu2 = 1^2*diag(p)
  
  n_1 = 15
  S_1 = matrix(c(2, 0, 0, 0.25), nrow=p, ncol=p, byrow=TRUE)*(n_1-p-1)
  n_2 = 50
  S_2 = matrix(c(10^2,0,0,10^2), nrow=p, ncol=p, byrow=TRUE)*(n_2-p-1)
  
  alpha_pi = 1
  beta_pi = 1
}

data_prior = list(mean_mu1=mean_mu1, var_mu1=var_mu1,
                  mean_mu2=mean_mu2, var_mu2=var_mu2,
                  n_1=n_1, n_2=n_2, S_1=S_1, S_2=S_2,
                  alpha_pi=alpha_pi, beta_pi=beta_pi, 
                  D=2, K=2)

output_prior = stan("naiveBayes_GMM_prior.stan", data=data_prior, chains=1, 
                    iter=1e4, warmup=0, thin=1)

prior_all= as.matrix(output_prior)
prior_cols = colnames(prior_all)
prior = prior_all[, !(grepl("_", prior_cols)|
                      grepl("comp", prior_cols)|
                      grepl("pred", prior_cols))]
priorpred = prior_all[, grepl("comp", prior_cols)|grepl("pred",prior_cols)]

write.table(prior, paste0("./OUTPUT/PCA_naiveBayes_GMM/allData_PRIOR.txt"), 
            row.names=FALSE, col.names=TRUE)
write.table(priorpred, paste0("./OUTPUT/PCA_naiveBayes_GMM/allData_PRIORPRED.txt"), 
            row.names=FALSE, col.names=TRUE)






