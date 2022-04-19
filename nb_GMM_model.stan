data{
  int D; //number of dimensions
  int K; //number of mixtures
  
  vector[D] mean_mu1;
  vector[D] mean_mu2;
  cov_matrix[D] var_mu1;
  cov_matrix[D] var_mu2;
  int<lower=D> n_1;
  int<lower=D> n_2;
  cov_matrix[D] S_1;
  cov_matrix[D] S_2;
  
  real<lower=0> alpha_pi;
  real<lower=0> beta_pi;
  
  int<lower=0> Nctrl; //number of control data points
  int<lower=0> Npat; //number of patient fibre points
  matrix[Nctrl, D] Yctrl; //data matrix for control data
  matrix[Npat, D] Ypat; //data matrix for patient data
}
parameters{
  vector[D] mu1;
  vector[D] mu2;
  cov_matrix[D] V1; //precision of component one
  cov_matrix[D] V2;
  real<lower=0, upper=1> probctrl; //'like-ctrl' patient proportion
}
model{
  vector[K] tmp;
  
  mu1 ~ multi_normal(mean_mu1, var_mu1);
  mu2 ~ multi_normal(mean_mu2, var_mu2);
  V1 ~ inv_wishart(n_1, S_1);
  V2 ~ inv_wishart(n_2, S_2);
  probctrl ~ beta(alpha_pi, beta_pi);
  for(i in 1:Nctrl){
    target += multi_normal_lpdf(Yctrl[i,] | mu1, V1);
  }
  for(j in 1:Npat){
    tmp[1] = exp(log(probctrl) + multi_normal_lpdf(Ypat[j,] | mu1, V1));
    tmp[2] = exp(log(1-probctrl) + multi_normal_lpdf(Ypat[j,] | mu2, V2));
    target += log(sum(tmp));    
  }
}
generated quantities{
  matrix[D,K] comp;
  vector[D] pred;
  int<lower=0, upper=1> z;
  int<lower=0, upper=1> classif[Npat];
  vector<lower=0, upper=1>[Npat] probvec;
  matrix<lower=0>[Npat, K] dens;
  
  for(i in 1:Npat){
    dens[i,1] = exp(log(probctrl) + multi_normal_lpdf(Ypat[i,] | mu1, V1)); 
    dens[i,2] = exp(log(1-probctrl) + multi_normal_lpdf(Ypat[i,] | mu2, V2));
    probvec[i] = dens[i,1]/sum(dens[i,]);
  }
  
  classif = bernoulli_rng(probvec);
  
  z = bernoulli_rng(probctrl);
  if(z){
    pred = multi_normal_rng(mu1, V1);
  } else {
    pred = multi_normal_rng(mu2, V2);
  }
  comp[,1] = multi_normal_rng(mu1, V1);
  comp[,2] = multi_normal_rng(mu2, V2);
}



