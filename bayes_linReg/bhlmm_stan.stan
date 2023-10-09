//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//
functions{
  // sample from truncated normal distribution - lower truncated
  real normal_lb_rng(real mu, real sigma, real lb){
    real p_lb = normal_cdf(lb, mu, sigma);
    real u = uniform_rng(p_lb, 1.0);
    real y = mu + sigma*Phi(u);
    return y;
  }
}
data{
  int D; // number of dimensions
  int K; // number of mixtures
  int M; // number of daughters in the hierarchy
  int nSyn; // number of predictive inputs
  
  vector[nSyn] xSyn;
  
  real mean_mu_m;
  real mean_mu_c;
  real<lower=0> prec_mu_m;
  real<lower=0> prec_mu_c;
  real<lower=0> shape_tau_m;
  real<lower=0> rate_tau_m;
  real<lower=0> shape_tau_c;
  real<lower=0> rate_tau_c;
  real<lower=0> shape_tau;
  real<lower=0> rate_tau;
  
  real slope_lb;
  
  real<lower=0> alpha_pi;
  real<lower=0> beta_pi;
  
  real tau_def;
  
  int<lower=0> nCtrl; //number of control data points
  int<lower=0> nPat; //number of patient fibre points
  int<lower=0> nPts; //number patients subjects in the dataset
  matrix[nCtrl, D] ctrl_mat; //data matrix for control data
  matrix[nPat, D] pat_mat; //data matrix for patient data
  int<lower=1, upper=M> ctrlIndex[nCtrl]; // vector indexing which observation belongs to which patient
}
parameters{
  vector<lower=0>[M] m; // slopes
  vector[M] c; // intercepts
  real<lower=0> tau_norm;
  
  real mu_m;
  real mu_c;
  real<lower=0> tau_m;
  real<lower=0> tau_c;
  
  real<lower=0, upper=1> probdiff;
}
model{
  real sigma_norm;
  real sigma_def;
  
  mu_m ~ normal(mean_mu_m, 1/sqrt(prec_mu_m));
  mu_c ~ normal(mean_mu_c, 1/sqrt(prec_mu_c));
  tau_m ~ gamma(shape_tau_m, rate_tau_m);
  tau_c ~ gamma(shape_tau_c, rate_tau_c);
  
  for( i in 1:M ){
    m[i] ~ normal(mu_m, 1/sqrt(tau_m)) T[slope_lb,];
    c[i] ~ normal(mu_c, 1/sqrt(tau_c));
  }
  
  probdiff ~ uniform(alpha_pi, beta_pi);
  tau_norm ~ gamma(shape_tau, rate_tau);
  
  sigma_norm = 1/sqrt(tau_norm);
  sigma_def = 1/sqrt(tau_def);
  
  for(i in 1:nCtrl){
    target += normal_lpdf(ctrl_mat[i,2] | m[ctrlIndex[i]]*ctrl_mat[i,1]+c[ctrlIndex[i]], sigma_norm) ;
  }
  for(j in 1:nPat){
    target += log_mix(probdiff,
                      normal_lpdf(pat_mat[j,2] | m[M]*pat_mat[j,1]+c[M], sigma_def), 
                      normal_lpdf(pat_mat[j,2] | m[M]*pat_mat[j,1]+c[M], sigma_norm));    
  }
}
generated quantities{
  // tmp quantities to reduce repeated calculations
  real log_probDef_tmp;
  real log_probHealthy_tmp;
  real sigma_norm_tmp;
  real sigma_def_tmp;
  
  // classifications and predictions
  int<lower=0, upper=1> classif[nPat];
  vector<lower=0, upper=1>[nPat] probvec;
  matrix<lower=0>[nPat, K] dens;
  vector[nSyn] yPred;
  real m_pred;
  real c_pred;

  // prior draws
  vector<lower=0>[M] m_prior; 
  vector[M] c_prior; 
  real<lower=0> tau_norm_prior;
  real mu_m_prior;
  real mu_c_prior;
  real<lower=0> tau_m_prior;
  real<lower=0> tau_c_prior;
  real<lower=0, upper=1> probdiff_prior;
  vector[nSyn] yPred_prior;
  real m_pred_prior;
  real c_pred_prior;
  
  // the log-likelihood
  real log_lik;
  
  // classifier
  log_probDef_tmp = log(probdiff);
  log_probHealthy_tmp = log( 1 - probdiff );
  sigma_norm_tmp = 1/sqrt(tau_norm);
  sigma_def_tmp = 1/sqrt(tau_def);
  
  for(j in 1:nPat){
    dens[j,1] = exp( log_probHealthy_tmp + normal_lpdf(pat_mat[j,2] | m[M]*pat_mat[j,1]+c[M], sigma_norm_tmp));
    dens[j,2] = exp( log_probDef_tmp + normal_lpdf(pat_mat[j,2] | m[M]*pat_mat[j,1]+c[M], sigma_def_tmp));
    probvec[j] = dens[j,2]/sum(dens[j,]);
  }
  
  classif = bernoulli_rng(probvec);
  
  // prior draws
  mu_m_prior = normal_rng(mean_mu_m, 1/sqrt(prec_mu_m));
  mu_c_prior = normal_rng(mean_mu_c, 1/sqrt(prec_mu_c));
  tau_m_prior = gamma_rng(shape_tau_m, rate_tau_m);
  tau_c_prior = gamma_rng(shape_tau_c, rate_tau_c);
  for( i in 1:M ){
    m_prior[i] = normal_lb_rng(mu_m_prior, 1/sqrt(tau_m_prior), slope_lb);
    c_prior[i] = normal_rng(mu_c_prior, 1/sqrt(tau_c_prior));
  }
  m_pred_prior = normal_lb_rng(mu_m_prior, 1/sqrt(tau_m_prior), slope_lb);
  c_pred_prior = normal_rng(mu_c_prior, 1/sqrt(tau_c_prior));
  tau_norm_prior = gamma_rng(shape_tau, rate_tau);
  probdiff_prior = uniform_rng(alpha_pi, beta_pi);
  
  // posterior parent slope and inter distribution
  m_pred = normal_lb_rng(mu_m, 1/sqrt(tau_m), slope_lb);
  c_pred = normal_rng(mu_c, 1/sqrt(tau_c));
  
  // preidictive model for patient data
  for(k in 1:nSyn){
    yPred[k] = normal_rng( m[M]*xSyn[k]+c[M], sigma_norm_tmp );
    yPred_prior[k] = normal_rng( m_prior[M]*xSyn[k]+c_prior[M], 1/sqrt(tau_norm_prior) );
  }
  
  // the log-likelihood
  log_lik = 0.0;
  for(i in 1:nCtrl){
    log_lik += normal_lpdf(ctrl_mat[i,2] | m[ctrlIndex[i]]*ctrl_mat[i,1]+c[ctrlIndex[i]], sigma_norm_tmp) ;
  }
  for(j in 1:nPat){
    log_lik += classif[j] ? normal_lpdf(pat_mat[j,2] | m[M]*pat_mat[j,1]+c[M], sigma_def_tmp) : normal_lpdf(pat_mat[j,2] | m[M]*pat_mat[j,1]+c[M], sigma_norm_tmp);
  }
}


