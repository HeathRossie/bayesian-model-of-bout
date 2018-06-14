data{
  
  int N_VI;
  int N_VR;
  int time_VI[N_VI];
  int resp_VI[N_VI];
  int time_VR[N_VR];
  int resp_VR[N_VR];
  int N_time_VI;
  int N_time_VR;
  
}

parameters{
  
  real lambda_VI_[N_time_VI];
  real lambda_VR_[N_time_VR];
  real<lower=0.0001> sigma1[2];
  real<lower=0.0001> sigma2[2];
  real q_VI_[N_time_VI];
  real q_VR_[N_time_VR];
  
}

transformed parameters{
  
  real<lower=0> lambda_VI[N_time_VI];
  real<lower=0> lambda_VR[N_time_VR];
  real<lower=0,upper=1> q_VI[N_time_VI];
  real<lower=0,upper=1> q_VR[N_time_VR];
  
  for(i in 1:N_time_VI){
    lambda_VI[i] = exp(lambda_VI_[i]);
  }
  
  for(i in 1:N_time_VR){
    lambda_VR[i] = exp(lambda_VR_[i]);
  }
  
  for(i in 1:N_time_VI){
    q_VI[i] = inv_logit(q_VI_[i]);
  }
  
  for(i in 1:N_time_VR){
    q_VR[i] = inv_logit(q_VR_[i]);
  }
  
}

model{
  
  // state model
  for(i in 1:N_time_VI){
    if(time_VI[i]==1){
      lambda_VI_[i] ~ normal(0, 100);
      q_VI_[i] ~ normal(0, 100);
    }else{
      lambda_VI_[i] ~ normal(lambda_VI_[i-1], sigma1[1]);
      q_VI_[i] ~ normal(q_VI_[i-1], sigma2[1]);
    }
  }
  
  for(i in 1:N_time_VR){
    if(time_VR[i]==1){
      lambda_VR_[i] ~ normal(0, 100);
      q_VR_[i] ~ normal(0, 100);
    }else{
      lambda_VR_[i] ~ normal(lambda_VR_[i-1], sigma1[2]);
      q_VR_[i] ~ normal(q_VR_[i-1], sigma2[2]);
    }
  }
  
  // observation model
  for(i in 1:N_VI){
    if(resp_VI[i]==0){
      target += log_sum_exp(bernoulli_lpmf(0|q_VI[time_VI[i]]),
                            bernoulli_lpmf(1|q_VI[time_VI[i]]) + 
                              poisson_lpmf(0|lambda_VI[time_VI[i]]));
    }else{
      target += bernoulli_lpmf(1|q_VI[time_VI[i]]) + 
        poisson_lpmf(resp_VI[i]|lambda_VI[time_VI[i]]);
    }
  }
  
  for(i in 1:N_VR){
    if(resp_VR[i]==0){
      target += log_sum_exp(bernoulli_lpmf(0|q_VR[time_VR[i]]),
                            bernoulli_lpmf(1|q_VR[time_VR[i]]) + 
                              poisson_lpmf(0|lambda_VR[time_VR[i]]));
    }else{
      target += bernoulli_lpmf(1|q_VR[time_VR[i]]) + 
        poisson_lpmf(resp_VR[i]|lambda_VR[time_VR[i]]);
    }
  }
}

generated quantities{
  int pred_VI[N_time_VI];
  int pred_VR[N_time_VR];
  
  // the reported vallues in the manuscript
  int pred_VI2[N_VI];
  int pred_VR2[N_VR];
  
  real mean_lambda_VI;
  real mean_lambda_VR;
  real mean_q_VI;
  real mean_q_VR;
  
  for(i in 1:N_time_VI){
    pred_VI[i] = bernoulli_rng(q_VI[i]) * 
      poisson_rng(lambda_VI[i]);
  }
  
  for(i in 1:N_time_VR){
    pred_VR[i] = bernoulli_rng(q_VI[i]) *
      poisson_rng(lambda_VR[i]);
  }
  
  for(i in 1:N_VI){
    pred_VI2[i] = bernoulli_rng(q_VI[time_VI[i]]) * 
      poisson_rng(lambda_VI[time_VI[i]]);
  }
  
  for(i in 1:N_VR){
    pred_VR2[i] = bernoulli_rng(q_VI[time_VR[i]]) *
      poisson_rng(lambda_VR[time_VR[i]]);
  }
  
  mean_lambda_VI = mean(lambda_VI);
  mean_lambda_VR = mean(lambda_VR);
  mean_q_VI = mean(q_VI);
  mean_q_VR = mean(q_VR);
  
}
