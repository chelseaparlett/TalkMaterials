library(tidyverse)
library(rstan)
library(bayesplot)
library(RCurl)

# load data
x <- getURL("https://raw.githubusercontent.com/cmparlettpelleriti/TalkMaterials/master/TalkMaterials/SCiP2021/zoib_data.csv")

df <- read.csv(text = x)

# standata
model_code <- 'data {
int<lower=1> N; //number of data points
int<lower=1> J; //number of participants
int<lower=1> W; //number of words
int N_counts[3]; //0s, 0-1, 1s

int word[N]; //word

int<lower=0,upper=10> binary_pred1[N]; //binary_pred1 condition
real continuous_pred1[N];
real continuous_pred2[N];

int<lower=0> id[N];

real<lower=0, upper=1> brier[N]; // brier scores


}

 parameters {
   simplex[3] lambda; // probabilities for different mixtures
   
   real intercept;
   
   real continuous_pred1_b;
   real continuous_pred2_b;
   real binary_pred1_b;
 
  vector[J] u; //subject intercepts
  vector[W] w; //item intercepts
 
  real<lower=0> sigma_u; //subj sd
  real<lower=0> sigma_w; //item sd
  real<lower=0> kappa; // error of beta
 
 }

 model {
 real mu;
 real mu_p;
 
 N_counts ~ multinomial(lambda); // counts for each component (0, 0-1, 1)
 
 
 
 //priors
 
 u ~ normal(0, sigma_u); //subj random effects
 w ~ normal(0, sigma_w); //item random effects
 
 lambda ~ dirichlet(rep_vector(1.0,3));
 intercept ~ normal(0,1);
 continuous_pred1_b ~ normal(0,2);
 continuous_pred2_b ~ normal(0,2);
 binary_pred1_b ~ normal(0,2);

 
 sigma_u ~ normal(0,2); //subj sd
 sigma_w ~ normal(0,2); //item sd
 kappa ~ normal(0,1); // error of beta, truncated at 0
 
 
 // likelihood
 for (i in 1:N) {
   
    if (brier[i] == 0) {
      // likelihood when score is exactly 0
      target += log(lambda[1]);
      
    } else if (brier[i] == 1) {
      // likelihood when score is exactly 1
      target += log(lambda[3]);
      
    } else {
       // likelihood when score is between 0-1
      
      mu = intercept +
       binary_pred1_b*binary_pred1[i] +
       continuous_pred1_b*continuous_pred1[i] +
       continuous_pred2_b*continuous_pred2[i] +
       u[id[i]] +
       w[word[i]]; 
       
       mu_p = inv_logit(mu); // predicted value
       
      target += log(lambda[2]) + beta_proportion_lpdf(brier[i] | mu_p, kappa); // using proportion parameterization
    }
}

 
 }
 
'

model_code2 <- 'data {
int<lower=1> N; //number of data points
int<lower=1> J; //number of participants
int<lower=1> W; //number of words
int N_counts[3]; //0s, 0-1, 1s

int word[N]; //word

int<lower=0,upper=10> binary_pred1[N]; //binary_pred1 condition
real continuous_pred1[N];
real continuous_pred2[N];

int<lower=0> id[N];

real<lower=0, upper=1> brier[N]; // brier scores


}

 parameters {
   //simplex[3] lambda; // probabilities for different mixtures
   real p01; //probability of being 0 or 1
   real p1g0; //probability of being 1 given it is 1 or 0
   
   real intercept;
   
   real continuous_pred1_b;
   real continuous_pred2_b;
   real binary_pred1_b;
 
  vector[J] u; //subject intercepts
  vector[W] w; //item intercepts
 
  real<lower=0> sigma_u; //subj sd
  real<lower=0> sigma_w; //item sd
  real<lower=0> kappa; // error of beta
 
 }

 model {
 real mu;
 real mu_p;
 
 N_counts ~ multinomial(lambda); // counts for each component (0, 0-1, 1)
 
 
 
 //priors
 
 u ~ normal(0, sigma_u); //subj random effects
 w ~ normal(0, sigma_w); //item random effects
 
 lambda ~ dirichlet(rep_vector(1.0,3));
 intercept ~ normal(0,1);
 continuous_pred1_b ~ normal(0,2);
 continuous_pred2_b ~ normal(0,2);
 binary_pred1_b ~ normal(0,2);

 
 sigma_u ~ normal(0,2); //subj sd
 sigma_w ~ normal(0,2); //item sd
 kappa ~ normal(0,1); // error of beta, truncated at 0
 
 
 // likelihood
 for (i in 1:N) {
   
    if (brier[i] == 0) {
      // likelihood when score is exactly 0
      target += log(lambda[1]);
      
    } else if (brier[i] == 1) {
      // likelihood when score is exactly 1
      target += log(lambda[3]);
      
    } else {
       // likelihood when score is between 0-1
      
      mu = intercept +
       binary_pred1_b*binary_pred1[i] +
       continuous_pred1_b*continuous_pred1[i] +
       continuous_pred2_b*continuous_pred2[i] +
       u[id[i]] +
       w[word[i]]; 
       
       mu_p = inv_logit(mu); // predicted value
       
      target += log(lambda[2]) + beta_proportion_lpdf(brier[i] | mu_p, kappa); // using proportion parameterization
    }
}

 
 }
 
'
data_list <- list(N = nrow(df),
                  J = length(unique(df$participant)),
                  W = length(unique(df$item)),
                  N_counts = c(sum(df$out == 0),
                               sum(df$out > 0 & df$out < 1),
                               sum(df$out == 1)),
                  order = df$order,
                  word = df$item,
                  binary_pred1 = df$binary_pred,
                  continuous_pred1 = df$cont_pred1,
                  continuous_pred2 = df$cont_pred2,
                  id = df$participant,
                  brier = df$out)

fitted_model <- stan(model_code = model_code,
                     data = data_list,
                     iter = 4000)