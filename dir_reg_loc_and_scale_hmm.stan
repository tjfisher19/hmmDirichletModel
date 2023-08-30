// This code implements the Hidden Markov Model with a Dirichlet response
//
// Here, the shape parameters from the Dirichlet distribution are 
// decomposed into a location and scale term, each is then modeled
// with a log-link linear model (Dirichlet regression idea).
//
// We constrain the Markov Model so that it will begin in State 1
// and stay there for P time units (one year in our data application
// where we model with predictors mapping to the seasonal terms).
// Thus, providing a burn-in type period for the initial parameters 
// to be set. 
//
// Snagged base of the code from here and then modified using Rstan manual:
//  https://github.com/stan-dev/stancon_talks/blob/master/2018/Contributed-Talks/04_damiano/stan/iohmm_reg.stan
//
// by: Thomas J. Fisher
//

functions {
  vector normalize(vector x) {
    return x / sum(x);
  }
}

data {
  int<lower=1> T;                   // number of observations (length)
  int<lower=1> M;                   // number of hidden states
  int<lower=1> P;                   // dimension of the input vector/covariates
  int<lower=1> K;                   // Dimension of simplex

  vector[K] y[T];                   // output (scalar)
  matrix[T,P] X;                    // inputs (design matrix)

  // priors for transition matrix
  real<lower = 0.0> rho_a;
  real<lower = 0.0> rho_b;
  
  // priors for variance of model coefficients
  real<lower = 0.0> b_loc_sig;
  real<lower = 0.0> b_scale_sig;
}

parameters {
  // Discrete state model
  vector<lower=0.0, upper=1.0>[M-1] prob_remain;   // prob of remaining in current state

  // Parameters for regression
  real b_loc[M,P,K];                        // regressors on covariates on location parameter
  real b_scale[M,P];                        // regressor for spread term

}

transformed parameters {
  vector[M] logalpha[T];                    // Forward algorithm state probabilities
  vector[M] logoblik[T];                    // Likelihood value at time t in state m
  simplex[M] A[M];                          // transition probabilities
  matrix<lower=0.0>[T,K] dir_shape[M];      // Parameters for Dirichlet distribution
  matrix<lower=0.0>[T,K] loc_parm[M];       // Overall shape parameter decomposed into
  matrix<lower=0.0>[M,T] scale_parm;        // location part and scale parts
  
  { // Transition probability setup          
    for (i in 1:(M-1) ) { 
      A[i] = rep_vector(0.0, M);            // Using a constrained matrix,
      A[i,i] = prob_remain[i];              // when you jump a state you cannot 
      A[i,i+1] = 1.0 - prob_remain[i];      // go back to the previous state
    }
    A[M] = rep_vector(0.0, M);
    A[M,M] = 1.0;
  }
  
  { // Observation likelihood

    // Get the Dirichlet shape parameters
 
    for(m in 1:M) {
      loc_parm[m] = exp(X*to_matrix(b_loc[m]) );
      scale_parm[m] = to_row_vector(exp(X*to_vector(b_scale[m])));
      for(t in 1:T) {
        loc_parm[m,t] = to_row_vector(normalize(to_vector(loc_parm[m,t]) ) );
        dir_shape[m,t] = scale_parm[m,t]*loc_parm[m,t];
      }
    }

    // Now the likelihood value
    for(t in 1:T) {
      for(m in 1:M) {
        logoblik[t,m] = dirichlet_lpdf(y[t] | to_vector(dir_shape[m,t]) );
      }
    }
  }

  { // Forward algorithm log p(z_t = j | x_{1:t})
    real accumulator[M];

   // For the first P terms, start in state 1
   // By forcing this, we essentially setup a 
   // "burn-in" period for our P regression
   // coefficients - This is appropriate
   // Given we have seasonal data - let the seasons
   // be seeded with the first observation
    // We essentially force the HMM to start in state 1
    // by adding (subtracting) a large penalty to higher states
    for(t in 1:P) {
      logalpha[t][1] = log(1) + logoblik[t][1];    // always start in stat 1
      for (j in 2:M) {                             // pi_1 = 1, pi_j = 0
        // The equation below should be logalpha[t][j] = -\intfy or negative_infinity(). 
        // Here, we use a large magnitude negative numbers so the gradient has a finite value
        logalpha[t][j] = -10^200;                  // on log scale, near negative infinity
      }
    }

    for (t in (P+1):T) {
      for (j in 1:M) { // j = current (t)
        for (i in 1:M) { // i = previous (t-1)
                         // Murphy (2012) Eq. 17.48
                         // belief state + transition prob + local evidence at t
          accumulator[i] = logalpha[t-1, i] + log(A[i,j]) + logoblik[t,j];
        }
        logalpha[t, j] = log_sum_exp(accumulator);
      }
    }
  } // Forward
}

model {

  to_array_1d(b_loc) ~ normal(0, b_loc_sig);
  to_array_1d(b_scale) ~ normal(0, b_scale_sig);

  prob_remain ~ beta(rho_a, rho_b);

  target += log_sum_exp(logalpha[T]); // Note: update based only on last logalpha
}

generated quantities {

  simplex[M] alpha[T];
  int<lower=1, upper=M> zstar[T];
  real logp_zstar;

  { // Forward algortihm
    for (t in 1:T)
      alpha[t] = softmax(logalpha[t]);
  } // Forward


  { // Viterbi algorithm
    int bpointer[T, M];             // backpointer to the most likely previous state on the most probable path
    real delta[T, M];               // max prob for the seq up to t
                                    // with final output from state k for time t

    for (j in 1:M)
      delta[1, j] = logoblik[1,j];

    for (t in 2:T) {
      for (j in 1:M) { // j = current (t)
        delta[t, j] = negative_infinity();
        for (i in 1:M) { // i = previous (t-1)
          real logp;
          logp = delta[t-1, i] + log(A[i, j]) + logoblik[t,j];
          if (logp > delta[t, j]) {
            bpointer[t, j] = i;
            delta[t, j] = logp;
          }
        }
      }
    }

    logp_zstar = max(delta[T]);

    for (j in 1:M)
      if (delta[T, j] == logp_zstar)
        zstar[T] = j;

    for (t in 1:(T - 1)) {
      zstar[T - t] = bpointer[T - t + 1, zstar[T - t + 1]];
    }
  }
}



