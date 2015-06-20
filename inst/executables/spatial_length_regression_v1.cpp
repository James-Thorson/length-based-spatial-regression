#include <TMB.hpp>

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// function for logistic transform
template<class Type>
Type plogis(Type x){
  return 1.0 / (1.0 + exp(-x));
}

// dlognorm
template<class Type>
Type dlognorm(Type x, Type meanlog, Type sdlog, int give_log=false){
  Type Return;
  if(give_log==false) Return = dnorm( log(x), meanlog, sdlog, false) / x;
  if(give_log==true) Return = dnorm( log(x), meanlog, sdlog, true) - log(x);
  return Return;
}

// dzinflognorm
template<class Type>
Type dzinflognorm(Type x, Type meanlog, Type encounter_prob, Type log_notencounter_prob, Type sdlog, int give_log=false){
  Type Return;
  if(x==0){
    if(give_log==false) Return = 1.0 - encounter_prob;
    if(give_log==true){
      if( isNA(log_notencounter_prob) ) Return = log(1.0 - encounter_prob);
      if( !isNA(log_notencounter_prob) ) Return = log_notencounter_prob;
    }
  }else{
    if(give_log==false) Return = encounter_prob * dlognorm( x, meanlog, sdlog, false );
    if(give_log==true) Return = log(encounter_prob) + dlognorm( x, meanlog, sdlog, true );
  } 
  return Return;
}

// dzinfgamma, shape = 1/CV^2, scale = mean*CV^2
template<class Type>
Type dzinfgamma(Type x, Type posmean, Type encounter_prob, Type log_notencounter_prob, Type cv, int give_log=false){
  Type Return;
  if(x==0){
    if(give_log==false) Return = 1.0 - encounter_prob;
    if(give_log==true){
      if( isNA(log_notencounter_prob) ) Return = log(1.0 - encounter_prob);
      if( !isNA(log_notencounter_prob) ) Return = log_notencounter_prob;
    }
  }else{
    if(give_log==false) Return = encounter_prob * dgamma( x, pow(cv,-2), posmean*pow(cv,2), false );
    if(give_log==true) Return = log(encounter_prob) + dgamma( x, pow(cv,-2), posmean*pow(cv,2), true );
  } 
  return Return;
}

// dzinfnorm
template<class Type>
Type dzinfnorm(Type x, Type posmean, Type encounter_prob, Type log_notencounter_prob, Type cv, int give_log=false){
  Type Return;
  if(x==0){
    if(give_log==false) Return = 1.0 - encounter_prob;
    if(give_log==true){
      if( isNA(log_notencounter_prob) ) Return = log(1.0 - encounter_prob);
      if( !isNA(log_notencounter_prob) ) Return = log_notencounter_prob;
    }
  }else{
    if(give_log==false) Return = encounter_prob * dnorm( x, posmean, posmean*cv, false );
    if(give_log==true) Return = log(encounter_prob) + dnorm( x, posmean, posmean*cv, true );
  } 
  return Return;
}

// Main function
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace R_inla;
  using namespace Eigen;
  using namespace density;
  
  // Options
  DATA_FACTOR( Options_vec );
  // Slot 0:  0=no Nu, 1=yes Nu
  // Slot 1: 0=no Omega, 1=yes Omega
  // Slot 2: 0=no Delta, 1=yes Delta
  // Slot 3: encounter function, 2=negative-exponential
  // Slot 4: Distribution for positive catches, 1=Gamma

  // Declared dimensions
  DATA_INTEGER(n_s);      // Number of sites
  DATA_INTEGER(n_b);      // Number of length bins

  // Data vectors
  DATA_VECTOR(c_i);         // Expanded counts for observation i 
  DATA_MATRIX(x_ij);         // Covariates for observation i 
  DATA_FACTOR(s_i);          // site for observation i
  DATA_FACTOR(b_i);          // Length bin for observation i
  
  // Extract dimensions
  int n_i = c_i.size();   // number of observations
  int n_j = x_ij.row(0).size();   // number of covariates

  // Aniso objects
  DATA_STRUCT(spde,spde_aniso_t);
  
  // Fixed effects
  PARAMETER_VECTOR(ln_H_input); // Anisotropy parameters
  PARAMETER_VECTOR(gamma_j);        // Covariate effect (presumably includes negative-exponential effect of length)
  PARAMETER_VECTOR(logeta_vec);         // 3 Slots: 0=Nu, 1=Omega, 2=Delta  (using vector to allow mirroring)
  PARAMETER_VECTOR(rho_vec);     // 2 Slots: 0=Nu, 1=Omega (using vector to allow mirroring)
  PARAMETER_VECTOR(logkappa_vec);   // 2 Slots: 0=Omega, 1=Delta  (using vector to allow mirroring)
  PARAMETER_VECTOR(theta_vec);   // 3 Slots: 0=encounter-prob-slope, 1=encounter-prob-asymptote, 1=CV-of-residual-errors
  
  // Random effects 
  PARAMETER_VECTOR(Nuinput_b);      // Variation by length-bin
  PARAMETER_ARRAY(Omegainput_sb);   // Interactive variation by site and length-bin  (site x length-bin)
  PARAMETER_VECTOR(Deltainput_s);     // Spatial variation by site
  
  // Objective function components
  Type jnll = 0;
  Type pi = 3.141592;  
  vector<Type> jnll_c(4); // Slots: 0=Nu, 1=Omega, 2=Delta, 3=data 
  jnll_c.setZero();
  vector<Type> nll_i(n_i);
  nll_i.setZero();
  
  // Random effect derived quantities
  vector<Type> logtau_vec( 3 );
  vector<Type> MargSigma_vec( 3 );
  vector<Type> Range_vec( 3 );
  // 0. Nu -- marginal length variation
  logtau_vec(0) = NA_REAL;
  MargSigma_vec(0) = 1 / exp(logeta_vec(0));    
  Range_vec(0) = NA_REAL;
  // 1. Omega -- length variation by site
  // 2. Delta -- marginal variation by site
  for(int z=0; z<2; z++){
    logtau_vec(z+1) = logeta_vec(z+1) - logkappa_vec(z);
    MargSigma_vec(z+1) = 1 / sqrt(4*pi*exp(2*logtau_vec(z+1))*exp(2*logkappa_vec(z)));    
    Range_vec(z+1) = sqrt(8) / exp( logkappa_vec(z) );
  }
  
  // Derived random effects
  vector<Type> Nu_b(n_b);
  vector<Type> Delta_s(n_s);
  matrix<Type> Omega_sb(n_s,n_b);
  matrix<Type> log_Lambda_sb(n_s,n_b);
  Nu_b = Nuinput_b * MargSigma_vec(0);
  vector<Type> eta_i(n_i);
  eta_i = x_ij * gamma_j.matrix();
  for(int s=0; s<n_s; s++){
    Delta_s(s) = Deltainput_s(s) / exp(logtau_vec(2));
    for(int b=0; b<n_b; b++){
      Omega_sb(s,b) = Omegainput_sb(s,b) / exp(logtau_vec(1));
      log_Lambda_sb(s,b) = Nu_b(b) + Delta_s(s) + Omega_sb(s,b);
    }
  }
  
  // Anisotropy elements
  matrix<Type> H(2,2);
  H(0,0) = exp(ln_H_input(0));
  H(1,0) = ln_H_input(1);
  H(0,1) = ln_H_input(1);
  H(1,1) = (1+ln_H_input(1)*ln_H_input(1)) / exp(ln_H_input(0));
  Type H_trace = H(0,0)+H(1,1);
  Type H_det = H(0,0)*H(1,1)-H(0,1)*H(1,0);

  // Calculate adjugate of H
  matrix<Type> adj_H(2,2);
  adj_H(0,0) = H(1,1);
  adj_H(0,1) = -1 * H(0,1);
  adj_H(1,0) = -1 * H(1,0);
  adj_H(1,1) = H(0,0);
  
  // Random field probability                                                                                                                              
  Eigen::SparseMatrix<Type> Q_omega;
  Q_omega = Q_spde(spde,exp(logkappa_vec(0)),H);
  Eigen::SparseMatrix<Type> Q_delta;
  Q_delta = Q_spde(spde,exp(logkappa_vec(1)),H);
  if(Options_vec(0)==1) jnll_c(0) += AR1(rho_vec(0))(Nuinput_b);
  if(Options_vec(1)==1) jnll_c(1) += SEPARABLE(AR1(rho_vec(1)),GMRF(Q_omega))(Omegainput_sb);  // site-by-year => Seperable(AR1, GMRF)
  if(Options_vec(2)==1) jnll_c(2) += GMRF(Q_delta)(Deltainput_s);

  // Probability of data
  // Likelihood
  vector<Type> logchat_i(n_i);
  vector<Type> encounterprob_i(n_i);
  Type log_notencounterprob = NA_REAL;
  for(int i=0; i<n_i; i++){
    // chat_i is the expectation for compound process including both encounter probability and positive catch rate components
    logchat_i(i) = eta_i(i) + log_Lambda_sb(s_i(i),b_i(i)); // eta_i(i) + Nu_b(b_i(i)) + Delta_s(s_i(i)) + Omega_sb(s_i(i),b_i(i));
    // Likelihood
    if( !isNA(c_i(i)) ){                
      // Calculate encounter probability (only used for delta-lognormal model)
      if( Options_vec(3)==0 ) encounterprob_i(i) = plogis( theta_vec(1) + theta_vec(0)*logchat_i(i) );
      if( Options_vec(3)==1 ) encounterprob_i(i) = plogis( theta_vec(1) ) * ( 1.0 - exp(-1 * exp(logchat_i(i)) * exp(theta_vec(0))) );
      if( Options_vec(3)==2 ){
        encounterprob_i(i) = ( 1.0 - exp(-1 * exp(logchat_i(i)) * exp(theta_vec(0))) );
        log_notencounterprob = -1 * exp(logchat_i(i)) * exp(theta_vec(0));
      }
      // probability of data
      if( Options_vec(4)==0 ) nll_i(i) = -1 * dzinflognorm(c_i(i), logchat_i(i)-log(encounterprob_i(i)), encounterprob_i(i), log_notencounterprob, exp(theta_vec(2)), true);
      if( Options_vec(4)==1 ) nll_i(i) = -1 * dzinfgamma(c_i(i), exp(logchat_i(i))/encounterprob_i(i), encounterprob_i(i), log_notencounterprob, exp(theta_vec(2)), true);
      if( Options_vec(4)==2 ) nll_i(i) = -1 * dzinfnorm(c_i(i), exp(logchat_i(i))/encounterprob_i(i), encounterprob_i(i), log_notencounterprob, exp(logchat_i(i))/encounterprob_i(i)*exp(theta_vec(2)), true);
    }
  }
  
  // Combine likelihood
  jnll_c(3) = nll_i.sum();
  jnll = jnll_c.sum();
  
  // Reporting
  REPORT( logchat_i );
  REPORT( encounterprob_i );
  REPORT( H );
  REPORT( nll_i );
  REPORT( jnll_c );
  REPORT( jnll );
  REPORT( eta_i );
  
  // 
  REPORT( Nu_b );
  REPORT( Omega_sb );
  REPORT( Delta_s );
  REPORT( log_Lambda_sb );
  REPORT( Range_vec );
  REPORT( logtau_vec );
  REPORT( MargSigma_vec );
  
  return jnll;
  
}
