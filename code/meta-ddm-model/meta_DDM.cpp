#include <Rcpp.h>
// [[Rcpp::depends(RcppZiggurat)]]
#include <Ziggurat.h>
using namespace Rcpp;

static Ziggurat::Ziggurat::Ziggurat zigg;

// [[Rcpp::export]]
NumericVector meta_DDM(
    double v, double a, double ter, double z, int ntrials, double s, double dt, 
    NumericVector t2distribution, double conf_add, double conf_mult, double vratio,
    double v_bias){
  
  // initialize output
  NumericMatrix DATA(ntrials,6);
  double sqrt_dt = s* sqrt(dt);
  
  // loop over trials
  for (int i = 0; i < ntrials; i++) {
    // initalize variables
    int resp = -1;
    int acc = 0;
    double evidence = a*z;
    double t = 0;
    double t2time = 0;
    int stim = 2*R::unif_rand(); // 0 or 1
    stim = stim*2 - 1;
    
    double v_trial_dt = v*stim*dt; // change the direction of v depending on stimulus
    
    // Decisional processing
    while (evidence>=0 && evidence<=a){
      
      t = t + dt;
      evidence = evidence + v_trial_dt + sqrt_dt * zigg.norm();
    }
    
    if (evidence >= a){
      resp = 1;
      evidence=a;
      if (stim > 0){
        acc = 1;
      }
    } else if (evidence <= 0) {
      resp = -1;
      evidence=0;
      if (stim < 0){
        acc = 1;
      }
    }
    
    DATA(i,0) = (t + ter);
    DATA(i,1) = resp;
    DATA(i,2) = acc;
    DATA(i,3) = stim;
    
    t2time = t2distribution(i);
    
    //Post-decisional processing
    double v_post = v*stim*vratio + v_bias*resp;

    // for (int j = 0; j < round(t2time/dt); j++){
      // t2 = t2 + dt;
      // evidence = evidence + v_post_dt + sqrt_dt * zigg.norm();
    // }
    evidence = evidence + t2time*v_post + sqrt(t2time)*zigg.norm();
    
    DATA(i,4) = evidence;

    double conf_ev;
    if (resp == 1){
      conf_ev = evidence-a;
    }
    if(resp == -1){
      conf_ev = -evidence;
    }
    
    DATA(i,5) = 1/(1+exp(-(conf_ev-conf_add)*conf_mult)); // logistic transform to get ratings between 0 and 1
    
  }
  
  return DATA; //RT, resp,accuracy, evidence2, rt2, confidence
  
}

