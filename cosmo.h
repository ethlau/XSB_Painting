// defines the cosmology and allows calculations of various cosmological params
#ifndef _COSMO_
#define _COSMO_

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>

//double angdiam_func(double x, void * p);


using namespace std;

class cosmo {

 protected:
  float H0, Omega_M, Omega_b, rho_crit, a, Omega_k, Omega_L, wt;
  float PI, m_sun, G, mpc, gnewt;

 public:

  cosmo(float inp1, float inp2, float inp3, float inp4, float inp5)
 {
    // cosmological params at z =0
    H0 = inp1;
    Omega_M = inp2;
    Omega_b = inp3;
    Omega_k = inp4;
    Omega_L = 1.0- Omega_M- Omega_k;
    wt= inp5;
    PI = 4.*atan(1.);
    m_sun = 1.98892e30;
    G = 6.67e-11/1.0e9; // in km^3 kg^-1 s^-2
    mpc = 3.0857e19; // in km
    gnewt = G*m_sun/mpc; //for cosm params (in km^2 Mpc msun^-1 s^-2)
  }

  float get_H0() {
    return H0;
  }

  float get_hubble70() {
    return H0/70.0;
  }

  float scale_fact(float redshift) {
    a = 1.0/(1.0+redshift);
    return a;
  }

  float hubblez(float redshift) {
    scale_fact(redshift);
    return H0 * Efact(redshift);
  }

  float calc_rho_crit(float redshift) {
    float Hz = hubblez(redshift);
    rho_crit = pow(Hz,2)*3./(8.*PI*gnewt);
    return rho_crit;
  }

  float Omega_Mz(float redshift) {
    scale_fact(redshift);
    float Hz = hubblez(redshift);
    return (Omega_M/pow(a,3))*pow(H0/Hz,2);
  }

  float Delta_vir(float redshift) {
    float x;
    x  = Omega_Mz(redshift) - 1.0;
    return (18*pow(PI,2)) + (82*x) - (39*pow(x,2));
 }

  float Efact(float redshift) {
    scale_fact(redshift);
    return sqrt(Omega_M/pow(a,3) + Omega_L/pow(a, 3*wt+3) + (1.0-Omega_M-Omega_L)/pow(a,2));
  }

  float ang_diam(float redshift) {
    // note rcutoff is in units of Rvir
    int dummy;
    float result;
    double x[1000], y[1000];
    double units = scale_fact(redshift)*3.0e5/H0;
   
    for(dummy=0;dummy<1000;dummy++){
      x[dummy]= 0.0+ redshift/999*dummy;
      y[dummy] = 1.0/(sqrt( Omega_M*pow(1+x[dummy],3) + Omega_L*pow(1+x[dummy], 3*wt+3) + (1.0-Omega_M-Omega_L)*pow(1+x[dummy],2)));
    }
    if (redshift==0.0) return 0.0;
    else {
        gsl_interp_accel *acc = gsl_interp_accel_alloc ();
 
        gsl_spline *spline 
            = gsl_spline_alloc (gsl_interp_cspline, 1000);

        gsl_spline_init (spline, x, y, 1000);
        result=units*gsl_spline_eval_integ (spline, x[0], x[dummy-1], acc); 
        gsl_spline_free (spline);
        gsl_interp_accel_free (acc);
    }
    return result;
 }
  float proper_dist(float redshift) {
      return ang_diam(redshift) * (1.0+redshift);
  }

  float lum_dist(float redshift) {
      return ang_diam(redshift) * (1.0+redshift)*(1.0+redshift);
  }
  float cosmic_time(float redshift) {
    // function to calculate cosmic time to redshift given cosmological model
    //good to ~1% for z<10

    // these are weight functions for solving the differential equation
    float x[16] = {0.0950,0.2816,0.4580,0.6179,0.7554,0.8656,0.9446,0.9894,-0.0950,-0.2816,-0.4580,-0.6179,-0.7554,-0.8656,-0.9446,-0.9894};
    float w[16] = {0.1895,0.1826,0.1692,0.1496,0.1246,0.0952,0.0623,0.0272,0.1895,0.1826,0.1692,0.1496,0.1246,0.0952,0.0623,0.0272};
   
    float a0 = 0.0, a1, Delta_a, h2, sum, cosm_h0, den0 = Omega_M, den1 = Omega_L, den2;
    float ss, ak, adot, cosm_t;
    float time_units, secperyr = 3.15569260e7, age_of_univ  = 8.93120195E+17/2.0; // in seconds
    int Nsteps, i, k;
    Nsteps = 10000;
    a1 = 1.0/(1.0+redshift);
    Delta_a = (a1-a0)/Nsteps;
    h2 = Delta_a/2.0;
    sum = 0.0;
    cosm_h0 = sqrt(8.0*PI/3.0);
    
    time_units = (age_of_univ*2/secperyr/1.0E9)*(1.0E2/H0); // converts to Gyr
    // now solve friedman eqn.
    for (i=0;i<Nsteps;i++) {
      ss = 0.0;
      for (k=0;k<16;k++) {
        ak = a0+h2*(1.0+x[k]);
        den2 = den1*dynrho(ak,wt);
        adot = sqrt( den0/ak + den2*ak*ak + 1.0-den0-den1 );
        ss = ss + w[k]/adot;
      }
      ss *= h2;
      a0 += Delta_a;
      sum += ss;
    }
    cosm_t = sum/cosm_h0;
    cosm_t = cosm_t*time_units;
    return cosm_t;
  }

  float dynrho(float a, float wdyn) {
    //       CALDWELL: This function RETURNs the energy density
    //       in the dynamical phi field, given the scale factor
    //       such that at a=1, dynrho = 1
    int n_dyn_q1,n_dyn_q2,n_dyn_q3;
    float drho;
    if (wdyn==-1.0) {
      n_dyn_q1 = 0;
      n_dyn_q2 = 0;
      n_dyn_q3 = 0;
    }
    else {
      n_dyn_q1 = 1;
      n_dyn_q2 = 1;
      n_dyn_q3 = 0;
    }
    //  did I ask for dynamical field? Leave this as an option for later
    if (n_dyn_q1==1) {
      if(n_dyn_q2==1) drho = pow(a,(-3*(1.0 + wdyn)));
    }
    //  otherwise, set the energy density to 1
    else  drho = 1.0;
    return drho;
  }
  





  friend class cluster;
  friend class growth; 
  //friend double angdiam_func(double x, void * p);
};


/*double angdiam_func(double x, void * p) {
  double *params = (double *) p;
  double OmegaM = params[0];
  double OmegaL = params[1];
  double wt= params[2];
  double y, a = 1.0/(1.0+x);
  y = 1.0/(sqrt( OmegaM/pow(a,3) + OmegaL/pow(a, 3*wt+3) + (1.0-OmegaM-OmegaL)/pow(a,2)));
  return y;
  }*/
#endif
