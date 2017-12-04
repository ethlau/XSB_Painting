#include <stdio.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_dht.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include "gas_model.h"

//--added by SF:
/*
double kappa_integrant(double r, void* p){ // see Ma et al 2015 Eq 3.2
    double *params = (double *) p;
    double ell = params[0];
    double chi = params[1];
    double Rs = params[2]; //NFW scale radius in Mpc
    double rhoS = params[3]; // NFW density at the scale radius
    double rho_NFW = rhoS/( (r/Rs) * pow(1+r/Rs,2.0) );
    double integrant;
    if(ell*r>0){
        integrant = 4.0*M_PI*pow(r,2) * sin(ell*r/chi) / (ell*r/chi) * rho_NFW;
    }
    else if(ell*r==0) {
        integrant = 4.0*M_PI*pow(r,2) * rho_NFW;
    }
    return integrant;
}
*/

double kappa_integrant(double x, void* p){ // see Ma et al 2015 Eq 3.2
    double *params = (double *) p;
    double ell = params[0];
    double chi = params[1];
    double Rs = params[2]; //NFW scale radius in Mpc
    double rhoS = params[3]; // NFW density at the scale radius
    double Rvir = params[4];
    double rho_NFW = rhoS/( (x*Rvir/Rs) * pow(1+x*Rvir/Rs,2.0) );
    double integrant;
    if(ell*x>0){
        integrant = 4.0*M_PI*pow(Rvir,3)*pow(x,2) * sin(ell*Rvir*x/chi) / (ell*Rvir*x/chi) * rho_NFW;
    }
    else if(ell*x==0) {
        integrant = 4.0*M_PI*pow(Rvir,3)*pow(x,2) * rho_NFW;
    }
    return integrant;
}

double ttx_func(double x, void * p) {
  struct my_func_params * params = (struct my_func_params *)p;
  float n = params->b;
  float C = params->c;
  double beta = params->d;
  float delta_rel = params->a;
  int pturbrad = params->e;
  float delta_rel_n= params-> f;
  gas_model gmod(delta_rel, n, C, 0.0, 0.0, 0.0, pturbrad, delta_rel_n);
  double ff = pow(gmod.theta(x, beta), n)*pow(x,2);
  return ff;
}

double tx_func(double x, void * p) {
  struct my_func_params * params = (struct my_func_params *)p;
  float n = params->b;
  float C = params->c;
  double beta = params->d;
  float delta_rel = params->a;
  int pturbrad = params->e;
  float delta_rel_n= params-> f;
  gas_model gmod(delta_rel, n, C, 0.0, 0.0, 0.0, pturbrad, delta_rel_n);
  double ff = pow(gmod.theta(x, beta),n+1.0)*pow(x,2);
  return ff;
}

double tx_func_p(double x, void * p) {
  struct my_func_params * params = (struct my_func_params *)p;
  float n = params->b;
  float C = params->c;
  double beta = params->d;
  float delta_rel = params->a;
  int pturbrad = params->e;
  float delta_rel_n= params-> f;
  gas_model gmod(delta_rel, n, C, 0.0, 0.0, 0.0, pturbrad, delta_rel_n);
  double ff = delta_rel*pow(gmod.theta(x,beta),n-1.0)*pow(x,2);
  return ff;
}


double ftx_func(double x, void * p) {
  struct my_func_params * params = (struct my_func_params *)p;
  float n = params->b;
  float C = params->c;
  double beta = params->d;
  float delta_rel = params->a;
  int pturbrad = params->e;
  float delta_rel_n= params-> f;
  gas_model gmod(delta_rel, n, C, 0.0, 0.0, 0.0, pturbrad, delta_rel_n);
  double ff = gmod.f(x)*pow(gmod.theta(x, beta),gmod.n)*pow(x,2);
  return ff;
}

double sx_func (double x, void * params) {
    double f = (x - (1.0+x)*log(1.0+x)) / (pow(x,3)*pow(1.0+x,3));
    return f;
}

double ss_func(double x, void * params) {
  double C = *(double *) params;
  gas_model gmod(0.0, 0.0, C, 0.0, 0.0, 0.0, 0, 0.0);
  double f = gmod.S_cx(x)*pow(x,2);
  return f;
}

double fx_func(double x, void * params) {
  double C = *(double *) params;
  gas_model gmod(0.0, 0.0, C, 0.0, 0.0, 0.0, 0, 0.0);
  double g = gmod.f(x)*x/pow(1.0+x,2);
  return g;
}

double yproj_func(double x, void * p) {
  double *params = (double *) p;
  float n = (float)params[1];
  float C = (float)params[2];
  float delta_rel = (float)params[0];
  double beta = params[3];
  double R = params[4];
  int pturbrad = (int)params[5];
  float delta_rel_n = params[6];
  float R500_div_ri = params[7];
  gas_model gmod(delta_rel, n, C, 0.0, 0.0, 0.0, pturbrad, delta_rel_n);
  double y;
  if ((x-R)<0.01) y = 0.0;
  else  y = pow(gmod.theta(x,beta),gmod.n+1.0)*x/pow(x*x - R*R,0.5);
  if (pturbrad==2) y = y*(1.0 - delta_rel*pow(x/R500_div_ri, delta_rel_n));
  return y;
}

double yfft_func(double x, void * p) {
    double *params = (double *) p;
    double n = (double)params[1];
    double C = (double)params[2];
    double delta_rel = (double)params[0];
    double beta = params[3];
    double ell = params[4];
    double elli = params[5];
    int pturbrad = (int)params[6];
    float delta_rel_n = params[7];
    float R500_div_ri = params[8];
    gas_model gmod(delta_rel, n, C, 0.0, 0.0, 0.0, pturbrad, delta_rel_n);
    double y;
    if(ell*x>0){
        y = pow(gmod.theta(x,beta),gmod.n+1.0)*x*x*sin(ell*x/elli)/(ell*x/elli);
    }
    else if(ell*x==0){
        y = pow(gmod.theta(x,beta),gmod.n+1.0)*x*x;
    }
    if (pturbrad==2) y = y*(1.0 - delta_rel*pow(x/R500_div_ri, delta_rel_n));
    return y;
    }



double arnaud_func(double x, void * p) {
  double *params = (double *) p;
  double c500 = (double)params[0];
  double alpha = (double)params[1];
  double beta = (double)params[2];
  double gamma = (double)params[3];
  double ell = (double)params[4];
  double elli = (double)params[5];
  double y;
  y = (1.0/ (pow(x,gamma)*pow(1.0 + pow(x,alpha),(beta-gamma)/alpha)))*x*x*sin(ell*x/elli)/(ell*x/elli);
  return y;
}

double yfft_func_k(double x, void * p) {
  double *params = (double *) p;
  double n = (double)params[1];
  double C = (double)params[2];
  double delta_rel = (double)params[0];
  double beta = params[3];
  double k = params[4];
  //double elli = params[5];
  int pturbrad = (int)params[5];
  float delta_rel_n = params[6];
  float R500_div_ri = params[7];
  gas_model gmod(delta_rel, n, C, 0.0, 0.0, 0.0, pturbrad, delta_rel_n);
  double y;
  y = pow(gmod.theta(x,beta),gmod.n+1.0)*x*x*sin(k*x)/(k*x);
  if (pturbrad==2) y = y*(1.0 - delta_rel*pow(x/R500_div_ri, delta_rel_n));
  return y;
}

double arnaud_func_k(double x, void * p) {
  double *params = (double *) p;
  double c500 = (double)params[0];
  double alpha = (double)params[1];
  double beta = (double)params[2];
  double gamma = (double)params[3];
  double k = (double)params[4];
  //double elli = (double)params[5];
  double y;
  y = (1.0/ (pow(x,gamma)*pow(1.0 + pow(x,alpha),(beta-gamma)/alpha)))*x*x*sin(k*x)/(k*x);
  return y;
}

double proj_arnaud_func(double x, void * p) {
  double *params = (double *) p;
  double c500 = (double)params[0];
  double alpha = (double)params[1];
  double beta = (double)params[2];
  double gamma = (double)params[3];
  double R = (double)params[4];
  double yint;
  if (sqrt(x*x-R*R)<0.00001) yint = 0.0;
  else yint = (1.0 / (pow(x,gamma)*pow(1.0 + pow(x,alpha),(beta-gamma)/alpha)))*x / pow(x*x-R*R,0.5);
  return yint;
}

double proj_KS02_func(double x, void * p) {
  double *params = (double *) p;
  double C = (double)params[0];
  double gamma;// = (double)params[1];
  double B;// = (double)params[2];
   double R = (double)params[3];
   double yint, ygas, eta;

  if (sqrt(x*x-R*R)<0.00001) yint = 0.0;
  else{
    gamma = 1.137 + 8.94*pow(10,-2)*log(C/5) - 3.68*pow(10,-3)*(C - 5);
    eta = 2.235 + 0.202*(C - 5)-1.16*pow(10,-3)*pow(C - 5,2);
    B = 3*pow(eta,-1)*(gamma - 1)/gamma*pow((log(1 + C)/C - 1/(1 + C)),-1);
    ygas = pow((1 - B*(1 - log(1 + x)/x)),1/(gamma - 1));
    yint = pow(ygas,gamma)*x / pow(x*x-R*R,0.5);}
  return yint;
}

double yint_func(double x, void * p) {
  // returns Pgas/P0 as a function of x=r/Rs
  double *params = (double *) p;
  float n = (float)params[1];
  float C = (float)params[2];
  float delta_rel = (float)params[0];
  double beta = params[3];
  double R = params[4];
  int pturbrad = (int)params[5];
  float delta_rel_n = params[6];
  float R500_div_ri = params[7];
  gas_model gmod(delta_rel, n, C, 0.0, 0.0, 0.0, pturbrad, delta_rel_n);
  double y;
  y = pow(gmod.theta(x,beta),gmod.n+1.0)*x*x;
  if (pturbrad==2) y = y*(1.0 - delta_rel*pow(x/R500_div_ri, delta_rel_n));
  return y;
}

double yint_func_sam(double x, void * p) {
  // returns Pgas/P0 as a function of x=r/R500 !!
  double *params = (double *) p;
  float n = (float)params[1];
  float C = (float)params[2];
  float delta_rel = (float)params[0];
  double beta = params[3];
  double R = params[4];
  int pturbrad = (int)params[5];
  float delta_rel_n = params[6];
  float R500_div_ri = params[7];
  gas_model gmod(delta_rel, n, C, 0.0, 0.0, 0.0, pturbrad, delta_rel_n);
  double y = pow(gmod.theta(x * R500_div_ri ,beta),gmod.n+1.0)*x*x;
  if (pturbrad==2) y *= (1.0 - delta_rel*pow(x, delta_rel_n));
  return y;
}


double mgas500_func(double x, void * p) {
  double *params = (double *) p;
  float n = (float)params[1];
  float C = (float)params[2];
  float delta_rel = (float)params[0];
  //float f_s = (float)params[4];
  double beta = params[3];
  //double R = params[4];
  int pturbrad = (int)(params[5]);
  float delta_rel_n = params[6];
  float R500_div_ri = params[7];

  gas_model gmod(0.0, n, C, 0.0, 0.0, 0.0, pturbrad, 0.0);
  double y;
  y = pow(gmod.theta(x,beta),gmod.n)*x*x;
  //if (pturbrad==2) y = y*(1.0 - delta_rel*pow(x/R500_div_ri, delta_rel_n));
  return y;
}

double mgas500_func_mod(double x, void * p) {
    // this is using the broken power-law model
    double *params = (double *) p;
    float n = (float)params[1];
    float C = (float)params[2];
    float delta_rel = (float)params[0];
    //float f_s = (float)params[4];
    double beta = params[3];
    //double R = params[4];
    int pturbrad = (int)(params[5]);
    float x_break_div_ri = params[6]; // breaking point in units of NFW scale radius
    float npoly_prime = params[7];

    gas_model gmod(0.0, n, C, 0.0, 0.0, 0.0, pturbrad, 0.0);
    double y;
    if(x>x_break_div_ri){
        y = pow(gmod.theta(x,beta),gmod.n)*x*x;
    }
    else if(x<=x_break_div_ri){
        y = pow(gmod.theta_mod(x,beta,x_break_div_ri,npoly_prime),npoly_prime)*x*x;
        y *= pow(gmod.theta(x_break_div_ri,beta),gmod.n) / pow(gmod.theta_mod(x_break_div_ri,beta,x_break_div_ri,npoly_prime),npoly_prime); // normalize at break
    }

  return y;
}



double gasproj_func(double x, void * p) {
  // SF: this gives the integrated gas density (it's theta^n; the pressure goes like theta^(n+1))
  // multiplied with some factor that comes from substitution:
  // int f(r) dl with dl = dl/dr dr and l=sqrt(r^2-R^2) this gives
  // = int f(r) r/(r^2-R^2) dr.
  double *params = (double *) p;
  float n = (float)params[1];
  float C = (float)params[2];
  float delta_rel = (float)params[0];
  //float f_s = (float)params[4];
  double beta = params[3];
  double R = params[4];
  int pturbrad = (int)(params[5]);
  float delta_rel_n = params[6];
  float R500_div_ri = params[7];

  gas_model gmod(0.0, n, C, 0.0, 0.0, 0.0, pturbrad, 0.0);
  double y;
  if ((x-R)<0.01) y = 0.0;
  else  y = pow(gmod.theta(x,beta),gmod.n)*x/pow(x*x - R*R,0.5);
  //if (pturbrad==2) y = y*(1.0 - delta_rel*pow(x/R500_div_ri, delta_rel_n));
  // suman commented above line out
  // my interpretation: the non-thermal component only affects the pressure, not the density.
  return y;
}

double gasproj_func_mod(double x, void * p) {
  // SF: this gives the integrated gas density (it's theta^n; the pressure goes like theta^(n+1))
  // multiplied with some factor that comes from substitution:
  // int f(r) dl with dl = dl/dr dr and l=sqrt(r^2-R^2) this gives
  // = int f(r) r/(r^2-R^2) dr.
  double *params = (double *) p;
  float n = (float)params[1];
  float C = (float)params[2];
  float delta_rel = (float)params[0];
  //float f_s = (float)params[4];
  double beta = params[3];
  double R = params[4];
  int pturbrad = (int)(params[5]);
  float delta_rel_n = params[6];
  float R500_div_ri = params[7];
  float x_break = params[8]; // x_break in units of R500
  x_break *= R500_div_ri; //now x_break in units of Rs
  float npoly_break = params[9];

  gas_model gmod(0.0, n, C, 0.0, 0.0, 0.0, pturbrad, 0.0);
  double y;
  if ((x-R)<0.01) {y = 0.0;}
  else if( (x-R)>=0.01 and x>x_break){
    y = pow(gmod.theta(x,beta),gmod.n)*x/pow(x*x - R*R,0.5);
  }
  else if( (x-R)>=0.01 and x<x_break){
    y = pow(gmod.theta_mod(x,beta,x_break,npoly_break),npoly_break)*x/pow(x*x - R*R,0.5);
    y *= pow(gmod.theta(x_break,beta),gmod.n) / pow(gmod.theta_mod(x_break,beta,x_break,npoly_break),npoly_break); //normalize at break
  }
  //if (pturbrad==2) y = y*(1.0 - delta_rel*pow(x/R500_div_ri, delta_rel_n));
  // suman commented above line out
  // my interpretation: the non-thermal component only affects the pressure, not the density.
  return y;
}

// SF: I added this for computing the proj. NFW profile
// this function returns rho_Nfw(x) * x/(x^2 + R^2)
double NFWproj_func(double x, void *p) {

  double *params = (double *) p;
  double R = params[0];
  double y = 1.0/(x*(1.0+x)*(1.0+x)); // that's the NFW profile without normalization
  y *= x/pow(x*x - R*R,0.5); // that factor comes from parameter substitution

  return y;
}



// function for solving model goes here

double gasmod_apply_bc(const gsl_vector * v, void *p) {
  double *params = (double *) p;
  float n = (float)params[1];
  float C = (float)params[2];
  float delta_rel = (float)params[0];
  double Aprime = params[4];
  double Bprime = params[5];
  float f_s = (float)params[6];
  int pturbrad = (int)params[7];
  float delta_rel_n= (float) params[8];
  gas_model gmod(delta_rel, n, C, Aprime, Bprime, f_s, pturbrad, delta_rel_n);
  double x0 = gsl_vector_get (v, 0); // beta
  double x1 = gsl_vector_get (v, 1); // Cf
  if (x0<0.01) x0 = 0.01;
  //if (x0>14) x0 = 14;
  if (x1<0.01) x1 = 0.01;
  //if (x1>14) x1 = 14;
  return sqrt(pow(gmod.energy_constraint(x0, x1),2)+pow(gmod.pressure_constraint(x0, x1),2));
}
