// classes defines a cluster object (mass, conc, redshift, etc).
#ifndef _CLUSTER_
#define _CLUSTER_

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
using namespace std;

class cluster {
 protected:
  float mass, radius, conc, redshift, rho_crit, ri, rhoi, B, C;
  float overden_id, relation, mass_overden, hubble;
  float PI, m_sun, G, mpc, gnewt;
  bool concset;
  //int relation;

 public:

  cluster(float inp1, float inp2, float inp3, float inp4, cosmo & cosm_model) {
    mass = inp1;
    redshift = inp2;
    overden_id = inp3;
    relation= inp4;
    conc = 0.0;
    concset = false;
    PI = 4.*atan(1.);
    m_sun = 1.98892e30;
    G = 6.67e-11/1.0e9; // in km^3 kg^-1 s^-2
    mpc = 3.0857e19; // in km
    gnewt = G*m_sun/mpc; //for cosm params (in km^2 Mpc msun^-1 s^-2)
    if (overden_id==-1) mass_overden =  cosm_model.Delta_vir(redshift);
	else if (overden_id==-3) mass_overden = 200.0;
	else if (overden_id==-4) mass_overden =cosm_model.Delta_vir(redshift);
    else mass_overden = overden_id;
    rho_crit = cosm_model.calc_rho_crit(redshift);
    hubble = cosm_model.get_H0()/100.0;
    calc_radius();
  }

  void reset_cluster(float inp1, float inp2, float inp3, float inp4, cosmo & cosm_model) {
    mass = inp1;
    redshift = inp2;
    overden_id = inp3;
    relation= inp4;
    conc = 0.0;
    concset = false;
    if (overden_id==-1) mass_overden =  cosm_model.Delta_vir(redshift);
    else mass_overden = overden_id;
    rho_crit = cosm_model.calc_rho_crit(redshift);
    hubble = cosm_model.get_H0()/100.0;
    calc_radius();
  }

  void calc_radius() {
    radius = mass/((4.0/3.0)*PI*mass_overden*rho_crit);
    radius = pow(radius, 0.333333);
  }

  float get_radius() {
    return radius;
  }

  float get_overden(){
    return mass_overden;
  }

  float get_ri() {
    if (!concset) {
      cout << "calculate concentration first!" << endl;
      return -1.0;
    }
    else {
    ri = radius/conc;
    return ri;
    }
  }

  float get_rhoi() {
    // must run concentration first
    float g;
    if (!concset) {
      cout << "calculate concentration first!" << endl;
      return -1.0;
    }
    else {
      get_ri();
      g = log(1.0+conc) - conc/(1.0+conc);
      rhoi = mass/(4.0*PI*pow(ri,3)*g); //in Msun/mpc^3
      return rhoi;
    }
  }

  float set_conc(float c) {
    conc = c;
    concset = true;
    get_ri();
    get_rhoi();
  }
  float bias(){
   if (overden_id == 200) {// Delta = 200


   }
  }

  double get_NFW_density(double r){
      double rho, rhos, x;
      if (!concset) {
          cout << "calculate concentration first!" << endl;
          return -1.0;
      }
      x = r/ri;
      rhos = ((double) rhoi)*(m_sun*1e-12)/mpc/mpc/mpc;
      rho = rhos/(x*pow(1.0+x, 2.0)); // total density in g/cm^3
      return rho;
  }

  float concentration(float conc_norm, float conc_mass_norm) {
    // uses mass, redshift and overdensity as initialised for the cluster
  float A, Mp;
    if (overden_id == 200) {// Delta = 200
      if (relation==1) {
	// mass-concentration relation for c = R200/ri, z = 0
	A = 5.74;
	B = -0.097;
	C = 0;
	Mp = 2e12; // Msol/h
      }
      else if (relation==2) {
	// relaxed params from Duffy et al. 08, z = 0
	A = 6.67;
	B = -0.092;
	C = 0;
	Mp = 2e12; // Msol/h
      }
      else if (relation==3) {
	// redshift varying 0->2, all
	A = 5.84;
	B = -0.084;
	C = -0.47;
	Mp = 2e12; // Msol/h
      }
      else if (relation==4) {
	// redshift varying, relaxed
	A = 6.71;
	B = -0.091;
	C = -0.44;
	Mp = 2e12; // Msol/h
      }
    }
    else if (overden_id == -1) {
    if (relation==1) {
      // Delta_vir relations
      // %%% mass-concentration relation for c = R200/ri, z = 0
      A = 7.96;
      B = -0.091;
      C = 0;
      Mp = 2e12; // Msol/h
    }
    else if (relation==2) {
      // relaxed params from Duffy et al. 08, z = 0
      A = 9.23;
      B = -0.089;
      C = 0;
      Mp = 2e12; // Msol/h
    }
    else if (relation==3) {
      // redshift varying 0->2, all
      A = 7.85;
      B = -0.081;
      C = -0.69;
      Mp = 2e12; // Msol/h
    }
    else if (relation==4) {
      // redshift varying 0->2, relaxed
      A = 9.23;
      B = -0.090;
      C = -0.44;
      Mp = 2e12; // Msol/h
    }
  }
  else if (overden_id==-2) {
      A = 4;
      Mp = 6.25e13;
      B = -0.13;
      C = 0;
   }
  else if(overden_id==-3){
	  conc= relation;
  }
  else if (overden_id==-4){
	  A=7.7;
	  B=-0.29;
	  C=-0.9;
	  Mp=5e13;
  }
   A *= conc_norm; // additional parameter allowing you to modify normalisation of mass-conc relation
   B *= conc_mass_norm; // modify c-M slope
   if(overden_id >= -3) conc =   A*pow((mass*hubble)/Mp,B)*pow(1.0+redshift,C);// calculate concentration
   if(overden_id==-4)  conc=pow(1+redshift,C*0.81)*A*pow(1.12*pow(mass*hubble/Mp,0.3)+0.53,B);
   //printf("B=%lf %lf %lf %le\n", B, C, A, Mp);
    concset = true;
    get_ri();
    get_rhoi();
    //cout<<conc<<endl;
    return conc;

 }

  float get_conc() {
    if (concset) return conc;
    else {
      cout << "set concentration using method 'cluster::concentration' first " << endl;
      return -1.0;
    }
  }

  float get_conc_slope(){
     if (concset) return B;
  }

float get_conc_z0(){
  //float con0= conc;
  //if (concset)

   return conc;
  }

  float get_mass_overden(float overden_to_calc) {

    double A=log(1+conc)- conc/(1+conc);
    double f= overden_to_calc/mass_overden*1.0/pow(conc, 3.0)*A;
    double a1= 0.5116;
    double a2= -0.4283;
    double a3= -3.13e-3;
    double a4= -3.52e-5;
    double p= a2+a3*log(f)+ a4*pow(log(f), 2.0);
    double x= pow(a1*pow(f,2.0*p) + pow(0.75, 4.0), -0.5) + 2*f;
    // cout<<conc<<endl;
    return mass*overden_to_calc/mass_overden*1.0/pow(conc*x,3);

  }

float get_rad_overden(float overden_to_calc){
   float massoverden= get_mass_overden(overden_to_calc);
   return pow(massoverden/(4./3.*PI*overden_to_calc*rho_crit) , 1.0/3.0);

}

};

#endif
