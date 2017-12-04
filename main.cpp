#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include "cosmo.h"
#include "cluster.h"
#include "gas_model.h"
#include "read_lightcone.h"
#include "xray.h"

#define MAXBINS 500

#define M500_THRESHOLD 1e13

double tarray[ntmax]; //keV
double zarray[nzmax]; //Solar unit
double rarray[nrmax]; 
double lambda_table[ntmax][nrmax];
double tres, zres, eres;
 
const double PI = 4.0*atan(1.0);
const double megapc = 3.0857e24; // in cm
float periodic(float x, float L);
void emission_projection (double* rbins, double* r_in, double* r_out, double* profile, double* proj_prof, int nbins);
void temperature_projection (double* rbins, double* r_in, double* r_out, double* profile, double* proj_prof, double* weight, int nbins);

using namespace std;
int main(int argc, char *argv[]){
    if ( argc!=3 ){
        fprintf(stderr,"usage: %s <directory containing halo catalog> <halo catalog id>\n", argv[0]);
        exit(1);
    } 

    //float H0 = 67.27, Omega_M = 0.3156, Omega_b = 0.04917;
    float H0 = 70.0, Omega_M = 0.2791, Omega_b = 0.04917;
    float wt = -1.0, Omega_k = 0.00;
    float conc_norm = 1.0, conc_mass_norm = 1.0;
    float delta_rel = 0.18, delta_rel_n = 0.8, delta_rel_zslope = 0.5;
    float ad_index = 5.0;
    float eps_fb = 3.97e-6;
    float eps_dm = 0.0;
    float fs_0 = 0.026;
    float fs_alpha = 0.12;
    int pturbrad = 2;
    bool verbose = false;
    float Mvir, Rvir, M500, R500, Rscale, conc, z, a, cosmic_t, cosmic_t0;
    float h = H0/100.0, E;
    // set cluster overdensity
    // this is the overdensity within which mass defined (i.e. \Delta)
    // set to -1.0 for virial radius, or 200 for M200 (rhocrit)
    float overden_id = -1.0; // 200 for delta=200 rho-c , -1 for delta=vir x rho-c
    int relation = 3; // concentration relation
    float rcutoff = 2.0;

    char halo_run[256];
    char root[1024], filename[1024], outprofname[1024], outhaloname[1024];
    float xpos, ypos, zpos, L;
    double x, Rmax, Pgas, rho, rhogas, Yanl, Ypart, m_p;
    gsl_integration_glfixed_table *t;
    double Time, Redshift;
    float pos[3];
    double UnitLength_in_cm= 3.085678e24;   /*  code length unit in cm/h */
    double UnitMass_in_g= 1.989e43;         /*  code mass unit in g/h */
    double UnitDensity_in_cgs = UnitMass_in_g/ pow(UnitLength_in_cm,3);

    int    nbins = 100;
    double delx;
    double r_in[MAXBINS], r_out[MAXBINS], dvol[MAXBINS];
    double rbins[MAXBINS], emiss_prof[MAXBINS], sb_prof[MAXBINS], emiss_measure[MAXBINS];
    double kT[MAXBINS], ngas[MAXBINS];
    double kT_2D[MAXBINS], ngas_2D[MAXBINS];
    double angbins[MAXBINS], ang_in[MAXBINS], ang_out[MAXBINS], solid_angle[MAXBINS];
    int i, j;
    FILE *outprof, *outhalo;

    sprintf(root, "%s", argv[1]);
    sprintf(halo_run, "%s", argv[2]);
    sprintf(outhaloname, "sbhalo_run%s.txt", halo_run);
    sprintf(outprofname, "sbprof_run%s.txt", halo_run);

    sprintf(filename, "%s/run%s.h_halo", root, halo_run);
    halo_list *halos;
    halo_struct *halo;

    halos = load_halo_run(filename);
    if(halos == NULL ){
        printf("!memory allocation failure\n");
        exit(1);
    }
    /* Set up X-ray emissivity table */
    printf("Setting up X-ray emission table\n");
    set_lambda_table(tarray,rarray,lambda_table); 

    outprof = fopen (outprofname, "w");
    outhalo = fopen (outhaloname, "w");

    fprintf(outhalo,"# halo_id lens_plane_id theta_x theta_y redshift M500 [Msun] R500 Rvir Rscale [Mpc]\n");
    fprintf(outprof,"# r_in r_mid r_out [Mpc] ang_in ang_mid ang_out [arcsecs] Sx [cts/s/cm^2/arcsec^2] kT_projected [keV]\n");
    for( i=0; i < halos->num_halos; i++){
    //for( i=0; i < 100; i++){
        halo = &(halos->list[i]);
        if(halo->M500c/h < M500_THRESHOLD) continue;
        
        for (j = 0; j < MAXBINS; j++) {
            rbins[j] = 0.0;
            r_in[j] = 0.0;
            r_out[j] = 0.0;
            dvol[j] = 0.0;
            sb_prof[j] = 0.0;
            emiss_prof[j] = 0.0;
            kT[j] = 0.0;
            ngas[j] = 0.0;
            kT_2D[j] = 0.0;
            ngas_2D[j] = 0.0;
            ang_in[j] = 0.0;
            ang_out[j] = 0.0;
            angbins[j] = 0.0;
        }

        printf("Computing SB profile for halo %d\n", i);
        Redshift = halo->redshift;

        cosmo cosm_model(H0, Omega_M, Omega_b, Omega_k, wt);
        cosmic_t = cosm_model.cosmic_time(Redshift);
        cosmic_t0 = cosm_model.cosmic_time(0.0);
        E = cosm_model.Efact(Redshift);

        Mvir = halo->Mvir/h;
        cluster nfwclus(Mvir, Redshift, overden_id, relation, cosm_model);

        // Here, use halo concentration from the halo catalog instead of the M-c relation from Duffy+08
        conc = halo->rvir/halo->rs;
        nfwclus.set_conc(conc);
        // set halo concentration using M-c relation of Duffy et al (2008)
        //conc = nfwclus.concentration(conc_norm, conc_mass_norm);
        //M500 = nfwclus.get_mass_overden(500.0); // M500 in Msol (for calculating stellar mass frac)

        M500 = halo->M500c/h;
        R500 = pow(M500/((4.0/3.0)*M_PI*500.0*cosm_model.calc_rho_crit(z)), 1.0/3.0);
        Rvir = nfwclus.get_radius();
        Rscale = halo->rs/1000;

        printf("Halo parameters: \n");
        printf("M500 = %e, Redshift = %f\n", halo->M500c/h, Redshift);
        printf("Rscale = %f, Rvir = %f, R500 = %f, Cvir = %f\n", Rscale, Rvir, R500,  conc);

        gas_model icm_mod(delta_rel, ad_index, eps_fb, eps_dm, fs_0, fs_alpha, pturbrad, delta_rel_zslope, delta_rel_n);

        icm_mod.calc_fs(M500, Omega_b/Omega_M, cosmic_t0, cosmic_t);
        icm_mod.evolve_pturb_norm(z, rcutoff);
        icm_mod.set_nfw_params(Mvir, Rvir, conc, nfwclus.get_rhoi(), R500);
        icm_mod.set_mgas_init(Omega_b/Omega_M);
        icm_mod.findxs();

        icm_mod.solve_gas_model(verbose, 1e-5);
        Rmax = icm_mod.thermal_pressure_outer_rad()*Rvir;
        Yanl = icm_mod.calc_Y(R500, Rvir, Rmax);
        for (j = 0; j < nbins; j++) {
            
            delx = (log10(5.0*Rvir)-log10(0.01*Rvir))/nbins;
            r_out[j] = pow(10.0, log10(0.01*Rvir) + (float)j*delx); // outer edge of the radial bin in Mpc
            rbins[j] = pow(10.0, log10(0.01*Rvir) + ((float)j-0.5)*delx); // midpt of the radial bin in Mpc
            if (j == 0) { 
                r_in[j] = 0.0; // inner edge of the radial bin in Mpc
                dvol[j] = (4.0/3.0)*M_PI*pow(r_out[j], 3.0); // volume of the radial bin in Mpc^3
            } else {
                r_in[j] = r_out[j-1];
                dvol[j] = (4.0/3.0)*M_PI*(pow(r_out[j], 3.0) - pow(r_out[j-1], 3.0));

            }

            //4 pi steradians in arcsec^2
            solid_angle[j] = 4.0*M_PI*180.0*3600.0/M_PI *180.0*3600.0/M_PI;

            emiss_prof[j] = icm_mod.calc_xray_emissivity(rbins[j], R500, Redshift);
            kT[j] = icm_mod.calc_gas_temperature (rbins[j], R500);
            ngas[j] = icm_mod.calc_gas_num_density (rbins[j], R500);
            angbins[j] = rbins[j] / cosm_model.ang_diam(Redshift) *180.0*3600.0/M_PI;// in arcsecs
            ang_in[j] = r_in[j] / cosm_model.ang_diam(Redshift) *180.0*3600.0/M_PI;// in arcsecs
            ang_out[j] = r_out[j] / cosm_model.ang_diam(Redshift) *180.0*3600.0/M_PI;// in arcsecs
            
        }
        // project to emissivity get surface brightness profile (sec 5.5.4 from Sarazin 88)
        emission_projection(rbins, r_in, r_out,  emiss_prof, sb_prof, nbins); 
        // 2D projected emission-weighted temperature 
        temperature_projection(rbins, r_in, r_out, kT, kT_2D, emiss_prof, nbins); 
        
        fprintf(outhalo,"%d %d %d %d %f %e %f %f %f\n", i, halo->lens_id, halo->theta_x, halo->theta_y, Redshift, M500, R500, Rvir, Rscale);
        fprintf(outprof,"# %d\n", i);
        for (j = 0; j < nbins; j++) {
            fprintf(outprof,"%f %f %f %f %f %f %e %e\n", r_in[j], rbins[j], r_out[j], ang_in[j], angbins[j], ang_out[j], sb_prof[j]/solid_angle[j], kT_2D[j]);
        }
    }
    fclose(outhalo);
    fclose(outprof);
    destroy_halo_list(halos);
    return 0;
}

void emission_projection (double* rbins, double* r_in, double* r_out, double* profile, double* proj_prof, int nbins){
    int i,j;
    double proj_R[MAXBINS];
    double proj_sum;
    double delr;


    for (i=0; i<MAXBINS; i++) {
        proj_prof[i] = 0.0;
        proj_R[i] = rbins[i];
    }

    // Eq 5.8 from Sarazin 88
    // I = 2 * int phi *r *dr/sqrt(r^2-R^2)
    for (j=0; j<nbins; j++) {
        proj_sum = 0.0;
        for (i=j; i<nbins; i++){
            delr = r_out[i] - r_in[i];
            proj_sum += 2.0 * rbins[i] * delr * profile[i]/ sqrt(r_out[i]*r_out[i]-proj_R[j]*proj_R[j]);
        }
        proj_prof[j] = proj_sum * megapc;
    }

}

void temperature_projection (double* rbins, double* r_in, double* r_out, double* profile, double* proj_prof, double* weight, int nbins){
    int i,j;
    double proj_R[MAXBINS];
    double proj_sum;
    double weight_sum;
    double delr;


    for (i=0; i<MAXBINS; i++) {
        proj_prof[i] = 0.0;
        proj_R[i] = rbins[i];
    }

    // Eq 5.8 from Sarazin 88
    // I = 2 * int phi *r *dr/sqrt(r^2-R^2)
    for (j=0; j<nbins; j++) {
        proj_sum = 0.0;
        weight_sum = 0.0;
        for (i=j; i<nbins; i++){
            delr = r_out[i] - r_in[i];
            proj_sum += 2.0 * rbins[i] * delr * profile[i] * weight[i]/ sqrt(r_out[i]*r_out[i]-proj_R[j]*proj_R[j]);
            weight_sum += 2.0 * rbins[i] * delr * weight[i]/ sqrt(r_out[i]*r_out[i]-proj_R[j]*proj_R[j]);

        }
        if ( weight_sum > 0 ) {
            proj_sum /= weight_sum;
        } else {
            proj_sum = 0.0;
        }
        if ( proj_sum > 0) proj_prof[j] = proj_sum ; 
    }

}

float periodic(float x, float L){
    float y = x;

    if(y>0.5*L) y-=L;

    if(y<-0.5*L) y+=L;

    return y;
}

