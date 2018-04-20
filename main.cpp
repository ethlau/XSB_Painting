#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "cosmo.h"
#include "cluster.h"
#include "gas_model.h"
#include "io/read_halo.h"
#include "xray.h"
#include "ConfigParser/ConfigParser.h"

#define MAXBINS 500

double tarray[ntmax]; //keV
double zarray[nzmax]; //Solar unit
double rarray[nrmax]; 
double lambda_table[ntmax][nrmax];
double tres, zres, eres;
 
const double PI = 4.0*atan(1.0);
const double megapc = 3.0857e24; // in cm
float periodic(float x, float L);
void luminosity_to_SB (double redshift, double* rbins, double* r_in, double* r_out, double* luminosity, double* surface_brightness, int nbins);
void emission_projection (double* rbins, double* r_in, double* r_out, double* profile, double* proj_prof, int nbins);
void temperature_projection (double* rbins, double* r_in, double* r_out, double* profile, double* proj_prof, double* weight, int nbins);
float rscale_from_mass(float m500, float z, float rhocrit, float h);
double vikhlinin_lum ( double M500, double redshift);

using namespace std;
int main(int argc, char *argv[]){

    if ( argc !=4 ){
        fprintf(stderr,"usage: %s <config file> <lightcone folder> <lightcone id>\n", argv[0]);
        exit(1);
    }
   
    printf("read_config ...");    
    read_config(argv[1]);
    printf("done\n");
    char *file_format = config_get_string("io","file_format");
    char root[512];
    char halo_run[512];
    //char *root = config_get_string("io","directory");
    //char *halo_run = config_get_string("io","identifier");
    float mass_threshold = config_get_float("halo","mass_threshold");
    float H0 = config_get_float("cosmo","H0");
    float Omega_M = config_get_float("cosmo","Om");
    float Omega_b = config_get_float("cosmo","Ob");


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

    char filename[512], outprofname[512], outhaloname[512];
    float xpos, ypos, zpos, L;
    double x, Rmax, Pgas, rho, rhogas, Yanl, Ypart, m_p;
    gsl_integration_glfixed_table *t;
    double Time, Redshift;
    float pos[3];
    double UnitLength_in_cm= 3.085678e24;   /*  code length unit in cm/h */
    double UnitMass_in_g= 1.989e43;         /*  code mass unit in g/h */
    double UnitDensity_in_cgs = UnitMass_in_g/ pow(UnitLength_in_cm,3);

    int nbins = 500;
    double delx;
    double r_in[MAXBINS], r_out[MAXBINS], dvol[MAXBINS];
    double ang_in[MAXBINS], ang_out[MAXBINS], angbins[MAXBINS], solid_angle[MAXBINS];
    double rbins[MAXBINS], emiss_prof[MAXBINS], sb_prof[MAXBINS], emiss_measure[MAXBINS];
    double kT[MAXBINS], ngas[MAXBINS], Lx[MAXBINS];
    double kT_2D[MAXBINS], ngas_2D[MAXBINS];
    double total_Lx, Lx_vik;
    int i, j;
    FILE *outprof, *outhalo;

    printf("file format is %s\n", file_format);
    printf("Halo mass threshod = %e [Msun]\n", mass_threshold);
    if ( strcmp(file_format,"simple") !=0 ){
        
        sprintf(root, "%s", argv[2]);
        sprintf(halo_run, "%s", argv[3]);

        sprintf(outhaloname, "sbhalo_run_%s.txt", halo_run);
    	sprintf(outprofname, "sbprof_run_%s.txt", halo_run);
 
        outprof = fopen (outprofname, "w");
        outhalo = fopen (outhaloname, "w");
		
	    fprintf(outprof,"# r_in r_mid r_out [Mpc] ang_in ang_mid ang_out [arcsecs] Lx [ergs/s] kT [keV] n_gas [cm^-3] SB [erg/s/cm^2/arcsec^2]\n");
    }


    halo_list *halos;
    halo_struct *halo;

    if ( strcmp(file_format,"rockstar" ) == 0) {
	    halos = load_halo_rs(root);
    } else if ( strcmp(file_format,"lightcone") == 0 ) {

	    sprintf(filename, "%s/run%s.h_halo", root, halo_run);
        printf("Reading lightcone file %s\n", filename);
	    halos = load_halo_run(filename);
    } else if ( strcmp(file_format,"simple") == 0 ) {
	    halos = load_halo_simple(root);
    } else {
	    fprintf(stderr,"Not a supported file format for halo catalog.\n");
        exit(1);
    }

    if( halos == NULL ){
        fprintf(stderr,"memory allocation failure!\n");
        exit(1);
    }

    if ( &mass_threshold == NULL) {
        mass_threshold = 1.e13;
    }
    printf("using mass threshold of M500 > %e Msun\n", mass_threshold);

    cout << "Done halo reading" << endl;

    /* Set up X-ray emissivity table */
    if (strcmp(file_format,"simple")!=0){
	    cout << "Setting up X-ray emission table" << endl;
	    set_lambda_table(tarray,rarray,lambda_table); 
    }

    if (strcmp(file_format,"rockstar")==0) {
	    fprintf(outhalo,"#ID PID X Y Z[Mpc] redshift Rvir Rs R500[kpc] Mvir M200c M500c[Msun] Xoff\n");
    }
    else if (strcmp(file_format,"lightcone")==0) {
	    fprintf(outhalo,"# halo_id lens_plane_id theta_x theta_y redshift M500 [Msun] R500 Rvir Rscale [Mpc] Lx Lx_Vik [erg/s]\n");
    }

    for( i=0; i < halos->num_halos; i++){
        // cout << "Computing gas profile for halo " << i << endl;
        halo = &(halos->list[i]);
        if( halo->M500c/h < mass_threshold ) { 
            //cout << "Mass is " << halo->M500c/h << endl;
	        //cout << "Skipped" <<endl;
	        continue;
        }

        for (j = 0; j < MAXBINS; j++) {
            rbins[j] = 0.0;
            r_in[j] = 0.0;
            r_out[j] = 0.0;
            dvol[j] = 0.0;
            sb_prof[j] = 0.0;
            emiss_prof[j] = 0.0;
            Lx[j] = 0.0;
            kT[j] = 0.0;
            ngas[j] = 0.0;
            kT_2D[j] = 0.0;
            ngas_2D[j] = 0.0;
            ang_in[j] = 0.0;
            ang_out[j] = 0.0;
            angbins[j] = 0.0;
        }

        Redshift = halo->redshift;

        cosmo cosm_model(H0, Omega_M, Omega_b, Omega_k, wt);
        cosmic_t = cosm_model.cosmic_time(Redshift);
        cosmic_t0 = cosm_model.cosmic_time(0.0);
        E = cosm_model.Efact(Redshift);
	    cluster nfwclus(0.0, Redshift, overden_id, relation, cosm_model);

	    if (strcmp(file_format,"simple")!=0){
	        if (halo->Mvir < 1) 
                Mvir = (4.0/3.0)*M_PI*cosm_model.Delta_vir(Redshift)*cosm_model.calc_rho_crit(Redshift)*pow(halo->rvir/(1000.0*h),3.0);
	        else 
                Mvir = halo->Mvir/h;

	        nfwclus.reset_cluster(Mvir, Redshift, overden_id, relation, cosm_model);
	        conc = halo->rvir/halo->rs;
	        nfwclus.set_conc(conc);
			
	        M500 = halo->M500c/h;
	        R500 = pow(M500/((4.0/3.0)*M_PI*500.0*cosm_model.calc_rho_crit(Redshift)), 1.0/3.0);
	        Rvir = halo->rvir/h/1000.0;
	        Rscale = halo->rs/1000.0;

	    } else {
	        M500 = halo->M500c/h;
	        Mvir = M500;
	        Rscale = halo->rs/1000;
	        //Rscale = rscale_from_mass(M500*h, Redshift, cosm_model.calc_rho_crit(Redshift), h);
	        R500 = pow(M500/((4.0/3.0)*M_PI*500.0*cosm_model.calc_rho_crit(Redshift)), 1.0/3.0);
	        conc = R500/Rscale;
	        nfwclus.reset_cluster(Mvir, Redshift, 500.0, relation, cosm_model);
	        nfwclus.set_conc(conc);
	        Rvir = halo->rvir/h/1000.0;
	        //Rvir = nfwclus.get_radius();
	    }
        /* Here, use halo concentration from the halo catalog instead of the M-c relation from Duffy+08
           set halo concentration using M-c relation of Duffy et al (2008) */
        //conc = nfwclus.concentration(conc_norm, conc_mass_norm);
        //M500 = nfwclus.get_mass_overden(500.0); // M500 in Msol (for calculating stellar mass frac)

        cout << "Halo parameters: " << endl;
        printf("M500 = %e, Redshift = %f, ", halo->M500c/h, Redshift);
	    printf("Rscale = %f, Rvir = %f, R500 = %f, Cvir = %f\n", Rscale, Rvir, R500,  conc);

        gas_model icm_mod(delta_rel, ad_index, eps_fb, eps_dm, fs_0, fs_alpha, pturbrad, delta_rel_zslope, delta_rel_n);

        icm_mod.calc_fs(M500, Omega_b/Omega_M, cosmic_t0, cosmic_t);
        icm_mod.evolve_pturb_norm(Redshift, rcutoff);
        icm_mod.set_nfw_params(Mvir, Rvir, conc, nfwclus.get_rhoi(), R500);
        icm_mod.set_mgas_init(Omega_b/Omega_M);
        icm_mod.findxs();

        icm_mod.solve_gas_model(verbose, 1e-5);
        Rmax = icm_mod.thermal_pressure_outer_rad()*Rvir;
        Yanl = icm_mod.calc_Y(R500, Rvir, Rmax);

        total_Lx = 0.0;
        
        Lx_vik = vikhlinin_lum ( M500, Redshift ); 
        
        for (j = 0; j < nbins; j++) {
            
            delx = (log10(5.0*R500)-log10(0.01*R500))/nbins;
            r_out[j] = pow(10.0, log10(0.01*Rvir) + (float)j*delx); // outer edge of the radial bin in Mpc
            rbins[j] = pow(10.0, log10(0.01*Rvir) + ((float)j-0.5)*delx); // midpt of the radial bin in Mpc

            if (j == 0) { 
                r_in[j] = 0.0; // inner edge of the radial bin in Mpc
                dvol[j] = (4.0/3.0)*M_PI*pow(r_out[j], 3.0); // volume of the radial bin in Mpc^3
            } else {
                r_in[j] = r_out[j-1];
                dvol[j] = (4.0/3.0)*M_PI*(pow(r_out[j], 3.0) - pow(r_out[j-1], 3.0));
            }

            kT[j] = icm_mod.calc_gas_temperature (rbins[j], R500);

	        if (strcmp(file_format,"simple")!=0){
	            ngas[j] = icm_mod.calc_gas_num_density (rbins[j], R500);
                solid_angle[j] = 4.0*M_PI*180.0*3600.0/M_PI *180.0*3600.0/M_PI;
                emiss_prof[j] = icm_mod.calc_xray_emissivity(rbins[j], R500, Redshift);
                Lx[j] = emiss_prof[j] * dvol[j] * pow(megapc,3.0);
                total_Lx += Lx[j];
                angbins[j] = rbins[j] / cosm_model.ang_diam(Redshift) *180.0*3600.0/M_PI;// in arcsecs
                ang_in[j] = r_in[j] / cosm_model.ang_diam(Redshift) *180.0*3600.0/M_PI;// in arcsecs
                ang_out[j] = r_out[j] / cosm_model.ang_diam(Redshift) *180.0*3600.0/M_PI;// in arcsecs
	        } else {
		        ngas[j] = icm_mod.calc_gas_density (rbins[j], R500);
	        }    
        }

	    if (strcmp(file_format,"simple")!=0){
            // project to emissivity get surface brightness profile (sec 5.5.4 from Sarazin 88)
            luminosity_to_SB(Redshift, rbins, r_in, r_out, Lx, sb_prof, nbins); 
	    }
    
	    if (strcmp(file_format,"rockstar")==0) {
            fprintf(outhalo,"%ld %d %f %f %f %f %f %f %f %e %e %e %f \n", 
                    halo->id,halo->pid,halo->x,halo->y,halo->z,halo->redshift,
                    Rvir,Rscale,R500,Mvir,halo->M200c/h,halo->M500c/h,halo->Xoff);
	    } else if (strcmp(file_format,"lightcone")==0) {
	        fprintf(outhalo,"%d %d %d %d %f %e %f %f %f %e %e\n", 
            i, halo->lens_id, halo->theta_x, halo->theta_y, Redshift, M500, R500, Rvir, Rscale, total_Lx, Lx_vik);
	    }

	    if (strcmp(file_format,"simple")==0){
	        cout << "writing output for simple " << endl;
	        sprintf(outprofname, "prof_run_%ld.txt", halo->id);
	        outprof = fopen (outprofname, "w");
	        fprintf(outprof,"# %ld\n", halo->id);
	        fprintf(outprof,"# M=%e[Msun] z=%f c=%f\n", M500, Redshift, conc);
	        fprintf(outprof,"# r_in r_mid r_out [Mpc] kT[keV] dgas[g/cm3] \n");
	    } else {
            fprintf(outprof,"# %ld\n", halo->id);
	    }
		
        for (j = 0; j < nbins; j++) {
	        if (strcmp(file_format,"simple")==0){
                fprintf(outprof,"%f %f %f %e %e\n", r_in[j], rbins[j], r_out[j], kT[j], ngas[j]);
	        } else {
		        fprintf(outprof,"%f %f %f %f %f %f %e %e %e %e\n", 
                  r_in[j], rbins[j], r_out[j], ang_in[j], angbins[j], ang_out[j], Lx[j], kT[j], ngas[j], sbprof[j]);
	        }
        }
    }
	if (strcmp(file_format,"simple")==0) {
        fclose(outprof);
    } else if (strcmp(file_format,"simple")!=0){
	    fclose(outhalo);
	    fclose(outprof);		
    }

    destroy_halo_list(halos);

    return 0;
}

void luminosity_to_SB (double redshift, double* rbins, double* r_in, double* r_out, double* luminosity, double* surface_brightness, int nbins){

    /* luminosity in erg/s in rest-frame of the source in the observer's waveband
     * output is surface brightness in the observer's frame; using Equation 3.53 from Mo, White, van den Bosch
     * unit of surface brightnes is erg /s /cm^2 
     */
    int i,j;
    double delr;


    for (i=0; i<MAXBINS; i++) {
        surface_brightness[i] = 0.0;
    }

    for (i=0; i<nbins; i++) {
        delr =(r_out[i] - r_in[i]) * megapc;
        surface_brightness[j] = luminosity[i]/(M_PI*M_PI*delr*delr)/pow((1.0+redshift),4); // erg/s/cm^2/radian^2
        surface_brightness[j] /= (180.0*3600.0/M_PI)*(180.0*3600.0/M_PI); //converting to erg/s/cm^2/arcsec^2
    
    }

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

float rscale_from_mass(float m500, float z, float rhocrit, float h){
	//use m200-c200 relation from dutton-maccio, and convert them to m500 and rs according to nfw profile
	float a, b, a1, a2, a3, a4, f, p;
	float deltadex = 0.02; //define the bin in interpolation
	float lgm200 = 10.0;
	float c200 = 0.0;
	float m500_calc=0.0;
	float lgm500_calc = 0.0;
	float lgm500_calcp, lgm200p, c200p; 
	float lgc200_c, lgm200_c;
	a = 0.52+(0.905-0.52)*exp(-0.617*pow(z,1.21));
	b = -0.101+0.026*z;
	a1 = 0.5116; 
	a2 = -0.4283;
	a3 = -3.13E-3;
	a4 = -3.52E-5;
	while (m500_calc<m500){
		c200p = c200;
		lgm200p = lgm200;
		lgm500_calcp = lgm500_calc;
		lgm200 += deltadex;
		
		c200 = pow(10.0,(b*(lgm200-12.0)+a));
		f = 2.5*(log(1+c200)-c200/(1+c200))/pow(c200,3.0);
		p = a2 + a3*log(f) + a4*log(f)*log(f);
		lgm500_calc = lgm200 + log10(2.5) - 3*log10(c200) -3*log10(pow(a1*pow(f,2*p)+0.75*0.75,-0.5)+2*f);
		m500_calc = pow(10.0,lgm500_calc);
	}
	
	lgc200_c = log10(c200p) + log10(c200/c200p)*(log10(m500)-lgm500_calcp)/(lgm500_calc-lgm500_calcp);
	lgm200_c = lgm200p + (lgm200-lgm200p)*(log10(m500)-lgm500_calcp)/(lgm500_calc-lgm500_calcp);
	
	return pow(3.0*(pow(10.0,lgm200_c)/h)/(200*4.0*M_PI*rhocrit),1.0/3)/pow(10.0,lgc200_c);
}

double vikhlinin_lum ( double M500, double redshift) {
    double Lx, Ez, aexp;

    Ez = sqrt(0.27 * pow((1.0 + redshift), 3.0) + 0.73);
    Lx = 47.392 + 1.61*log(M500) + 1.850 * log(Ez) - 0.39*log(0.70/0.72);

    Lx = exp(Lx);
    
    return Lx;

}
