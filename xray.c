#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "xray.h"
#include "Apec/Apec.h"

double tarray[ntmax]; //keV
double zarray[nzmax]; //Solar unit
double rarray[nrmax];
double lambda_table[ntmax][nrmax];
double tres, zres, eres, rres;
double keV2erg = 1.602177e-9;

void set_lambda_table ( double tarray[ntmax], double rarray[nrmax], double lambda_table[ntmax][nrmax] );
double int_lambda_table (double temp, double redshift, double tarray[ntmax], double rarray[nrmax], double lambda_table[ntmax][nrmax] );


void set_lambda_table ( double tarray[ntmax], double rarray[nrmax], double lambda_table[ntmax][nrmax] ){
    /* 
     * Tabulate X-ray emission integrated in receiver's energy range of 0.5-2.0 keV 
     * for a range of temperature and abundance   
     * Original spectrum in units of photon*cm**3/s (into 4pi angle) in emitter's frame
     * Converted to ergs*cm**3/s by multiplying E in keV and keV2erg unit conversion 
     */

    char filename[256];
    int neff = 0;
    double *energy = NULL;
    double *effarea = NULL;

    double temp, metal, redshift;
    double ebin[nemax+1],spec[nemax];
    double area[nemax];
    double total_area;

    double lambda;
    double velocity,dem,density;
    int ne, qtherm;

    int ie, je, it, iz, ir, i, j;

    init_apec();
    printf("Number of bins in arrays of Energy, Temperature, Redshift = %d, %d, %d\n", nemax, ntmax, nrmax);
    printf("Energy range = [%f,%f] keV\n", emin, emax);
    printf("Temperature range = [%f,%f] keV\n", tmin, tmax);
    printf("Redshift range = [%f,%f] \n", rmin, rmax);
    eres = (emax-emin)/nemax;
    printf("Energy resolution = %f [keV]\n", eres);
    for (ie = 0; ie < nemax+1; ie++) {
	ebin[ie] = emin + (double)ie*eres;
    }

    tres = (log10(tmax)-log10(tmin))/(ntmax);
    for (it = 0; it < ntmax; it++) {
	tarray[it] = pow(10.0,log10(tmin)+(double)it*tres);
	//printf("%e %e %e\n", tarray[it], tmin, pow(10.0,log10(tmin)+(double)it*tres));
    }

    zres = (zlmax-zlmin)/(nzmax);
    for (iz = 0; iz < nzmax; iz++) {
	zarray[iz] = pow(10.0,zlmin+(double)iz*zres);
	//printf("%e %e %e\n", zarray[iz], zlmin, pow(10.0,zlmin+(double)iz*zres));
    }

    rres = (log10(rmax)-log10(rmin))/(nrmax);
    for (ir = 0; ir < nrmax; ir++) {
	rarray[ir] = pow(10.0, log10(rmin)+(double)ir*rres);
	//printf("%e %e %e\n", zarray[iz], zlmin, pow(10.0,zlmin+(double)iz*zres));
    }
   

    velocity = 0.0;
    qtherm = 1;
    dem = 1.0;
    density = 1.0;
    ne = nemax;

    printf("Tabulating lambda... ");

    for (i = 0; i < ntmax; i++) {
	for (j = 0; j < nrmax; j++) {
            for (ie = 0; ie < nemax; ie++) {
	        spec[ie] = 0.0;
            }
 
	    temp = tarray[i];
	    //metal = zarray[j];
            metal = ABUNDANCE;
            redshift = rarray[j];
            //spec is in units of photons cm^3/s/bin  in receiver's frame
            //note that the photon arrival rate is already redshifted
            //ebin is already in receiver's frame
            //still need to redshift photon energy when converting from photon counts to flux
            apec ( ebin, ne, metal, temp, redshift, spec );
	    lambda = 0.0;
	    for (ie = 0; ie < nemax; ie++) {
                //printf("%e", spec[ie]);
                // photons -> ergs
		lambda += spec[ie]*(0.5*(ebin[ie]+ebin[ie+1]))*keV2erg/(1.+redshift);
	    }
            //printf("temp, z, lambda = %e %f %e\n", temp, redshift, lambda);
            
	    lambda_table[i][j] = lambda;
            assert(lambda >= 0);
	}
    }
    printf("done.\n");
    return;
}

double int_lambda_table (double temp, double redshift, double tarray[ntmax], double rarray[nrmax], double lambda_table[ntmax][nrmax] ) {

    int it,iz,ir;
    double lt1, lt2, lr1, lr2, lz1, lz2, f11, f12, f21, f22;
    double q1, q2;
    double emiss = 0.0;

    double ltemp, lredshift;

    ltemp = log10(temp);
    lredshift = log10(redshift);
    //lmetal = log10(metal);

    it = (int)((ltemp-log10(tmin))/tres);
    //iz = (int)((lmetal-zlmin)/zres);
    ir = (int)((lredshift-log10(rmin))/rres);


    //if (iz < 0 ) {
	//iz = 0;
    //}

    if (ir < 0 ) {
	ir = 0;
    }

    if (it < 0 ) {
        it = 0;
    }

    if (ir > nrmax -1) ir = nrmax - 1;
    if (it > ntmax -1) ir = ntmax - 1;


    if ( it >= 0 && it < ntmax-1  && ir >= 0 && ir < nrmax-1 ) { 

	//printf("%d %e %e %e\n",it,tarray[it],temp,tarray[it+1]);
	//printf("%d %e %e %e\n",iz,zarray[iz],metal,zarray[iz+1]);

	lt1 = log10(tarray[it]);
	lt2 = log10(tarray[it+1]);
	lr1 = log10(rarray[ir]);
	lr2 = log10(rarray[ir+1]);
	
	// bilinear interpolation
        //
        
        f11 = log10(lambda_table[it][ir]);
	f12 = log10(lambda_table[it][ir+1]);
	f21 = log10(lambda_table[it+1][ir]);
	f22 = log10(lambda_table[it+1][ir+1]);

	emiss = f11/((lt2-lt1)*(lr2-lr1))*((+lt2-ltemp)*(+lr2-lredshift))
	        + f21/((lt2-lt1)*(lr2-lr1))*((-lt1+ltemp)*(+lr2-lredshift))
		+ f12/((lt2-lt1)*(lr2-lr1))*((+lt2-ltemp)*(-lr1+lredshift))
        	+ f22/((lt2-lt1)*(lr2-lr1))*((-lt1+ltemp)*(-lr1+lredshift));

	emiss = pow(10.0,emiss);

    } else if (it == ntmax-1 && ir >= 0 && ir < nrmax-1) {

	lr1 = log10(rarray[ir]);
	lr2 = log10(rarray[ir+1]);

	f11 = log10(lambda_table[it][ir]);
	f12 = log10(lambda_table[it][ir+1]);
	
	emiss = f11 + (f12-f11)/(lr2-lr1)*(lredshift-lr1); 
	emiss = pow(10.0,emiss);

    } else if (it >= 0 && it < ntmax-1 && ir == nrmax-1) {

	lt1 = log10(tarray[it]);
	lt2 = log10(tarray[it+1]);

	f11 = log10(lambda_table[it][ir]);
	f21 = log10(lambda_table[it+1][ir]);
	
	emiss = f11 + (f21-f11)/(lt2-lt1)*(ltemp-lt1); 
	emiss = pow(10.0,emiss);

    } else if ( it == ntmax-1 && ir == nrmax-1 ) {

	emiss = (lambda_table[it][ir]);
	//emiss = pow(10.0,emiss);
    }

    //printf("%d %d %e %e %e %e\n", it, ir, tarray[it], rarray[ir], temp, redshift);
    //printf("%d %d %e %e %e\n", it, ir, lambda_table[it][ir], lambda_table[it+1][ir+1], emiss);
    //assert(emiss >= 0 && emiss < 1e-3);
    if (emiss != emiss) emiss = 1.0e-70;

    return emiss;
}


