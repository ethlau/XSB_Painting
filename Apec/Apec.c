#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <limits.h>
#include "Apec.h"


struct EMISSION *apec_data;

void init_apec () {

    /* Initialize apec calculation by reading in the required fits files */

    char linefile[MAXSTRLEN];
    char cocofile[MAXSTRLEN];

    //strncpy(linefile, "$ATOMDB/apec_line.fits",MAXSTRLEN);
    //strncpy(cocofile, "$ATOMDB/apec_coco.fits",MAXSTRLEN);
    //expand_env(linefile);
    //expand_env(cocofile);
    strncpy(linefile, LINE_FITS_FILE, MAXSTRLEN);
    strncpy(cocofile, COCO_FITS_FILE, MAXSTRLEN);

    readapec_getdata_extern(linefile, cocofile);

    return;

}

void close_apec() {
    if (apec_data != NULL) free(apec_data);
}

void apec ( double *ebins, int nbins, double abun, double temp, double redshift, double *spec ) {

    /* Main routine of computing APEC spectrum
     * Inputs:
     *  ebins: energy bins with size nbins+1
     *  nbins: size of the spectrum
     *  abun: metallcity of the as in solar units (assuming Anders & Grevesse 89) 
     *  temp: temperture of the gas in keV
     * Output:
     *  spec: spectrum in photons cm^3 / sec / bin
     */

    int i; 
    int nelements = NUMZ;
    int byion = 0;
    int nearest = 0;

    int Zlist[NUMZ];
    double abundlist[NUMZ];
    double defaultabund = 1.0;
    double spec_sum = 0.0;

    double ebins_emit[nbins+1];

    struct EMISSION_LIST *emission_list=NULL;

    for (i = 0; i < NUMZ; i++){
        Zlist[i] = i+1;
        if ( i > 1 ) {
            abundlist[i] = abun;
        } else {
            abundlist[i] = 1.0;
        }
    }

    // blueshift the energy bin from detector's frame to emitter's frame
    for (i = 0; i <= nbins; i++){
        ebins_emit[i] = ebins[i]*(1.0+redshift);
    }

    readapec_calc_allion_spectrum(apec_data, nbins, ebins_emit, temp, nearest, nelements, Zlist, byion, &emission_list);
    readapec_calc_total_emission_abundance (emission_list, nelements, Zlist, abundlist, defaultabund, spec);
 
    if (emission_list != NULL) {
        if ( emission_list->ion_emission != NULL){
            free(emission_list->ion_emission);
        }
        free(emission_list);
     }
     
    for ( i = 0; i < nbins; i++) {
        //redshifting spectrum to account for time dilation in photon arrival rate 1/dt 
        // spec[i] = spec[i]/(1.0+redshift);   
        spec_sum += spec[i];
    }
    
    //printf("spec_sum =  %e\n", spec_sum);

    return;
}

