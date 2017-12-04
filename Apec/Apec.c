#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <limits.h>
#include "Apec.h"

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

    readapec_getdata_extern(linefile,cocofile);

    //readapec_getdata(linefile,cocofile,&apec_data);

    return;

}

void apec ( double *ebins, int nbins, double abun, double redshift, double temp, double *spec ) {

    /* Main routine of computing APEC spectrum
     * Inputs:
     *  ebins: energy bins with size nbins+1
     *  nbins: size of the spectrum
     *  abun: metallcity of the as in solar units (assuming Anders & Grevesse 89) 
     *  redshift: redshift of the gas
     *  temp: temperture of the gas in keV
     * Output:
     *  spec: spectrum
     */

    int i; 
    int nelements = NUMZ;
    int byion = 0;
    int nearest = 0;

    int Zlist[NUMZ];
    double abundlist[NUMZ];
    double defaultabund = 1.0;

    struct EMISSION_LIST *emission_list=NULL;

    for (i = 0; i < NUMZ; i++){
        Zlist[i] = i+1;
        if ( i > 1 ) {
            abundlist[i] = abun;
        } else {
            abundlist[i] = 1.0;
        }
    }

    readapec_calc_allion_spectrum(apec_data, nbins, ebins, temp, nearest, nelements, Zlist, byion, redshift, &emission_list);
    readapec_calc_total_emission_abundance (emission_list, nelements, Zlist, abundlist, redshift, defaultabund, spec);

    if (emission_list != NULL) {
        if ( emission_list->ion_emission != NULL){
            free(emission_list->ion_emission);
        }
        free(emission_list);
     }

    //for ( i = 0; i < nbins; i++) printf("%d: %f %e\n",i, ebins[i], spec[i]);

    return;
}

