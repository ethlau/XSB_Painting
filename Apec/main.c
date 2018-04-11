#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <limits.h>
#include "readapec.h"

#define emax             (2.000000) //keV
#define emin             (0.500000) //keV
#define eres             (1.00e-3) //keV
#define nsize            3

typedef struct {
       float *ebin;
       float *(pd[nsize][nsize]);
} struct_spectrum;


void init_apec();
void close_apec();
void apec ( double *ebins, int nbins, double abun, double temp, double redshift, double spec[] );

int set_energy_bins(double **ebin){

    int i;
    int ne, iemin, iemax;
    double e, logstep, linstep;
    double eold,z;

    ne = 0;
    e = emin;
    
    iemin = ne;

    while ( e < emax ) {
	eold = e;
	(*ebin) = (double *)realloc((*ebin),(ne+1)*sizeof(double));
        if ((*ebin) == NULL) {
	    printf("3:ebin is NULL!\n");
	    exit(1);
        }
        (*ebin)[ne] = eold;
        e = e + eres;
        ne++;
    }
    ne--;
    iemax = ne;
    printf("ne = %d\n", ne);
    printf("iemin, emin = %d %e\n", iemin, (*ebin)[iemin]);
    printf("iemax, emax = %d %e\n", iemax, (*ebin)[iemax]);
    return ne;

}

int main() {

    int ne;
    int i, j, k,ie, je;
    double *ebin = NULL;
    double delta_e;

    double abun, inputt;

    double p;
    double sum;
    double z_emit;

    double *spec;
        
    int num_grid = nsize;
    int num_mesh;

    double *mesh;

    struct_spectrum *spectrum;	

    FILE *output;

    char filename[1024];

    /* Set up energy bins. ne is the number of energy bins. Note ebins has size of ne+1.  */
    ne = set_energy_bins(&ebin);

    num_mesh = num_grid*num_grid*ne;

    /* input temperature of 10 keV, metallicity of 0.2 Solar */
    inputt = 10.0;
    abun = 0.20;
    z_emit = 0.06; // set emission redshift at 0.06


    for (i = 0; i < nsize; i++) {
        init_apec();
        printf("%d\n", i);    
        spec = (double *)malloc(ne*sizeof(double));
        
        for (ie = 0; ie < ne; ie++) {
            spec[ie] = 0.0;
        }

        /* spectrum in restframe */
        apec (ebin, ne, abun, inputt, z_emit, spec);
        
        // spectrum in counts / sec / cm^3 / keV
        for (ie=0; ie < ne; ie++) {
            delta_e = (ebin[ie+1]-ebin[ie]); 
	    spec[ie] = spec[ie]/delta_e;
        }

        // printing out the spectrum into an ascii file 
        sprintf(filename, "spec_%d.dat", i);
        output = fopen(filename, "w");
        for (ie=0;ie<ne;ie++) {
            fprintf(output,"%e %e\n",0.5*(ebin[ie]+ebin[ie+1]),spec[ie]);
        }

        fclose (output);

        free(spec);
 
        
    }

    free(ebin);
    return 0;

}

