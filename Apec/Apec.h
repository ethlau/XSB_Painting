#ifndef _APEC_
#define _APEC_

#include "readapec.h"

#define NUMZ 30

struct EMISSION *apec_data;

void init_apec ();
void apec ( double *ebins, int nbins, double abun, double redshift, double temp, double *spec );

extern void readapec_getdata_extern(char *linefile,
                       char *cocofile);


extern void readapec_getdata(char *linefile,
                       char *cocofile,
                       struct EMISSION **apec_data);



extern void readapec_simple_spectrum(const struct EMISSION *apec_data,
                          int nbins,
                          double *ebins,
                          int hdu,
                          double *spectrum);

extern void readapec_calc_hdu_elem_spectrum(const struct EMISSION *apec_data,
                          int nbins,
                          double *ebins,
                          int hdu,
                          int Z,
                          double *spectrum);

extern void readapec_calc_hdu_ion_spectrum(const struct EMISSION *apec_data,
                          int nbins,
                          double *ebins,
                          int hdu,
                          int Z,
                          int rmJ,
                          double *spectrum);

extern void readapec_calc_hdu_ion_spectrum_part(const struct EMISSION *apec_data,
                          int nbins,
                          double *ebins,
                          int hdu,
                          int Z,
                          int rmJ,
                          int doline,
                          int docont,
                          double *spectrum);

extern int readapec_find_kT_hdu(const struct EMISSION *apec_data,
                double kT, int log);

extern void readapec_calc_ion_spectrum(const struct EMISSION *apec_data,
                       int nbins,
                       double *ebins,
                       double kT,
                       int Z,
                       int rmJ,
                       int doline,
                       int docont,
                       int nearest,
                       double *spectrum);

extern void readapec_calc_allion_spectrum(const struct EMISSION *apec_data,
                          int nbins,
                          double *ebins,
                          double kT,
                          int nearest,
                          int nelements,
                          int *Zlist,
                          int byion,
                          struct EMISSION_LIST **emission_list);

extern void readapec_calc_total_emission_abundance(struct EMISSION_LIST *emission_list,
                                   int nelements,
                                   int *Zlist,
                                   double *abundlist,
                                   double defaultabund,
                                   double *spectrum);


#endif
