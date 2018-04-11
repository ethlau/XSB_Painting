/* This is the header file for the readapec library
 * All publicly accessible routines will be in readapec.c
 * Everything else is (in theory) private, obviously they
 * can be linked to but there is no guarantee I won't mess with
 * their interfaces */
#include <fitsio.h>
#include <math.h>
#ifndef MAXSTRLEN
#define MAXSTRLEN 1024
#endif

#ifndef INTERNAL
/* This avoids redefining the same code in the external and internal linking
  There is clearly a better way to do this...*/
  
  struct EMISSION {
    int nblock; /* number of blocks */ 
    int *nline; /* number of lines */
    int *nelement ; /* number of lines in continuum */
    int *ncont ; /* number of lines in continuum */
    int *npseudo ; /* number of lines in continuum */
    double *kT ; /* number of lines in continuum */
    double *eDensity ; /* number of lines in continuum */
    double *time ; /* number of lines in continuum */
    struct EMISSION_BLOCK *emdat; /* pointer to the first block */
  };
  struct EMISSION_LIST {
    int nelement ; /* number of lines in continuum */
    int *Zlist; /* number of blocks */ 
    double kT ; /* temperature of data (keV) */
    double eDensity ; /*  density of data (cm-3)*/
    double time ; /* time of data (s) */
    int nbins ; /* number of bins for the data */
    double *binedges; /* the edges of each bins. Should be nbins+1*/
    struct EMISSION_ION *ion_emission; /* Block with the emission from an ion */
  };
#endif

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


