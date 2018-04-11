#include <fitsio.h>
#define MAXSTRLEN 1024
#define READONLY  0    /* options when opening a file */
#define MAXARRFLOAT 10000
#define INTERNAL
struct EMISSION_LINES {
  double lambda;
  //double lambda_err;
  double epsilon;
  //double Plambda;
  int Z;
  int rmJ;
  int ul;
  int ll;
  //struct EMISSION_LINES *next;
};



struct COMPRESS_CONT { /* Holds one ion's continuum only */
  int     Z;
  int     rmJ;
  int     Ncont;
  double *E_cont;  /* in keV */
  double *cont;    /* in photons cm^3 / s / keV */
  int     Npseudo;
  double *E_pseudo;/* in keV */
  double *pseudo;  /* in photons cm^3 / s / keV */
  struct COMPRESS_CONT *next;
};
  
struct EMISSION_BLOCK {
  double kT; /* temperature of this block, keV */
  double Ne; /* density of this block, cm-3 */ 
  double time; /* time of this block, 0 if steady state */
  int nelement ; /* number of elements in continuum */
  int nline; /* number of lines */
  int ncont ; /* number of continuum data points */
  int npseudo ; /* number of pseduocontinuum data points */
  struct EMISSION_LINES *lines; /* pointer to the first line */
  struct COMPRESS_CONT *coco; /* pointer to the first continuum */
  struct EMISSION_BLOCK *next; /* pointer to the next emission block */
};
  

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
#define malloc_safe(d) (malloc_or_die(d , __FILE__))
char * malloc_or_die(size_t size, char *routine);

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

struct EMISSION_ION {
  int Z; /* element atomic number */
  int rmJ; /* ion number. 1=neutral, 2=+1 etc. 
            * 0 means all ions of the element together */
  double *emissivity; /* emissivity if ph cm3 s-1 bin-1 */
  struct EMISSION_ION *next;
};

extern const struct EMISSION ** const apec_data_global_pointer;
#define apec_data_memstore (*apec_data_global_pointer)
