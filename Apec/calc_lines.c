#include "readapec_spectrum.h"

#define keV_erg 1.60217657e-9
#define c_cm_s  2.99792458e10
#define amu_g   1.66053892e-24

//keV to Angstrom
#define keVtoA  12.398425

void calcGaussianLine(double *energyArray, int nbins, double ecenter, double width, double lineflux, double crtsig, double *fluxArray);

double getAtomicMass(int Z) {

    double atomM[30] = {1.00794, 4.002602, 6.941, 9.012182, 10.811,
                        12.0107, 14.0067, 15.9994, 18.9984032, 20.1797,
                        22.989770, 24.3050, 26.981538, 28.0855, 30.973761,
                        32.065, 35.4527, 39.948, 39.0983, 40.078,
                        44.955910, 47.867, 50.9415, 51.9961, 54.938049,
                        55.845, 58.933200, 58.6934, 63.456, 65.39};

    return atomM[Z-1];

}

void calc_lines(struct EMISSION *apec_data, int nbins, 
                          double *ebins, int hdu, int Z, int rmJ, 
                          double redshift,
                          double *spectrum) {
  /* note that hdu is indexed from 1, so 1 is the first one. This will
   * actually be #2 in the fits file, due to the parameters hdu.
   * */
   
  int ihdu;
  int iline;
  int ibin;
  int icont;
  struct EMISSION_BLOCK *emblock;
  double tmplineenergy;
  
  double width;
  double temperature;
  double mass;
  double z1, z15;
  double crtsig = 12.0;
  double velocity = 0.0;

  temperature = apec_data->kT[hdu-1];
  z1 = 1.0 + redshift;
  z15 = 0.5*z1;

  /* zero the spectrum */
  for (ibin=0;ibin<nbins;ibin++) {
    spectrum[ibin] = 0.0;
  }
  
  /* point emblock to the first block */
  emblock = apec_data->emdat;
  
  /* cycle through the blocks until we get to the desired one */
  for (ihdu=1;ihdu<hdu;ihdu++) {
    emblock = emblock->next;
  }
  
  /* add every line to the appropriate bin */
  if (rmJ==0) {
    /* this means do all the ions of element Z */
    for (iline=0;iline<emblock->nline;iline++) {
      tmplineenergy=keVtoA/emblock->lines[iline].lambda;
      //tmplineenergy /= z1;
    
      if (tmplineenergy < ebins[0]) continue;
      if (tmplineenergy > ebins[nbins]) continue;
      if (emblock->lines[iline].Z != Z) continue;
    
    /* find the appropriate bin */
      ibin = 0;
      while (ebins[ibin] < tmplineenergy) ibin++;
    
    /* the appropriate bin is ibin-1, add in the emissivity here */
      //spectrum[ibin-1] += emblock->lines[iline].epsilon; 
      
      width = tmplineenergy * (sqrt( temperature * keV_erg /
                            getAtomicMass(Z) / amu_g ) + velocity*1.0e5) / c_cm_s;  

      calcGaussianLine(ebins, nbins, tmplineenergy, width, emblock->lines[iline].epsilon, crtsig, spectrum);
    
    } 
  } else {
    for (iline=0;iline<emblock->nline;iline++) {
      tmplineenergy=keVtoA/emblock->lines[iline].lambda;
      //tmplineenergy /= z1;
    
      if (tmplineenergy < ebins[0]) continue;
      if (tmplineenergy > ebins[nbins]) continue;
      if (emblock->lines[iline].Z != Z) continue;
      if (emblock->lines[iline].rmJ != rmJ) continue;
    
    /* find the appropriate bin */
      ibin = 0;
      while (ebins[ibin] < tmplineenergy) ibin++;
    
    /* the appropriate bin is ibin-1, add in the emissivity here */
      //spectrum[ibin-1] += emblock->lines[iline].epsilon;

    
      width = tmplineenergy * (sqrt( temperature * keV_erg /
                            getAtomicMass(Z) / amu_g ) + velocity*1.0e5) / c_cm_s;

      calcGaussianLine(ebins, nbins, tmplineenergy, width, emblock->lines[iline].epsilon, crtsig, spectrum);

    
    } 
  }
};
