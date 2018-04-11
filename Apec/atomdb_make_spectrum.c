#include "readapec_spectrum.h"
#include <stdlib.h>
#include <stdio.h>


void calc_continuum(struct EMISSION *apec_data,
                    int nbins, double *ebins,  int hdu, int Z, int rmJ, 
                    double *spectrum);
void atomdb_make_spectrum(struct EMISSION *apec_data, int nbins, 
                          double *ebins, int hdu, double *spectrum) {
  /* note that hdu is indexed from 1, so 1 is the first one. This will
   * actually be #2 in the fits file, due to the parameters hdu.
   * */
   
  int ihdu;
  int iline;
  int ibin;
  int icont;
  struct EMISSION_BLOCK *emblock;
  double tmplineenergy;
  struct COMPRESS_CONT *cont;
  double tmpcoco[nbins];
  // move to the right HDU
  
  /* point emblock to the first block */
  emblock = apec_data->emdat;
  
  /* cycle through the blocks until we get to the desired one */
  for (ihdu=1;ihdu<hdu;ihdu++) {
    emblock = emblock->next;
  }
  
  /* debug*/
  /* add every line to the appropriate bin */
  for (iline=0;iline<emblock->nline;iline++) {
    tmplineenergy=12.398425/emblock->lines[iline].lambda;
    /*
     * 
     *  check if outside the range completely */
    
    
    if (tmplineenergy < ebins[0]) continue;
    if (tmplineenergy > ebins[nbins]) continue;
    
    /* find the appropriate bin */
    ibin = 0;
    while (ebins[ibin] < tmplineenergy) ibin++;
    
    /* the appropriate bin is ibin-1, add in the emissivity here */
    spectrum[ibin-1] += emblock->lines[iline].epsilon;
  }
     
  /* now do the continuum */
  cont = emblock->coco;
  
  for (icont=0;icont<emblock->nelement;icont++) {
    calc_continuum(apec_data,
                   nbins, ebins, hdu, cont->Z, 0, 
                   tmpcoco);
    cont=cont->next;
    for (ibin=0;ibin<nbins;ibin++){
      spectrum[ibin]+=tmpcoco[ibin];
    }

  } 
  
  
  
};
