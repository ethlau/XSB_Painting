#include "readapec_spectrum.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void calc_continuum_process(struct COMPRESS_CONT *cont, int nbins, 
                            double *ebins, int ncont, double *econt,
                            double *ccont, double redshift, double *spectrum);


void calc_continuum(struct EMISSION *apec_data,
                    int nbins, double *ebins,  int hdu, int Z, int rmJ, 
                    double redshift,
                    double *spectrum){
  int nallbins;
  struct EMISSION_BLOCK *emblock;
  struct COMPRESS_CONT *cont;
  int ibin, ihdu;

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

  /* set cont to point to the first continuum in this hdu*/  
  cont = emblock->coco;
  
  /* now find the right continuum  for this Z & rmJ */
  
  while ((cont->Z != Z) || (cont->rmJ != rmJ)) {
    cont=cont->next;
    if (cont == NULL) {
      printf("No matching continuum data found for hdu=%i, Z=%i, rmJ=%i\n",
             hdu,Z,rmJ);
      return;
    }
  }
  
  /* get the continuum */
  calc_continuum_process(cont, nbins, ebins, cont->Ncont, cont->E_cont,
                         cont->cont, redshift, spectrum);
  /* get the pseudocontinuum */
  calc_continuum_process(cont, nbins, ebins, cont->Npseudo, cont->E_pseudo,
                         cont->pseudo, redshift, spectrum);

  /* and we're done */

};


void calc_continuum_process(struct COMPRESS_CONT *cont, int nbins, 
                            double *ebins, int ncont, double *econt,
                            double *ccont, double redshift, double *spectrum){

  /* first, get the number of bins in total */
  int nallbins;
  double *eallbins;
  double *wallbins;
  double *callbins;
  double *callbinedges;
  int *indbins;
  int itmpcont, itmpbins, blank;
  int i,j;
  
  double z1, z15;

  z1 = (1.0+redshift);
  z15 = 0.5*z1;

  /* Redshift energy bin in emitter's frame */
  //for (i = 0; i < ncont; i++) {
  //    econt[i] /= z1;
  //}
  
  nallbins = nbins + ncont;
  /* malloc a big pile of temporary arrays */  
  eallbins = (double *) malloc_safe(sizeof(double)*(nallbins+1));
  wallbins = (double *) malloc_safe(sizeof(double)*(nallbins));
  callbins = (double *) malloc_safe(sizeof(double)*(nallbins));
  callbinedges = (double *) malloc_safe(sizeof(double)*(nallbins+1));
  indbins = (int *) malloc_safe(sizeof(int)*(nbins+1));
  itmpcont=0;
  itmpbins=0;
  
/*  for (i=0;i<nbins;i++){
    spectrum[i]=0.0;
  }*/
  for (i=0; i<nallbins+1;i++) {
    if (itmpcont >= ncont) {
      /* we are at the end of the continuum bins already */
      eallbins[i] = ebins[itmpbins];
      indbins[itmpbins] = i;
      itmpbins++;
    } else if (itmpbins > nbins) {
      eallbins[i] = econt[itmpcont];
      itmpcont++;
    } else {
      if (econt[itmpcont] > ebins[itmpbins]) {
        eallbins[i] = ebins[itmpbins];
        indbins[itmpbins] = i;
        itmpbins++;
      } else {
        eallbins[i] = econt[itmpcont];
        itmpcont++;
      }
    }
  }
  
  /* get widths */
  for (i=0; i<nallbins;i++) {
    wallbins[i]=eallbins[i+1]-eallbins[i];
  }
  
  /* exception case: if there are only 2 bins, and they are all zeros, then pass */
  for (i=0;i<nallbins;i++){
    callbins[i]=0.0;
  }
  
  blank = 0;
  
  if (ncont == 2) {
    if ((ccont[0]==0.0) & (ccont[1]==0.0)) {
      blank = 1;
    }
  }
  
  
  if (blank==0) {
    itmpbins=0;
    itmpcont=0;
    for (i=0;i<nallbins+1;i++) {
      if (itmpbins >=nbins){
         /* this means we have a match: no need to interpolate */
        callbinedges[i] = ccont[itmpcont];
        itmpcont++;     
      }  else if (indbins[itmpbins]!=i) {
        /* this means we have a match: no need to interpolate */
        callbinedges[i] = ccont[itmpcont];
        itmpcont++;
      } else {
/*        callbinedges[i] = exp(((log(eallbins[i])-log(econt[itmpcont-1]))/
                       (log(econt[itmpcont])-log(econt[itmpcont-1])))*
                       (log(ccont[itmpcont])-log(ccont[itmpcont-1]))+
                       log(ccont[itmpcont-1]));*/
                       
        callbinedges[i] =(eallbins[i]-econt[itmpcont-1])/
                       (econt[itmpcont]-econt[itmpcont-1])*
                       (ccont[itmpcont]-ccont[itmpcont-1])+
                       ccont[itmpcont-1];
                       
        itmpbins++;               
      }
    }
    
  /* ok, so now we have the emissivity at each bin edge. 
   * Calculate the flux in each bin*/
    for (i=0;i<nallbins;i++) {
      callbins[i] = (callbinedges[i]+callbinedges[i+1])*wallbins[i]/2.0;
    }
  
  /* sum bins appropriately into spectrum*/
  /* this will keep track of which bin in the spectrum to
               * put the data in */
               
    for (i=0;i<nbins;i++) {
      for (j=indbins[i];j<indbins[i+1];j++) {
        spectrum[i] += callbins[j];
      }
    }
  
  }
  free(eallbins);
  free(wallbins); 
  free(callbins);
  free(callbinedges);
  free(indbins);
  return;
};
   
