/* this file contains all the publicly accessible routines in
 * libreadapec */

 /* version 0.1, started by Adam Foster 29-Oct-2013 */
 /* version 0.2, Adam Foster 01-Apr-2014 */

#include "readapec_spectrum.h"
#include <math.h>
#include <float.h>
#include "readapec.h"


void readapec_getdata(char *linefile,
                       char *cocofile,
                       struct EMISSION **apec_data) {

  /* This is just a wrapper for read_fits_spectrum, which reads in
   * all of the apec outputs
   *
   * INPUTS
   * linefile     filename of _line.fits apec file
   * cocofile     filename of _comp.fits or _coco.fits apec continuum
   *              file
   *
   * OUTPUTS
   * apec_data    all of the emissivity data
   * 
   * Version 0.1: original version
   *              ARF 04-Nov-2013
   * */

   read_fits_spectrum(1, linefile, cocofile, apec_data);

};


void readapec_getdata_extern(char *linefile,
                       char *cocofile) {

  /* This is just a wrapper for read_fits_spectrum_extrnal. It reads in all the data
   * from the emissivity files and puts it into a global variable,apec_data_memstore
   *
   * INPUTS
   * linefile     filename of _line.fits apec file
   * cocofile     filename of _comp.fits or _coco.fits apec continuum
   *              file
   *
   * OUTPUTS
   * none (though apec_data_memstore now contains emissivity data)
   * 
   * Version 0.1: original version
   *              ARF 22-Sep-2014
   * */
   
   read_fits_spectrum_external(1, linefile, cocofile);

};


void readapec_simple_spectrum(const struct EMISSION *apec_data,
                          int nbins,
                          double *ebins,
                          int hdu,
                          double *spectrum) {
  /* This routine calculates the spectrum for a specific HDU.
   * It does not separate the emissivity by ion or element,
   * and assumes an abundance of 1.0 for all elements.
   * It is the quick & dirty version
   *
   * INPUTS
   * apec_data    all of the emissivity data (returned by
   *              read_apec_outputs)
   * nbins        the number of bins the output spectrum
   * *ebins       the edges of the bins, in keV. Note that the array
   *              length is nbins+1. Must be monotonically increasing.
   * hdu          Which HDU in the fits file we are interested. Note
   *              that they are indexed from 1, with 1 being the first
   *              "EMISSIVITY" HDU. This is typically HDU 2 in most
   *              FITS files.
   *
   * OUTPUTS
   * *spectrum    The emissivity in ph cm3 s-1 bin-1
   * 
   * Version 0.1: original version
   *              ARF 04-Nov-2013
   * */

  /* Sort out the emissivity data. 
   * The files are obtained by a call to readapec_getdata_extern
   * or readapec_getdata. 
   * 
   * If using the former, a global variable 
   * apec_data_memstore will be populated, and apec_data should be NULL
   * For the latter, the returned EMISSION structure should be supplied
   * as apec_data.*/
   
  if (apec_data==NULL) {
    /* if no data supplied, we need to see if there is preloaded data */
    if (apec_data_memstore == NULL) {
      /* nope = problem */
      printf ("\n***ERROR*** %s: %s Line %d:\n", __func__, __FILE__, __LINE__);
      printf ("   must load emissivity data first (see readapec_getdata or readapec_getdata_extern)\n");
    } else {
      /* yes = use this */
      apec_data=apec_data_memstore;
    }
  }

  atomdb_make_spectrum(apec_data, nbins,
                      ebins, hdu, spectrum);
};

void readapec_calc_hdu_elem_spectrum(const struct EMISSION *apec_data,
                          int nbins,
                          double *ebins,
                          int hdu,
                          int Z,
                          double *spectrum) {


  /* This routine calculates the spectrum of one element for a specific 
   * HDU. It assumes an abundance of 1.0 solar for all elements.
   * It is the quick & dirty version
   *
   * INPUTS
   * apec_data    all of the emissivity data (returned by
   *              read_apec_outputs)
   * nbins        the number of bins the output spectrum
   * *ebins       the edges of the bins, in keV. Note that the array
   *              length is nbins+1. Must be monotonically increasing.
   * hdu          Which HDU in the fits file we are interested. Note
   *              that they are indexed from 1, with 1 being the first
   *              "EMISSIVITY" HDU. This is typically HDU 2 in most
   *              FITS files.
   * Z            Atomic number of element of interest (e.g. 6=carbon)
   *
   * OUTPUTS
   * *spectrum    The emissivity in ph cm3 s-1 bin-1
   * 
   * Version 0.1: original version
   *              ARF 04-Nov-2013
   */

  int i;
  /* Sort out the emissivity data. 
   * The files are obtained by a call to readapec_getdata_extern
   * or readapec_getdata. 
   * 
   * If using the former, a global variable 
   * apec_data_memstore will be populated, and apec_data should be NULL
   * For the latter, the returned EMISSION structure should be supplied
   * as apec_data.*/
   
  if (apec_data==NULL) {
    /* if no data supplied, we need to see if there is preloaded data */
    if (apec_data_memstore == NULL) {
      /* nope = problem */
      printf ("\n***ERROR*** %s: %s Line %d:\n", __func__, __FILE__, __LINE__);
      printf ("   must load emissivity data first (see readapec_getdata or readapec_getdata_extern)\n");
    } else {
      /* yes = use this */
      apec_data=apec_data_memstore;
    }
  }

  /* zero the data */
  for (i=0;i<nbins; i++) {
    spectrum[i]=0.0;
  }

  /* get the line data */
  calc_lines(apec_data, nbins, ebins, hdu, Z, 0, spectrum);
  /* get the continuum data */
  calc_continuum(apec_data, nbins, ebins, hdu, Z, 0, spectrum);

  /* and we're done */
};

void readapec_calc_hdu_ion_spectrum(const struct EMISSION *apec_data,
                          int nbins,
                          double *ebins,
                          int hdu,
                          int Z,
                          int rmJ,
                          double *spectrum) {


  /* This routine calculates the spectrum of an ion for a specific HDU.
   * It assumes an abundance of 1.0 for all elements.
   *
   * INPUTS
   * apec_data    all of the emissivity data (returned by
   *              read_apec_outputs)
   * nbins        the number of bins the output spectrum
   * *ebins       the edges of the bins, in keV. Note that the array
   *              length is nbins+1. Must be monotonically increasing.
   * hdu          Which HDU in the fits file we are interested. Note
   *              that they are indexed from 1, with 1 being the first
   *              "EMISSIVITY" HDU. This is typically HDU 2 in most
   *              FITS files.
   * Z            Atomic number of element of interest (e.g. 6=carbon)
   * rmJ          The ionization stage of interest (1=neutral)
   *
   * OUTPUTS
   * *spectrum    The emissivity in ph cm3 s-1 bin-1
   * 
   * Version 0.1: original version
   *              ARF 04-Nov-2013
   * */

  int i;
  double lspec[nbins], cspec[nbins];
  int doline, docont;
  doline=1;
  docont=1;
  
  readapec_calc_hdu_ion_spectrum_part(apec_data,
                          nbins,
                          ebins,
                          hdu,
                          Z,
                          rmJ,
                          doline,
                          docont,
                          spectrum);
  
  
  /* and we're done */
};


void readapec_calc_hdu_ion_spectrum_part(const struct EMISSION *apec_data,
                          int nbins,
                          double *ebins,
                          int hdu,
                          int Z,
                          int rmJ,
                          int doline,
                          int docont,
                          double *spectrum) {


  /* This routine calculates the spectrum of an ion for a specific HDU.
   * It assumes an abundance of 1.0 for all elements.
   *
   * INPUTS
   * apec_data    all of the emissivity data (returned by
   *              read_apec_outputs)
   * nbins        the number of bins the output spectrum
   * *ebins       the edges of the bins, in keV. Note that the array
   *              length is nbins+1. Must be monotonically increasing.
   * hdu          Which HDU in the fits file we are interested. Note
   *              that they are indexed from 1, with 1 being the first
   *              "EMISSIVITY" HDU. This is typically HDU 2 in most
   *              FITS files.
   * Z            Atomic number of element of interest (e.g. 6=carbon)
   * rmJ          The ionization stage of interest (1=neutral)
   * doline       If != 0, include line data in spectrum
   * docont       If != 0, include continuum data in spectrum
   *
   * OUTPUTS
   * *spectrum    The emissivity in ph cm3 s-1 bin-1
   * 
   * Version 0.1: original version
   *              ARF 04-Nov-2013
   * */

  int i;
  double lspec[nbins], cspec[nbins];

  /* Sort out the emissivity data. 
   * The files are obtained by a call to readapec_getdata_extern
   * or readapec_getdata. 
   * 
   * If using the former, a global variable 
   * apec_data_memstore will be populated, and apec_data should be NULL
   * For the latter, the returned EMISSION structure should be supplied
   * as apec_data.*/
   
  if (apec_data==NULL) {
    /* if no data supplied, we need to see if there is preloaded data */
    if (apec_data_memstore == NULL) {
      /* nope = problem */
      printf ("\n***ERROR*** %s: %s Line %d:\n", __func__, __FILE__, __LINE__);
      printf ("   must load emissivity data first (see readapec_getdata or readapec_getdata_extern)\n");
    } else {
      /* yes = use this */
      apec_data=apec_data_memstore;
    }
  }

  /* zero the data */
  for (i=0;i<nbins; i++) {
    spectrum[i]=0.0;
  }
  
  /* get the line data */
  if (doline != 0) {
    calc_lines(apec_data, nbins, ebins, hdu, Z, rmJ, lspec);
  } else {
    for (i=0; i<nbins; i++) {
      lspec[i] = 0.0;
    }
  }
  
  /* get the continuum data */
  if (docont != 0) {
    calc_continuum(apec_data, nbins, ebins, hdu, Z, rmJ, cspec);
  } else {
    for (i=0; i<nbins; i++) {
      cspec[i] = 0.0;
    }
  } 
  
  for (i=0; i<nbins; i++) {
    spectrum[i] = cspec[i]+lspec[i];
  }
  
  /* and we're done */
};



int readapec_find_kT_hdu(const struct EMISSION *apec_data,
                double kT, int log) {


  /* This routine finds the nearest HDU index for a given temperature
   *
   * INPUTS
   * apec_data    all of the emissivity data (returned by
   *              read_apec_outputs)
   * kT           target temperature (keV)
   * log          if log=1, then find closest on a log10 scale
   *
   * RETURNS
   *      The hdu nearest to the parameter of interest (in this case kT)
   * 
   * Version 0.1: original version
   *              ARF 04-Nov-2013
   * */

  int i, dolog;
  double deltakt;
  double bestdeltakt;
  int closestind;
  double tmpkT;

  /* Sort out the emissivity data. 
   * The files are obtained by a call to readapec_getdata_extern
   * or readapec_getdata. 
   * 
   * If using the former, a global variable 
   * apec_data_memstore will be populated, and apec_data should be NULL
   * For the latter, the returned EMISSION structure should be supplied
   * as apec_data.*/
   
  if (apec_data==NULL) {
    /* if no data supplied, we need to see if there is preloaded data */
    if (apec_data_memstore == NULL) {
      /* nope = problem */
      printf ("\n***ERROR*** %s: %s Line %d:\n", __func__, __FILE__, __LINE__);
      printf ("   must load emissivity data first (see readapec_getdata or readapec_getdata_extern)\n");
    } else {
      /* yes = use this */
      apec_data=apec_data_memstore;
    }
  }

  if (log==1) {
    dolog = 1;
  } else if (log==0) {
    dolog = 0;
  } else {
    dolog = 0;
  }
  closestind=-1;

  bestdeltakt=log10(1e20);

  if (dolog) {
    tmpkT = log10(kT);
  } else {
    tmpkT=kT;
  }

  if (dolog) {
    for (i=0;i<apec_data->nblock; i++) {
      deltakt = fabs(log10(apec_data->kT[i])-tmpkT);
      if (deltakt < bestdeltakt) {
        bestdeltakt = deltakt;
        closestind = i;
      }
    }
  } else {
    for (i=0;i<apec_data->nblock; i++) {
      deltakt = fabs(apec_data->kT[i]-tmpkT);
      if (deltakt < bestdeltakt) {
        bestdeltakt = deltakt;
        closestind = i;
      }
    }
  }
  /* return the index, adding 1 to index from 1 instead of zero */
  return(closestind+1);

};

void readapec_calc_ion_spectrum(const struct EMISSION *apec_data,
                       int nbins,
                       double *ebins,
                       double kT,
                       int Z,
                       int rmJ,
                       int nearest,
                       int doline,
                       int docont,
                       double *spectrum) {

  /* This routine calculates the spectrum at a specific temperature
   * Setting nearest = 1, it will calculate the spectrum of the nearest
   * HDU, setting it =0, it will interpolate between the 2 adjacent
   * HDUs in log space.
   *
   * INPUTS
   * apec_data    all of the emissivity data (returned by
   *              read_apec_outputs)
   * nbins        the number of bins the output spectrum
   * *ebins       the edges of the bins, in keV. Note that the array
   *              length is nbins+1. Must be monotonically increasing.
   * hdu          Which HDU in the fits file we are interested. Note
   *              that they are indexed from 1, with 1 being the first
   *              "EMISSIVITY" HDU. This is typically HDU 2 in most
   *              FITS files.
   * Z            Atomic number of element of interest (e.g. 6=carbon)
   * rmJ          The ionization stage of interest (1=neutral). If 0,
   *              do all the ions of the element.
   * nearest      If 0, interpolate to get exact results.
   *              If 1, just return the sepctrum of the nearest HDU
   * doline       1:Include line emissions 0: skip
   * docont       1:Include continuum emissions 0: skip
   * nearest      If 0, interpolate to get exact results.
   *              If 1, just return the sepctrum of the nearest HDU
   * OUTPUTS
   * *spectrum    The emissivity in ph cm3 s-1 bin-1
   * 
   * Version 0.1: original version
   *              ARF 04-Nov-2013
   * Version 0.2: Added doline and docont variables.
   *              ARF 30-Mar-2015
   * */

  double speclo[nbins], spechi[nbins];
  double interpfrac;
  double ktlo, kthi;
  int ihdulo, ihduhi, i;

  /* Sort out the emissivity data. 
   * The files are obtained by a call to readapec_getdata_extern
   * or readapec_getdata. 
   * 
   * If using the former, a global variable 
   * apec_data_memstore will be populated, and apec_data should be NULL
   * For the latter, the returned EMISSION structure should be supplied
   * as apec_data.*/
   
  if (apec_data==NULL) {
    /* if no data supplied, we need to see if there is preloaded data */
    if (apec_data_memstore == NULL) {
      /* nope = problem */
      printf ("\n***ERROR*** %s: %s Line %d:\n", __func__, __FILE__, __LINE__);
      printf ("   must load emissivity data first (see readapec_getdata or readapec_getdata_extern)\n");
    } else {
      /* yes = use this */
      apec_data=apec_data_memstore;
    }
  }
  
    /* get the nearest temperature HDU */
  ihdulo = readapec_find_kT_hdu(apec_data, kT, 1);

  if ((nearest==1) || (apec_data->kT[ihdulo-1] == kT)) {
    readapec_calc_hdu_ion_spectrum_part(apec_data,
                          nbins,
                          ebins,
                          ihdulo,
                          Z,
                          rmJ,
                          doline,
                          docont,
                          spectrum);
    return;

  }

  if (apec_data->kT[ihdulo-1] > kT) {
    /* in this case, this is actually ihduhi */
    ihduhi = ihdulo;
    ihdulo = ihdulo-1;
  } else  {
    ihduhi = ihdulo+1;
  }

  /* get the spectra */
  readapec_calc_hdu_ion_spectrum_part(apec_data,
                        nbins,
                        ebins,
                        ihdulo,
                        Z,
                        rmJ,
                        doline,
                        docont,
                        speclo);
  /* get the spectra */
  readapec_calc_hdu_ion_spectrum_part(apec_data,
                        nbins,
                        ebins,
                        ihduhi,
                        Z,
                        rmJ,
                        doline,
                        docont,
                        spechi);
  /* calculate the interpolation fraction */
  interpfrac = (log(kT)-log(apec_data->kT[ihdulo-1]))/
               (log(apec_data->kT[ihduhi-1])-log(apec_data->kT[ihdulo-1]));

  /* calculate the new values */
  for (i=0;i<nbins;i++) {
    if (spechi[i] == speclo[i]) {
      /* skip interpolation if they are identical */
      spectrum[i] = speclo[i];
    } else if ((spechi[i] < FLT_MIN) | (speclo[i] < FLT_MIN)) {
      /* If less than minimum float value, set to zero */
      spectrum[i] = 0.0;
    } else {
      /* Otherwise, we have "regular" data - interpolate */
      spectrum[i] = exp((log(spechi[i])-log(speclo[i]))*interpfrac
                        + log(speclo[i]));
    }
  }
  return;
};

void readapec_calc_allion_spectrum(const struct EMISSION *apec_data,
                          int nbins,
                          double *ebins,
                          double kT,
                          int nearest,
                          int nelements,
                          int *Zlist,
                          int byion,
                          struct EMISSION_LIST **emission_list) {
  /* This routine calculates the spectrum at a specific temperature
   * for all the ions of all the elements in Zlist. It then returns 
   * these in an EMISSIONLIST structure, which stores them as separate
   * emissivities. This is useful if you want to later change abundances
   * /
   * 
   * INPUTS
   * apec_data    all of the emissivity data (returned by
   *              readapec_getdata). NOTE that if used 
   *              readapec_getdata_extern this should be NULL.
   * nbins        the number of bins the output spectrum
   * *ebins       the edges of the bins, in keV. Note that the array
   *              length is nbins+1. Must be monotonically increasing.
   * kT           The temperature requested.
   * nearest      If 0, interpolate to get exact results.
   *              If 1, just return the sepctrum of the nearest HDU
   * nelements    Number of elements to get spectra for
   * Zlist        Atomic numbers of elements of interest (e.g. 6=carbon)
   * byion        If 0, return results for each element summed.
   *              If 1, return results for each ion separately.
   * OUTPUTS
   * emission_list  All of the emissivities in ph cm^3 s-1 bin-1
   * 
   * Version 0.1: original version
   *              ARF 04-Nov-2013
   *
   * Version 0.2: Corrected interpolation to handle low emissivity values:
   *              was previously returning Nans instead of zeros.
   *              ARF 01-Apr-2014
   * */
  
  int i, iZ, irmJ, ihdulo;
/*  struct EMISSION_ION emission_ion;*/
  struct EMISSION_ION *emission_ion_ptr;
  struct EMISSION_LIST *emission_list_tmp;
  double tmpspectrum[nbins];
  /* ok, let's do this */

/* *********************************************************************** */
//  emission_list_tmp = (struct EMISSION *) malloc_safe(sizeof(struct EMISSION));
//  emission->nblock = Nhdu;
// 
//  /* declare the temperature, density, time, nelement, ncont, npseudo & nline */
//  emission->kT = (double *) malloc_safe(sizeof(double)*Nhdu);
//  emission->eDensity = (double *) malloc_safe(sizeof(double)*Nhdu);
//  emission->time = (double *) malloc_safe(sizeof(double)*Nhdu);
//  emission->nelement = (int *) malloc_safe(sizeof(int)*Nhdu);
//  emission->ncont = (int *) malloc_safe(sizeof(int)*Nhdu);
//  emission->npseudo = (int *) malloc_safe(sizeof(int)*Nhdu);
//  emission->nline = (int *) malloc_safe(sizeof(int)*Nhdu);
// 
//  *apec_data = emission;
/* *********************************************************************** */

  /* Sort out the emissivity data. 
   * The files are obtained by a call to readapec_getdata_extern
   * or readapec_getdata. 
   * 
   * If using the former, a global variable 
   * apec_data_memstore will be populated, and apec_data should be NULL
   * For the latter, the returned EMISSION structure should be supplied
   * as apec_data.*/
   
  if (apec_data==NULL) {
    /* if no data supplied, we need to see if there is preloaded data */
    if (apec_data_memstore == NULL) {
      /* nope = problem */
      printf ("\n***ERROR*** %s: %s Line %d:\n", __func__, __FILE__, __LINE__);
      printf ("   must load emissivity data first (see readapec_getdata or readapec_getdata_extern)\n");
    } else {
      /* yes = use this */
      apec_data=apec_data_memstore;
    }
  }
  
  /* free up some memory for the first block */
  emission_list_tmp = (struct EMISSION_LIST *) malloc_safe(sizeof(struct EMISSION_LIST));
  /* populate the data */
  ihdulo = readapec_find_kT_hdu(apec_data, kT, 1);
  emission_list_tmp->kT = apec_data->kT[ihdulo-1];
  emission_list_tmp->eDensity = apec_data->eDensity[ihdulo-1];
  emission_list_tmp->time = apec_data->time[ihdulo-1];
  emission_list_tmp->nelement = nelements;
  emission_list_tmp->Zlist = (int *) malloc_safe(sizeof(int)*nelements);
  for (i=0;i<nelements;i++) {
    emission_list_tmp->Zlist[i]=Zlist[i];
  }
  emission_list_tmp->ion_emission=NULL;
  emission_list_tmp->nbins = nbins;
  *emission_list = emission_list_tmp;
  
  for (iZ=0; iZ < nelements;iZ++) {
    if (byion) {
      for (irmJ=1;irmJ<=Zlist[iZ];irmJ++) {
        /* do this by ion */
        for (i=0;i<nbins;i++) {
          tmpspectrum[i]=0.0;
        }
        readapec_calc_ion_spectrum(apec_data,
                          nbins,
                          ebins,
                          kT,
                          Zlist[iZ],
                          irmJ,
                          nearest,
                          1,
                          1,
                          tmpspectrum);
        
        /* now allocate this to a new emission_ion block */

        if (emission_list_tmp->ion_emission==NULL) {
          emission_list_tmp->ion_emission =  (struct EMISSION_ION *)
              malloc_safe(sizeof(struct EMISSION_ION));
          emission_ion_ptr = emission_list_tmp->ion_emission;
        } else {
          emission_ion_ptr->next = (struct EMISSION_ION *)
               malloc_safe(sizeof(struct EMISSION_ION));
          emission_ion_ptr = emission_ion_ptr->next;
        }
        emission_ion_ptr->next = NULL;
        
        /* put the data in emission_ion_ptr */
        emission_ion_ptr->Z = Zlist[iZ];
        emission_ion_ptr->rmJ = irmJ;
        for (i=0;i<nbins;i++) {
          emission_ion_ptr->emissivity[i]=tmpspectrum[i];
        }
          
      }
    } else {
      for (i=0;i<nbins;i++) {
        tmpspectrum[i]=0.0;
      }
      irmJ = 0;
      readapec_calc_ion_spectrum(apec_data,
                        nbins,
                        ebins,
                        kT,
                        Zlist[iZ],
                        irmJ,
                        nearest,
                        1,
                        1,
                        tmpspectrum);
      
      /* now allocate this to a new emission_ion block */
      if (emission_list_tmp->ion_emission==NULL) {
        emission_list_tmp->ion_emission = (struct EMISSION_ION *) 
                               malloc_safe(sizeof(struct EMISSION_ION));
        emission_ion_ptr = emission_list_tmp->ion_emission;
      } else {
        emission_ion_ptr->next = (struct EMISSION_ION *) 
                               malloc_safe(sizeof(struct EMISSION_ION));
        emission_ion_ptr = emission_ion_ptr->next;
      }
      emission_ion_ptr->next = NULL;
        
        /* put the data in emission_ion_ptr */
      emission_ion_ptr->Z = Zlist[iZ];
      emission_ion_ptr->rmJ = irmJ;
      /* malloc some space for the new emissivity array */
      emission_ion_ptr->emissivity = (double *) malloc_safe(sizeof(double) * nbins);
      for (i=0;i<nbins;i++) {
        emission_ion_ptr->emissivity[i]=tmpspectrum[i];
      }
    }
  }
};

void readapec_calc_total_emission_abundance(struct EMISSION_LIST *emission_list,
                                   int nelements,
                                   int *Zlist,
                                   double *abundlist,
                                   double defaultabund,
                                   double *spectrum) {

   /* This routine takes a calculated emission list (separated by ion
    * or element) and applies a relative abundance to them. There are
    * already abundances for the solar photosphere from Andres and 
    * Grevesse 1989 built into the emissivities, so new abundances
    * should be entered relative to these values */
    
  /* 
   * INPUTS
   * emission_list  All of the emissivities in ph cm^3 s-1 bin-1
   * nelements      The number of elements for which abundances will be
   *                provided. Any elements which are in "emission_list"
   *                but not specified in Zlist will be assumed to have
   *                abundance as specified in unlistedabund
   * Zlist          the elements for which abundances are provided
   * abundlist      the abundances of these elements
   * defaultabund  The default abundances for elements not in abundlist
   * OUTPUTS
   * spectrum       The total spectrum
   * 
   * Version 0.1: original version
   *              ARF 04-Nov-2013
   * */
   
   
   /* go through every block in the emission_list */
   
   struct EMISSION_ION *ion_emission;
   int Z, i, ind;
   double abfrac;
   
   ion_emission = emission_list->ion_emission;
   /* zero the return spectrum */
   for (i=0; i< emission_list->nbins; i++) {
     spectrum[i] = 0.0;
   }
   ind = -1;
   /* cycle through each block */
   while (ion_emission != NULL) {
     
     Z = ion_emission->Z;
     /* find matching abundance */
     for (i=0; i<nelements; i++) {
       if (Zlist[i] == Z) {
         ind = i;
       }
     }
     /* if no match, use default abund */
     if (ind < 0) {
       abfrac = defaultabund;
     } else {
       abfrac = abundlist[ind];
     }
     
     /* modification for electron-electron bremsstrahlung:
      * since Ne = 1.2NH, multiply by 1.2. This can/will be
      * replaced by full Ne calculation at a later date*/
      
     if (Z==0) {
       abfrac *= 1.2;
     }

     if (abfrac > 0.0) {
       for (i=0; i < emission_list->nbins; i++) {
         spectrum[i] += ion_emission->emissivity[i]*abfrac;
       }
     }
     
     if (ion_emission->emissivity != NULL) free(ion_emission->emissivity);
       
     ion_emission = ion_emission->next;
   }
 
   return;
   
};










