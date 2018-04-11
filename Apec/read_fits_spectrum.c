#include <stdarg.h>
#include "readapec_spectrum.h"
#include <stdlib.h>
#include <stdio.h>

#include "messages.h"
void read_line_data(fitsfile *fptr, struct EMISSION_LINES *lines, int nrows);
void read_continuum_data(fitsfile *fptr, struct COMPRESS_CONT *coco, int nelements);
void read_parameters(fitsfile *fptrline, fitsfile *fptrcoco, 
                     struct EMISSION *emission, int Nhdu);



static struct EMISSION *apec_data_static;


/*static VariableType static_variable;
const VariableType * const global_variable_pointer = & static_variable;
*/
const struct EMISSION ** const apec_data_global_pointer = &apec_data_static;

unsigned long read_fits_spectrum(int require_checksum,
           char *linefile, char *cocofile,
           struct EMISSION **apec_data) {
 /* inputs: 
  * require_checksum: check the file's checksum is correct or exit with
  *                   error
  * linefile: fits linefile
  * cocofile: fits cocofile
  * 
  * outputs:
  * lev_dat_ptr
  * 
  */

  fitsfile *fptrline,*fptrcoco;  
  int Nhdu, ihdu;
  char value[MAXSTRLEN];
  char comment[MAXSTRLEN];
  int status = 0;
  int *NULLINT = 0;
  int hdutype = 0;
  struct EMISSION *emission;
//  struct EMISSION_BLOCK *emblock;
  struct EMISSION_BLOCK *emblock_ptr=NULL;
  struct EMISSION_LINES *lines_ptr;
  struct EMISSION_BLOCK *emblock_test=NULL;
  struct COMPRESS_CONT *coco_ptr;
  
  /* First, check that the is file OK and return the checksum and datasum */
  /* if (require_checksum) datasum = check_fits_file(fptr); */

  /* Figure out how many values we have to load. */
  /*fits_read_key(fptr, TINT, "NAXIS2", &Nrows, comment, &status);*/
  
  if (fits_open_file(&fptrline, linefile, READONLY, &status)) {
    fitsmess("fits_open_file",status,"Error opening file %s.",*linefile);
  }

  if (fits_open_file(&fptrcoco, cocofile, READONLY, &status)) {
    fitsmess("fits_open_file",status,"Error opening file %s.",*cocofile);
  }


  fits_get_num_hdus(fptrline, &Nhdu, &status);
  /* There are Nhdu-2 useful HDUs */
  Nhdu -= 2;
  
  /* now move to parameters HDU in line & coco file*/
  fits_movnam_hdu(fptrline, ANY_HDU, "PARAMETERS", 0, &status);
  fits_movnam_hdu(fptrcoco, ANY_HDU, "PARAMETERS", 0, &status);
  
  /* claim memory for the emission block, which is the master structure*/
  emission = (struct EMISSION *) malloc_safe(sizeof(struct EMISSION));
  emission->nblock = Nhdu;
  
  /* declare the temperature, density, time, nelement, ncont, npseudo & nline */
  emission->kT = (double *) malloc_safe(sizeof(double)*Nhdu);
  emission->eDensity = (double *) malloc_safe(sizeof(double)*Nhdu);
  emission->time = (double *) malloc_safe(sizeof(double)*Nhdu);
  emission->nelement = (int *) malloc_safe(sizeof(int)*Nhdu);
  emission->ncont = (int *) malloc_safe(sizeof(int)*Nhdu);
  emission->npseudo = (int *) malloc_safe(sizeof(int)*Nhdu);
  emission->nline = (int *) malloc_safe(sizeof(int)*Nhdu);
  
  apec_data_static = emission;
  
  /*get the data*/
  read_parameters(fptrline, fptrcoco, emission, Nhdu);

  
  fits_get_hdu_num(fptrcoco, &ihdu);
  /* now go through the remaining Nhdu and populate with line data */
  
  for (ihdu=0;ihdu<Nhdu;ihdu++) {
    /* move forward 1 HDU */
    fits_movrel_hdu(fptrline, 1, &hdutype, &status);
    fits_movrel_hdu(fptrcoco, 1, &hdutype, &status);
    /* check that we are in the right kind of HDU */
    fits_read_keyword(fptrline, "EXTNAME", value, comment,&status);
    
    
    
    if (ihdu==0){
      /* special case for 1st block in terms of pointer handling*/
      
       
      
      emission->emdat = (struct EMISSION_BLOCK *) 
                   malloc_safe(sizeof(struct EMISSION_BLOCK));
      emblock_ptr = emission->emdat;
    } else {
      /* set emblockptr and emblock*/
      
      emblock_ptr->next = (struct EMISSION_BLOCK *) 
                    malloc_safe(sizeof(struct EMISSION_BLOCK));
      emblock_ptr = emblock_ptr->next;
      
    }
    /* set so there is no "next" block for the moment */
    emblock_ptr->next = NULL;
    
    
    /* read in numbers from header */
    fits_read_key(fptrline, TDOUBLE, "HIERARCH TEMPERATURE", &emblock_ptr->kT, comment, &status); /* temperature of this block, keV */
    fits_read_key(fptrline, TDOUBLE, "DENSITY", &emblock_ptr->Ne, comment, &status); /* temperature of this block, keV */
    fits_read_key(fptrline, TDOUBLE, "TIME", &emblock_ptr->time, comment, &status); /* temperature of this block, keV */
    fits_read_key(fptrline, TINT, "N_LINES", &emblock_ptr->nline, comment, &status); /* temperature of this block, keV */
  
    
//    emblock_ptr->nline=ihdu;
    /* malloc some space for the lines */
    //*lines_ptr=NULL;
    /* malloc space for a new emission lines block */
    lines_ptr = (struct EMISSION_LINES *)
                         malloc_safe(emblock_ptr->nline*sizeof(struct EMISSION_LINES));
    /* set pointer for line data to point to this emission_lines block */
    emblock_ptr->lines=lines_ptr;
    read_line_data(fptrline, emblock_ptr->lines, emblock_ptr->nline);

    coco_ptr = (struct COMPRESS_CONT *)
                         malloc_safe(sizeof(struct COMPRESS_CONT));
    fits_read_keyword(fptrcoco, "EXTNAME", value, comment,&status);
                         
    /* set pointer for coco data to point to this COMPRESS_CONT block */
    emblock_ptr->coco=coco_ptr;
    
    read_continuum_data(fptrcoco, emblock_ptr->coco,*emission->nelement);
    emblock_ptr->nelement=*emission->nelement;
    
    
    
    
  }
  /* so now we have ALL THE DATAS*/
  
  /* some tests... */
  
  emblock_test = emission->emdat;
  
  for (ihdu=0;ihdu<Nhdu;ihdu++) {
    emblock_test = emblock_test->next;
  }
  
  
  //
   return status;
};





unsigned long read_fits_spectrum_external(int require_checksum,
           char *linefile, char *cocofile) {
 /* inputs: 
  * require_checksum: check the file's checksum is correct or exit with
  *                   error
  * linefile: fits linefile
  * cocofile: fits cocofile
  * 
  * outputs:
  * lev_dat_ptr
  * 
  */

  fitsfile *fptrline,*fptrcoco;  
  int Nhdu, ihdu;
  char value[MAXSTRLEN];
  char comment[MAXSTRLEN];
  int status = 0;
  int *NULLINT = 0;
  int hdutype = 0;
  struct EMISSION *emission;
//  struct EMISSION_BLOCK *emblock;
  struct EMISSION_BLOCK *emblock_ptr=NULL;
  struct EMISSION_LINES *lines_ptr;
  struct EMISSION_BLOCK *emblock_test=NULL;
  struct COMPRESS_CONT *coco_ptr;
  
  /* First, check that the is file OK and return the checksum and datasum */
  /* if (require_checksum) datasum = check_fits_file(fptr); */

  /* Figure out how many values we have to load. */
  /*fits_read_key(fptr, TINT, "NAXIS2", &Nrows, comment, &status);*/
  
  if (fits_open_file(&fptrline, linefile, READONLY, &status)) {
    fitsmess("fits_open_file",status,"Error opening file %s.",*linefile);
  }

  if (fits_open_file(&fptrcoco, cocofile, READONLY, &status)) {
    fitsmess("fits_open_file",status,"Error opening file %s.",*cocofile);
  }


  fits_get_num_hdus(fptrline, &Nhdu, &status);
  /* There are Nhdu-2 useful HDUs */
  Nhdu -= 2;
  
  /* now move to parameters HDU in line & coco file*/
  fits_movnam_hdu(fptrline, ANY_HDU, "PARAMETERS", 0, &status);
  fits_movnam_hdu(fptrcoco, ANY_HDU, "PARAMETERS", 0, &status);
  
  /* claim memory for the emission block, which is the master structure*/
  emission = (struct EMISSION *) malloc_safe(sizeof(struct EMISSION));
  emission->nblock = Nhdu;
  
  /* declare the temperature, density, time, nelement, ncont, npseudo & nline */
  emission->kT = (double *) malloc_safe(sizeof(double)*Nhdu);
  emission->eDensity = (double *) malloc_safe(sizeof(double)*Nhdu);
  emission->time = (double *) malloc_safe(sizeof(double)*Nhdu);
  emission->nelement = (int *) malloc_safe(sizeof(int)*Nhdu);
  emission->ncont = (int *) malloc_safe(sizeof(int)*Nhdu);
  emission->npseudo = (int *) malloc_safe(sizeof(int)*Nhdu);
  emission->nline = (int *) malloc_safe(sizeof(int)*Nhdu);
  
  apec_data_static = emission;
  
  /*get the data*/
  read_parameters(fptrline, fptrcoco, emission, Nhdu);

  
  fits_get_hdu_num(fptrcoco, &ihdu);
  /* now go through the remaining Nhdu and populate with line data */
  
  for (ihdu=0;ihdu<Nhdu;ihdu++) {
    /* move forward 1 HDU */
    fits_movrel_hdu(fptrline, 1, &hdutype, &status);
    fits_movrel_hdu(fptrcoco, 1, &hdutype, &status);
    /* check that we are in the right kind of HDU */
    fits_read_keyword(fptrline, "EXTNAME", value, comment,&status);
    
    
    
    if (ihdu==0){
      /* special case for 1st block in terms of pointer handling*/
      
       
      
      emission->emdat = (struct EMISSION_BLOCK *) 
                   malloc_safe(sizeof(struct EMISSION_BLOCK));
      emblock_ptr = emission->emdat;
    } else {
      /* set emblockptr and emblock*/
      
      emblock_ptr->next = (struct EMISSION_BLOCK *) 
                    malloc_safe(sizeof(struct EMISSION_BLOCK));
      emblock_ptr = emblock_ptr->next;
      
    }
    /* set so there is no "next" block for the moment */
    emblock_ptr->next = NULL;
    
    
    /* read in numbers from header */
    fits_read_key(fptrline, TDOUBLE, "HIERARCH TEMPERATURE", &emblock_ptr->kT, comment, &status); /* temperature of this block, keV */
    fits_read_key(fptrline, TDOUBLE, "DENSITY", &emblock_ptr->Ne, comment, &status); /* temperature of this block, keV */
    fits_read_key(fptrline, TDOUBLE, "TIME", &emblock_ptr->time, comment, &status); /* temperature of this block, keV */
    fits_read_key(fptrline, TINT, "N_LINES", &emblock_ptr->nline, comment, &status); /* temperature of this block, keV */
  
    
//    emblock_ptr->nline=ihdu;
    /* malloc some space for the lines */
    //*lines_ptr=NULL;
    /* malloc space for a new emission lines block */
    lines_ptr = (struct EMISSION_LINES *)
                         malloc_safe(emblock_ptr->nline*sizeof(struct EMISSION_LINES));
    /* set pointer for line data to point to this emission_lines block */
    emblock_ptr->lines=lines_ptr;
    read_line_data(fptrline, emblock_ptr->lines, emblock_ptr->nline);
    
    coco_ptr = (struct COMPRESS_CONT *)
                         malloc_safe(sizeof(struct COMPRESS_CONT));
    fits_read_keyword(fptrcoco, "EXTNAME", value, comment,&status);
                         
    /* set pointer for coco data to point to this COMPRESS_CONT block */
    emblock_ptr->coco=coco_ptr;
    
    read_continuum_data(fptrcoco, emblock_ptr->coco,*emission->nelement);
    
    emblock_ptr->nelement=*emission->nelement;
    
    
    
    
  }
  /* so now we have ALL THE DATAS*/
  
   return status;
};
