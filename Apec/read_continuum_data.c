#include "readapec_spectrum.h"
void read_continuum_data(fitsfile *fptr, struct COMPRESS_CONT *coco, int nelements){

  int iCol, iRow, anynul,i;
  int status=0;
  float var_float, fnul_val;
  int var_int;
  char *colnames[8] = {"Z", "rmJ", 
                       "N_Cont", "E_Cont", "Continuum", 
                       "N_Pseudo", "E_Pseudo","Pseudo"};
                       
  int colnum[8]; /* 6 names listed above */
  int Ncols = 8;
  char comment[MAXSTRLEN];
  struct COMPRESS_CONT *coco_ptr;
  float arr_float[MAXARRFLOAT];
  unsigned long datasum = 0;
  
  /* Get the column numbers */
  for(iCol=0;iCol<Ncols;iCol++) {
    fits_get_colnum(fptr, CASEINSEN, colnames[iCol], &(colnum[iCol]), &status);
    if (status!=0) {
      fitsmess("read_fits_datafiles",status,"Error reading %s column number",
         colnames[iCol]);
    }
  }

  /* read in the data */
  for (iRow=1;iRow<nelements+1;iRow++) {
    /* move coco pointer to a fresh structure */
    if (iRow > 1) {
      coco->next =  (struct COMPRESS_CONT *) 
                    malloc_safe(sizeof(struct COMPRESS_CONT));    
      coco = coco->next;
      coco->next = NULL;
    } else {
      coco->next = NULL;
    }
    /* get the data */
    fits_read_col(fptr, TINT, colnum[0], iRow, 1, 1, 0, &var_int, &anynul,
      &status);
    if (status != 0) {
      fitsmess("read_fits_datafiles",status,"Error 1reading %s data",
         colnames[0]);
    } else {
      coco->Z = (int) var_int;
    }
     
    fits_read_col(fptr, TINT, colnum[1], iRow, 1, 1, 0, &var_int, &anynul,
      &status);
    if (status != 0) {
      fitsmess("read_fits_datafiles",status,"Error 2reading %s data",
         colnames[1]);
    } else {
      coco->rmJ = (int) var_int;
    }
    
    fits_read_col(fptr, TINT, colnum[2], iRow, 1, 1, 0, &var_int, &anynul,
      &status);
    if (status != 0) {
      fitsmess("read_fits_datafiles",status,"Error 3reading %s data",
         colnames[2]);
    } else {
      coco->Ncont = (int) var_int;
    }

    fits_read_col(fptr, TINT, colnum[5], iRow, 1, 1, 0, &var_int, &anynul,
      &status);
    if (status != 0) {
      fitsmess("read_fits_datafiles",status,"Error 4reading %s data",
         colnames[5]);
    } else {
      coco->Npseudo = (int) var_int;
    }

    /* allocate memory for the E_Cont, Continuum, 
     * E_Pseudo & Pseudo arrays */


    coco->E_cont = (double *) malloc_safe(sizeof(double)*coco->Ncont);
    coco->cont = (double *) malloc_safe(sizeof(double)*coco->Ncont);
    coco->E_pseudo = (double *) malloc_safe(sizeof(double)*coco->Npseudo);
    coco->pseudo = (double *) malloc_safe(sizeof(double)*coco->Npseudo);
    fits_read_col(fptr, TFLOAT, colnum[3], iRow, 1, coco->Ncont, 
                  0, &arr_float, &anynul,
      &status);
    if (status != 0) {
      fitsmess("read_fits_datafiles",status,"Error 5reading %s data, iRow=%i",
         colnames[3], iRow);
    } else {
      for (i=0;i<coco->Ncont;i++) coco->E_cont[i]=(double) arr_float[i];
    }

    fits_read_col(fptr, TFLOAT, colnum[4], iRow, 1, coco->Ncont, 
                  0, &arr_float, &anynul,
      &status);
    if (status != 0) {
      fitsmess("read_fits_datafiles",status,"Error 6reading %s data",
         colnames[4]);
    } else {
      for (i=0;i<coco->Ncont;i++) coco->cont[i]=(double) arr_float[i];
    }


    fits_read_col(fptr, TFLOAT, colnum[6], iRow, 1, coco->Npseudo, 
                  0, &arr_float, &anynul,
      &status);
    if (status != 0) {
      fitsmess("read_fits_datafiles",status,"Error 7reading %s data",
         colnames[6]);
    } else {
      for (i=0;i<coco->Npseudo;i++) coco->E_pseudo[i]=(double) arr_float[i];
    }

    fits_read_col(fptr, TFLOAT, colnum[7], iRow, 1, coco->Npseudo, 
                  0, &arr_float, &anynul,
      &status);
    if (status != 0) {
      fitsmess("read_fits_datafiles",status,"Error 8reading %s data",
         colnames[7]);
    } else {
      for (i=0;i<coco->Npseudo;i++) coco->pseudo[i]=(double) arr_float[i];
    }
    
  }
  };
