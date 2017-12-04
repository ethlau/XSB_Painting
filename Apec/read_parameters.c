#include "readapec_spectrum.h"

void read_parameters(fitsfile *fptrline, fitsfile *fptrcoco, 
                     struct EMISSION *emission, int Nhdu) {
  int iCol, iRow, anynul;
  int status=0;
  float var_float, fnul_val;
  int var_int;
  char *colnames[7] = {"kT", "eDensity", "Time", "NCont", 
                       "NPseudo", "NElement","NLine"};
  int colnum[7]; /* 6 names listed above */
  int Ncols = 7;
  char comment[MAXSTRLEN];
  //struct EMISSION_LINES *line_dat;
  unsigned long datasum = 0;

  for(iCol=0;iCol<6;iCol++) {
    fits_get_colnum(fptrcoco, CASEINSEN, colnames[iCol], &(colnum[iCol]), &status);
    if (status!=0) {
      fitsmess("read_fits_datafiles",status,"Error reading %s column number",
         colnames[iCol]);
    }
  }

  iCol=6;
  fits_get_colnum(fptrline, CASEINSEN, colnames[iCol], &(colnum[iCol]), &status);
  if (status!=0) {
    fitsmess("read_fits_datafiles",status,"Error reading %s column number",
       colnames[iCol]);
  }


  for (iRow=1;iRow<Nhdu+1;iRow++) {
    fits_read_col(fptrcoco, TFLOAT, colnum[0], iRow, 1, 1, 0, &var_float, &anynul,
      &status);
    if (status != 0) {
      fitsmess("read_fits_datafiles",status,"Error reading %s data",
         colnames[0]);
    } else {
      emission->kT[iRow-1] = (double) var_float;
    }
  }
  
  for (iRow=1;iRow<Nhdu+1;iRow++) {
    fits_read_col(fptrcoco, TFLOAT, colnum[1], iRow, 1, 1, 0, &var_float, &anynul,
      &status);
    if (status != 0) {
      fitsmess("read_fits_datafiles",status,"Error reading %s data",
         colnames[1]);
    } else {
      emission->eDensity[iRow-1] = (double) var_float;
    }
  }
  
  for (iRow=1;iRow<Nhdu+1;iRow++) {
    fits_read_col(fptrcoco, TFLOAT, colnum[2], iRow, 1, 1, 0, &var_float, &anynul,
      &status);
    if (status != 0) {
      fitsmess("read_fits_datafiles",status,"Error reading %s data",
         colnames[2]);
    } else {
      emission->time[iRow-1] = (double) var_float;
    }
  }

  for (iRow=1;iRow<Nhdu+1;iRow++) {
    fits_read_col(fptrcoco, TINT, colnum[3], iRow, 1, 1, 0, &var_int, &anynul,
      &status);
    if (status != 0) {
      fitsmess("read_fits_datafiles",status,"Error reading %s data",
         colnames[3]);
    } else {
      emission->ncont[iRow-1] = (int) var_int;
    }
  }

  for (iRow=1;iRow<Nhdu+1;iRow++) {
    fits_read_col(fptrcoco, TINT, colnum[4], iRow, 1, 1, 0, &var_int, &anynul,
      &status);
    if (status != 0) {
      fitsmess("read_fits_datafiles",status,"Error reading %s data",
         colnames[4]);
    } else {
      emission->npseudo[iRow-1] = (int) var_int;
    }
  }
  
  for (iRow=1;iRow<Nhdu+1;iRow++) {
    fits_read_col(fptrcoco, TINT, colnum[5], iRow, 1, 1, 0, &var_int, &anynul,
      &status);
    if (status != 0) {
      fitsmess("read_fits_datafiles",status,"Error reading %s data",
         colnames[5]);
    } else {
      emission->nelement[iRow-1] = (int) var_int;
    }
  }


  for (iRow=1;iRow<Nhdu+1;iRow++) {
    fits_read_col(fptrline, TINT, colnum[6], iRow, 1, 1, 0, &var_int, &anynul,
      &status);
    if (status != 0) {
      fitsmess("read_fits_datafiles",status,"Error reading %s data",
         colnames[6]);
    } else {
      emission->nline[iRow-1] = (int) var_int;
    }
  } 

  
                         
};

