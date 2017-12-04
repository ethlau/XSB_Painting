#include "readapec_spectrum.h"
void read_line_data(fitsfile *fptr, struct EMISSION_LINES *lines, int nrows){

  int iCol, iRow, anynul;
  int status=0;
  float var_float, fnul_val;
  int var_int;
  char *colnames[6] = {"Lambda", "Epsilon", "Element", "Ion", 
                       "UpperLev", "LowerLev"};
  int colnum[6]; /* 6 names listed above */
  int Ncols = 6;
  char comment[MAXSTRLEN];
  //struct EMISSION_LINES *line_dat;
  unsigned long datasum = 0;

  \

  /* Get the column numbers */
  for(iCol=0;iCol<Ncols;iCol++) {
    fits_get_colnum(fptr, CASEINSEN, colnames[iCol], &(colnum[iCol]), &status);
    if (status!=0) {
      fitsmess("read_fits_datafiles",status,"Error reading %s column number",
         colnames[iCol]);
    }
  }
  /* read in the data */
  for (iRow=1;iRow<nrows+1;iRow++) {
    fits_read_col(fptr, TFLOAT, colnum[0], iRow, 1, 1, 0, &var_float, &anynul,
      &status);
    if (status != 0) {
      fitsmess("read_fits_datafiles",status,"Error reading %s data",
         colnames[0]);
    } else {
      lines[iRow-1].lambda = (double) var_float;
    }
    
    fits_read_col(fptr, TFLOAT, colnum[1], iRow, 1, 1, 0, &var_float, &anynul,
      &status);
    if (status != 0) {
      fitsmess("read_fits_datafiles",status,"Error reading %s data",
         colnames[1]);
    } else {
      lines[iRow-1].epsilon = (double) var_float;
    }

    fits_read_col(fptr, TINT, colnum[2], iRow, 1, 1, 0, &var_int, &anynul,
      &status);
    if (status != 0) {
      fitsmess("read_fits_datafiles",status,"Error reading %s data",
         colnames[2]);
    } else {
      lines[iRow-1].Z = (double) var_int;
    }

    fits_read_col(fptr, TINT, colnum[3], iRow, 1, 1, 0, &var_int, &anynul,
      &status);
    if (status != 0) {
      fitsmess("read_fits_datafiles",status,"Error reading %s data",
         colnames[3]);
    } else {
      lines[iRow-1].rmJ = (double) var_int;
    }

    fits_read_col(fptr, TINT, colnum[4], iRow, 1, 1, 0, &var_int, &anynul,
      &status);
    if (status != 0) {
      fitsmess("read_fits_datafiles",status,"Error reading %s data",
         colnames[4]);
    } else {
      lines[iRow-1].ul = (double) var_int;
    }

    fits_read_col(fptr, TINT, colnum[5], iRow, 1, 1, 0, &var_int, &anynul,
      &status);
    if (status != 0) {
      fitsmess("read_fits_datafiles",status,"Error reading %s data",
         colnames[5]);
    } else {
      lines[iRow-1].ll = (double) var_int;
    }

  }  
};
