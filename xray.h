#ifndef _XRAY_
#define _XRAY_

/* Array lengths for lambda tables */
#define nemax           100
#define emin		0.5
#define emax		2.0
#define ntmax		300
#define tmin		0.01
#define tmax		30.0
#define nzmax		10
#define zlmin		-6.00
#define zlmax		2.0
#define rlmin           1.e-6
#define rlmax           2.2
#define nrmax           50

extern double tarray[ntmax]; //keV
extern double zarray[nzmax]; //Solar unit
extern double rarray[nrmax]; //Solar unit
extern double lambda_table[ntmax][nrmax];
extern double tres, zres, rres, eres;

#ifdef __cplusplus
/*
extern "C" void set_lambda_table ( double zobs, double tarray[ntmax], double zarray[nzmax], double lambda_table[ntmax][nzmax] );
extern "C" void set_lambdaeff_table ( double zobs, char *filename,  double tarray[ntmax], double zarray[nzmax], double lambdaeff_table[ntmax][nzmax] );
extern "C" int read_effarea ( char *filename, double **energy, double **effarea );
extern "C" double int_lambda_table (double temp, double metal, double tarray[ntmax], double zarray[nzmax], double lambda_table[ntmax][nzmax] );
I*/
extern "C" void set_lambda_table ( double tarray[ntmax], double rarray[nrmax], double lambda_table[ntmax][nrmax] );
extern "C" double int_lambda_table (double temp, double redshift, double tarray[ntmax], double rarray[nrmax], double lambda_table[ntmax][nrmax] );


#endif

#endif
