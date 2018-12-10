#include <fstream>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <math.h>
#include <fitsio.h>
#include "read_halo.h"


using namespace std;

void printerror( int status );

halo_list *load_halo_simple ( const char *filename ) {

    halo_list *halos;
    return halos;
}

halo_list *load_halo_rs ( const char *filename ) {
    int i;
    int halo_count = 0;
    int comment_lines = 0;
    string line;

    long id, pid;
    double rvir, x, y, z, M200c, M500c, Xoff, rs, redshift, Mvir;

    halo_list *halos;

    cout << "Reading from File: " << filename << endl;
    ifstream file(filename);

    if (file.is_open()) {
        while (getline(file, line)){
	    if (line[0]!='#'){
		++halo_count;
	    }
	    else {
		++comment_lines;
	    }
	}
        cout << "Number of halos: " << halo_count << endl;
        halos = (halo_list *)malloc(sizeof(halo_list));
        halos->num_halos = halo_count;

        halos->list = (halo_struct *) malloc (halo_count*sizeof(halo_struct));

        file.clear();
        file.seekg(0, ios::beg);
	for (i = 0; i < comment_lines; i++ ) getline(file, line); 
        for (i = 0; i < halos->num_halos; i++) {
            file >> id >> pid >> redshift >> x >> y >> z >> rvir >> rs >> Mvir >> M200c >> M500c >> Xoff;
            halos->list[i].id = id;
            halos->list[i].pid = pid;
            halos->list[i].x = x;
            halos->list[i].y = y;
            halos->list[i].z = z;
	    halos->list[i].redshift = redshift;
            halos->list[i].rvir = rvir;
            halos->list[i].rs = rs;
            halos->list[i].Mvir = Mvir;
	    halos->list[i].M500c = M500c;
            halos->list[i].M200c = M200c;
	    halos->list[i].Xoff = Xoff;
        }

        file.close();

    } else {
        cout << "Unable to open file : " <<  filename << endl;
        exit(1);
    }

    return halos;
}

halo_list *load_halo_run ( const char *filename ) {
    int i;
    int halo_count = 0;
    string line;

    int lens_id, theta_x, theta_y, debug;
    double redshift, pec_vel, rvir, rs, Mvir, M200b, M200c, M500c, M2500c;

    halo_list *halos;

    ifstream file(filename);

    if (file.is_open()) {
        while (getline(file, line))
            ++halo_count;

        cout << "Number of halos: " << halo_count << endl;
        halos = (halo_list *)malloc(sizeof(halo_list));
        halos->num_halos = halo_count;

        halos->list = (halo_struct *) malloc (halo_count*sizeof(halo_struct));

        file.clear();
        file.seekg(0, ios::beg);

        for (i = 0; i < halos->num_halos; i++) {
            file >> lens_id >> theta_x >> theta_y >> redshift >> pec_vel >> rvir >> rs >> Mvir >> M200b >> M200c >> M500c >> M2500c >> debug;
            halos->list[i].id = i;
            halos->list[i].lens_id = lens_id;
            halos->list[i].theta_x = theta_x;
            halos->list[i].theta_y = theta_y;
            halos->list[i].redshift = redshift;
            halos->list[i].pec_vel = pec_vel;
            halos->list[i].rvir = rvir;
            halos->list[i].rs = rs;
            halos->list[i].Mvir = Mvir;
            halos->list[i].M200b = M200b;
            halos->list[i].M200c = M200c;
            halos->list[i].M500c = M500c;
            halos->list[i].M2500c = M2500c;
            
        }

        file.close();

    } else {
        cout << "Unable to open file" << endl;
        exit(1);
    }
    return halos;
}

halo_list *load_halo_fits ( const char *filename ) {

    int status, hdunum, hdutype, anynull;
    double doublenull;
    long nrows;

    int i, ii;
    int halo_count = 0;
    string line;

    int lens_id, theta_x, theta_y, debug;
    double ra, dec, redshift_R, redshift_S, galaxy_stellar_mass, 
           M200c, M500c, M2500c, Mvir, Mpeak,
           Acc_Rate_100Myr, Acc_Rate_1Tdyn, Acc_Rate_2Tdyn, 
           NFWconc,
           b_to_a, b_to_a_500c,
           c_to_a, c_to_a_500c,
           rs, rvir, vmax;

    halo_list *halos;

    fitsfile *fptr;

    status = 0; 
    if ( fits_open_file(&fptr, filename, READONLY, &status) )
         printerror( status );
    /* move to the HDU of the binary table */

    hdunum = 2;
    if ( fits_movabs_hdu(fptr, hdunum, &hdutype, &status) )
        printerror( status );

    if (hdutype == BINARY_TBL) {
        printf("\nReading binary table in HDU %d:\n", hdunum);
    } else {
        printf("Error: this HDU is not a binary table\n");
        printerror( status );
    }

    if (fits_get_num_rows (fptr, &nrows, &status)) 
        printerror( status );

    halo_count = (int)nrows;

    cout << "Number of halos: " << halo_count << endl;
    halos = (halo_list *)malloc(sizeof(halo_list));
    halos->num_halos = halo_count;
    halos->list = (halo_struct *) malloc (halo_count*sizeof(halo_struct));

    for (i = 0; i < halos->num_halos; i++) {
      /*  read the columns */
      /*
      Column Name                Format     Dims       Units     TLMIN  TLMAX
      1 ra                         D                   degree
      2 dec                        D                   degree
      3 redshift_R                 D                   real space
      4 redshift_S                 D                   redshift space
      5 galaxy_stellar_mass        D                   M_sun
      6 M200c                      D                   log10(M/[M_sun])
      7 M500c                      D                   log10(M/[M_sun])
      8 M2500c                     D                   log10(M/[M_sun])
      9 Mvir                       D                   log10(M/[M_sun])
     10 Mpeak                      D                   log10(M/[M_sun])
     11 Acc_Rate_100Myr            D                   Msun/yr
     12 Acc_Rate_1Tdyn             D                   Msun/yr
     13 Acc_Rate_2Tdyn             D                   Msun/yr
     14 NFWconcentration           D                   -
     15 b_to_a                     D                   -
     16 b_to_a_500c                D                   -
     17 c_to_a                     D                   -
     18 c_to_a_500c                D                   -
     19 rs                         D                   kpc/h
     20 rvir                       D                   kpc/h
     21 vmax                       D                   km/s
*/

        if (fits_read_col(fptr, TDOUBLE, 1, i+1, 1, 1, 0, &ra, &anynull, &status)) printerror( status );
        if (fits_read_col(fptr, TDOUBLE, 2, i+1, 1, 1, 0, &dec, &anynull, &status)) printerror( status );
        if (fits_read_col(fptr, TDOUBLE, 3, i+1, 1, 1, 0, &redshift_R, &anynull, &status)) printerror( status );
        if (fits_read_col(fptr, TDOUBLE, 4, i+1, 1, 1, 0, &redshift_S, &anynull, &status)) printerror( status );
        if (fits_read_col(fptr, TDOUBLE, 5, i+1, 1, 1, 0, &galaxy_stellar_mass, &anynull, &status)) printerror( status );
        if (fits_read_col(fptr, TDOUBLE, 6, i+1, 1, 1, 0, &M200c, &anynull, &status)) printerror( status );
        if (fits_read_col(fptr, TDOUBLE, 7, i+1, 1, 1, 0, &M500c, &anynull, &status)) printerror( status );
        if (fits_read_col(fptr, TDOUBLE, 8, i+1, 1, 1, 0, &M2500c, &anynull, &status)) printerror( status );
        if (fits_read_col(fptr, TDOUBLE, 9, i+1, 1, 1, 0, &Mvir, &anynull, &status)) printerror( status );
        if (fits_read_col(fptr, TDOUBLE, 10, i+1, 1, 1, 0, &Mpeak, &anynull, &status)) printerror( status );
        if (fits_read_col(fptr, TDOUBLE, 11, i+1, 1, 1, 0, &Acc_Rate_100Myr, &anynull, &status)) printerror( status );
        if (fits_read_col(fptr, TDOUBLE, 12, i+1, 1, 1, 0, &Acc_Rate_1Tdyn, &anynull, &status)) printerror( status );
        if (fits_read_col(fptr, TDOUBLE, 13, i+1, 1, 1, 0, &Acc_Rate_2Tdyn, &anynull, &status)) printerror( status );
        if (fits_read_col(fptr, TDOUBLE, 14, i+1, 1, 1, 0, &NFWconc, &anynull, &status)) printerror( status );
        if (fits_read_col(fptr, TDOUBLE, 15, i+1, 1, 1, 0, &b_to_a, &anynull, &status)) printerror( status );
        if (fits_read_col(fptr, TDOUBLE, 16, i+1, 1, 1, 0, &b_to_a_500c, &anynull, &status)) printerror( status );
        if (fits_read_col(fptr, TDOUBLE, 17, i+1, 1, 1, 0, &c_to_a, &anynull, &status)) printerror( status );
        if (fits_read_col(fptr, TDOUBLE, 18, i+1, 1, 1, 0, &c_to_a_500c, &anynull, &status)) printerror( status );
        if (fits_read_col(fptr, TDOUBLE, 19, i+1, 1, 1, 0, &rs, &anynull, &status)) printerror( status );
        if (fits_read_col(fptr, TDOUBLE, 20, i+1, 1, 1, 0, &rvir, &anynull, &status)) printerror( status );
        if (fits_read_col(fptr, TDOUBLE, 21, i+1, 1, 1, 0, &vmax, &anynull, &status)) printerror( status );

        halos->list[i].id = i;
        halos->list[i].ra = ra;
        halos->list[i].dec = dec;
        halos->list[i].redshift = redshift_R;
        halos->list[i].redshift_S = redshift_S;
        halos->list[i].galaxy_stellar_mass = galaxy_stellar_mass;
        halos->list[i].Mvir = pow(10., Mvir);
        halos->list[i].M200c = pow(10.,M200c);
        halos->list[i].M500c = pow(10.,M500c);
        halos->list[i].M2500c = pow(10.,M2500c);
        halos->list[i].Mpeak = pow(10.,Mpeak);
        halos->list[i].Acc_Rate_100Myr = Acc_Rate_100Myr;
        halos->list[i].Acc_Rate_1Tdyn = Acc_Rate_1Tdyn;
        halos->list[i].Acc_Rate_2Tdyn = Acc_Rate_2Tdyn; 
        halos->list[i].NFWconc = NFWconc;
        halos->list[i].b_to_a = b_to_a;
        halos->list[i].b_to_a_500c = b_to_a_500c;
        halos->list[i].c_to_a = c_to_a;
        halos->list[i].c_to_a_500c = c_to_a_500c;
        halos->list[i].rs = rs;
        halos->list[i].rvir = rvir;
        halos->list[i].vmax = vmax;
 
    }

    if ( fits_close_file(fptr, &status) ) 
         printerror( status );

    return halos;
}

void destroy_halo_list (halo_list *halos){

    int i;
    int num_halos = halos->num_halos;
    free(halos->list);
    free(halos);
}

void printerror( int status) {
    /*****************************************************/
    /* Print out cfitsio error messages and exit program */
    /*****************************************************/
    if (status) {
       fits_report_error(stderr, status); /* print error report */

       exit( status );    /* terminate the program, returning error status */
    }
    return;
}
