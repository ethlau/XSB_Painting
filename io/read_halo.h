#ifndef _READ_HALO_
#define _READ_HALO_

typedef struct {
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

    int id;
    double ra;
    double dec;
    double redshift_S;
    double galaxy_stellar_mass;
    double M200c, M500c, M2500c, Mvir, Mpeak; // in Msun/h
    double Acc_Rate_100Myr, Acc_Rate_1Tdyn, Acc_Rate_2Tdyn;
    double NFWconc;
    double b_to_a, b_to_a_500c, c_to_a, c_to_a_500c;
    double rs, rvir, vmax;

    //common
    double redshift;
    
    //lightcone
    int lens_id;
    int theta_x, theta_y; //pixel unitis in ray-tracing simulation
    double pec_vel; //LOS peculiar velocity in km/s
    double M200b; // in Msun/h
    
    //rockstar
    long pid;
    double x, y, z, Xoff; //in Mpc/h

} halo_struct ;

typedef struct {
    int num_halos;
    halo_struct *list;
} halo_list;

halo_list *load_halo_fits ( const char *filename );
halo_list *load_halo_run ( const char *filename );
halo_list *load_halo_rs ( const char *filename );
halo_list *load_halo_simple ( const char *filename );
void destroy_halo_list (halo_list *halos);

#endif
