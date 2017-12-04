#ifndef _READ_LIGHTCONE_
#define _READ_LIGHTCONE_

typedef struct {
    int id;
    int lens_id;
    int theta_x, theta_y; //pixel unitis in ray-tracing simulation
    float redshift;
    float pec_vel; //LOS peculiar velocity in km/s
    float rs, rvir, Mvir; //scale radius, virial radius in comoving kpc/h, virial mass in Msun/h
    float M200b, M200c, M500c, M2500c; // in Msun/h
} halo_struct ;

typedef struct {
    int num_halos;
    halo_struct *list;
} halo_list;

halo_list *load_halo_run ( const char *filename );
void destroy_halo_list (halo_list *halos);

#endif
