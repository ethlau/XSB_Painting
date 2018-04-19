#ifndef _READ_HALO_
#define _READ_HALO_

typedef struct {
    //common
    long id;
    float redshift;
    float rs, rvir, Mvir; //scale radius, virial radius in comoving kpc/h, virial mass in Msun/h 
    float M200c, M500c;
    
    //lightcone
    int lens_id;
    int theta_x, theta_y; //pixel unitis in ray-tracing simulation
    float pec_vel; //LOS peculiar velocity in km/s
    float M200b, M2500c; // in Msun/h
    
    //rockstar
    long pid;
    float x, y, z, Xoff; //in Mpc/h
} halo_struct ;

typedef struct {
    int num_halos;
    halo_struct *list;
} halo_list;

halo_list *load_halo_run ( const char *filename );
halo_list *load_halo_rs ( const char *filename );
halo_list *load_halo_simple ( const char *filename );
void destroy_halo_list (halo_list *halos);

#endif
