#include <fstream>
#include <iostream>
#include <string>
#include <stdlib.h>
#include "read_halo.h"


using namespace std;

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
    float rvir, x, y, z, M200c, M500c, Xoff, rs, redshift, Mvir;

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
    float redshift, pec_vel, rvir, rs, Mvir, M200b, M200c, M500c, M2500c;

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

void destroy_halo_list (halo_list *halos){

    int i;
    int num_halos = halos->num_halos;
    free(halos->list);
    free(halos);
}
