#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "para.h"


int load_header(char *fname, struct io_header_1 *header){
    FILE *fd;
    char buf[200];
    int dummy;

    sprintf(buf,"%s.%d",fname,0);
    if(!(fd=fopen(buf,"rb"))){
        printf("can't open file `%s`\n",buf);
        exit(0);
    }
//    printf("loading header from `%s`\n", buf);

    fread(&dummy, sizeof(dummy), 1, fd);
    fread(header, sizeof(struct io_header_1), 1, fd);
    fread(&dummy, sizeof(dummy), 1, fd);

//    NumPart = (long long) header1.npartTotal[1] + (((long long)header1.npartTotal[2]) << 32);
    fclose(fd);

    return 0;
}

void load_denshsml(char *path, int snapno, int chunk, struct particle_data *P){
    FILE *fp;
    char fname[1024];
    long long nids;
    float buf;

    sprintf(fname, "%s/density_%03d.%d", path, snapno, chunk);
    if(!(fp = fopen(fname, "rb"))){
        printf("can't open file `%s`\n", fname);
        exit(0);
    }

    fread(&nids, sizeof(long long), 1, fp);
    for(long long i=0;i<nids;i++){
        fread(&buf, sizeof(float), 1, fp);
        P[i].Rho = buf;
    }
    fclose(fp);

    printf("nids:%lld\n", nids);
    sprintf(fname, "%s/hsml_%03d.%d", path, snapno, chunk);
    if(!(fp = fopen(fname, "rb"))){
        printf("can't open file `%s`\n", fname);
        exit(0);
    }

    fread(&nids, sizeof(long long), 1, fp);
    for(long long i=0;i<nids;i++){
        fread(&buf, sizeof(float), 1, fp);
        P[i].hsml = buf;
    }
    fclose(fp);

    return;
}
