#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "read_rockstar.h"

/*
int main(){
    struct binary_output_header bheader;
    struct halo *halos;
    long long *part_ids;
    char root[1024];
    long long tothalos, totids;
    long long i;

    sprintf(root,"/work/osatokn/tSZ/DM/output/snapdir_000/rockstar_000");
    load_rockstar(root, 0, &bheader, &halos, &part_ids, &tothalos, &totids);

    printf("#total halos:%lld total ids:%lld\n", tothalos, totids);
    printf("#id m num_p p_start\n");
    for(i=0;i<10;i++){
    //for(i=0;i<tothalos;i++){
        printf("%lld %g %lld %lld\n", halos[i].id, halos[i].m, halos[i].num_p, halos[i].p_start);
    }
    //printf("#pid\n");
    //for(i=0;i<100;i++) printf("%lld\n", part_ids[i]);

    return 0;
}
*/

void load_rockstar(const char *root, int snap, struct binary_output_header *bheader, struct halo **halos, long long **part_ids, float **part_pos, long long *TotHalo, long long *TotIds){
    long long nhalo, nids, tothalo, totids;
    int i;

    nhalo = 0, nids = 0, tothalo = 0, totids = 0;

    for(i=0;i<TOTAL_CHUNK;i++){
        load_header_rockstar(root, snap, i, bheader);
        tothalo += bheader->num_halos;
        totids += bheader->num_particles;
    }

    *TotHalo = tothalo, *TotIds = totids;
    printf("Total halos:%lld Total particles:%lld\n", tothalo, totids);

    *halos = (struct halo *) malloc(sizeof(struct halo)*tothalo);
    *part_ids = (long long *) malloc(sizeof(long long)*totids);
    *part_pos = (float *) malloc(sizeof(float)*totids*3);

    for(i=0;i<TOTAL_CHUNK;i++){
        load_halos_rockstar(root, snap, i, bheader, halos, part_ids, part_pos, nhalo, nids);
        nhalo += bheader->num_halos;
        nids += bheader->num_particles;
        printf("chunk:%d -> halos:%lld ids:%lld loaded\n", i, nhalo, nids);
    }

}

void load_header_rockstar(const char *root, int snap, int chunk, struct binary_output_header *bheader){
    char name[1024];
    FILE *input;

    sprintf(name, "%s/halos_%d.%d.bin", root, snap, chunk);
    if((input = fopen(name, "rb")) == NULL){
        printf("file open error:%s\n", name);
        exit(1);
    }

    fread(bheader, sizeof(struct binary_output_header), 1, input);
    fclose(input);
}

void load_halos_rockstar(const char *root, int snap, int chunk, struct binary_output_header *bheader, struct halo **halos, long long **part_ids, float **part_pos, long long nhalo, long long nids){
    char name[1024];
    long long i;
    FILE *input;

    sprintf(name, "%s/halos_%d.%d.bin", root, snap, chunk);
    if((input = fopen(name, "rb")) == NULL){
        printf("file open error:%s\n", name);
        exit(1);
    }

    fread(bheader, sizeof(struct binary_output_header), 1, input);

    fread(*halos+nhalo, sizeof(struct halo), bheader->num_halos, input);
    fread(*part_ids+nids, sizeof(long long), bheader->num_particles, input);
    fread(*part_pos+3*nids, 3*sizeof(float), bheader->num_particles, input);

    for(i=0;i<bheader->num_halos;i++) (*halos+nhalo)[i].p_start += nids;

    fclose(input);
}

void load_parent(const char *root, int snap, struct halo *halos, long long tothalos){
    char name[1024];
    long long i, id, parent_id, key;
    long long *id_map;
    FILE *input;

    sprintf(name, "%s/parent_%d.list", root, snap);
    if((input = fopen(name, "r")) == NULL){
        printf("file open error:%s\n", name);
        exit(1);
    }
    printf("parent ID loaded from %s\n", name);

    id_map = (long long *) malloc(sizeof(long long)*tothalos);

    for(i=0;i<tothalos;i++){
        id_map[halos[i].id] = i;
    }

    for(i=0;i<tothalos;i++){
        fscanf(input, "%lld %lld", &id, &parent_id);
        key = id_map[id];
        /*
        if(halos[i].id != id){
            printf("id descrepancy found!-> %lld and %lld\n", halos[i].id, id);
            exit(1);
        }
        */
        /* tentatively I use "desc" for parent ID*/
        halos[key].desc = parent_id;
    }
    free(id_map);
    fclose(input);
}

