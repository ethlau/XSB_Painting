#ifndef __PARAM__
#define __PARAM__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>

#define FORCE_RES (0.0244140625)
#define SQR(x) ((x)*(x))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

int write_snapshot(struct particle_data *P, struct io_header_1 header, char *path, char *basename, int snapno, int chunk);
int load_header(char *fname, struct io_header_1 *header1);
void load_denshsml(char *path, int snapno, int chunk, struct particle_data *P);

struct io_header_1
{
  int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  int      npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */
};

struct particle_data 
{
  float  Pos[3];
//  float  Vel[3];
  long long Id;
  float  Mass;
  int    Type;
  float  Rho, Pre, hsml;
  //float  Rho, U, hsml, Pot, Disp, Mean;
  //float  Ne, Temp, Pre;
};

#endif
