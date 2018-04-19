#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include "para.h"

using namespace std;

int write_snapshot(struct particle_data *P, struct io_header_1 header, char *path, char *basename, int snapno, int chunk)
{
  int blklen, i, k;
  int nfile, N_tot;
  int SizeMass;
  char outname[1024];
  struct io_header_1 Header;

  Header = header;
  N_tot = Header.npart[1];

  //output name
  sprintf(outname, "%s/snapdir_%03d/%s_%03d.%d", path, snapno, basename, snapno, chunk);
  cout << "Writing snapshot:" << outname << endl;
  cout << "Number of particles:" << N_tot << endl;
  ofstream outshot(outname);

  // the number of mass array
  for(i=0,SizeMass=0;i<6;i++){
    if(Header.mass[i] > 0)
      SizeMass += 0;
    else
      SizeMass += Header.npart[i];
  }

  blklen=sizeof(Header);
  outshot.write((char *)&blklen, sizeof(int));
  outshot.write((char *)&Header, sizeof(Header));
  outshot.write((char *)&blklen, sizeof(int));


  // prepare arrays
  float * buf_f = new float[N_tot*3];

  //cout<<"Position..."<<endl;
  for(i=0;i<N_tot;i++){
      for(k=0;k<3;k++) buf_f[3*i+k] = P[i].Pos[k];
  }
  blklen=N_tot*3*sizeof(float);
  outshot.write((char *)&blklen, sizeof(int));
  outshot.write((char *)buf_f, N_tot*3*sizeof(float));
  outshot.write((char *)&blklen, sizeof(int));


  //cout<<"Velocity..."<<endl;
  for(i=0;i<N_tot;i++){
      for(k=0;k<3;k++) buf_f[3*i+k] = 0.0;
  }
  blklen=N_tot*3*sizeof(float);
  outshot.write((char *)&blklen, sizeof(int));
  outshot.write((char *)buf_f, N_tot*3*sizeof(float));
  outshot.write((char *)&blklen, sizeof(int));

  delete [] buf_f;


  //cout<<"Index..."<<endl;
  long long * index_out = new long long[N_tot];
  for(i=0;i<N_tot;i++) index_out[i] = P[i].Id;
  blklen=N_tot*sizeof(long long);
  outshot.write((char *)&blklen, sizeof(int));
  outshot.write((char *)index_out, N_tot*sizeof(long long));
  outshot.write((char *)&blklen, sizeof(int));

  delete [] index_out;


  if(SizeMass>0){
  float * buf_mass = new float[N_tot];

  //write mass block for split runs
  //cout<<"Mass..."<<endl;

  blklen=N_tot*sizeof(float);

  outshot.write((char *)&blklen,     sizeof(int));
  outshot.write((char *)buf_mass,    N_tot*sizeof(float));
  outshot.write((char *)&blklen,     sizeof(int));

  delete [] buf_mass;
  }

  float * buf1 = new float[N_tot];

  //cout<<"Density..."<<endl;

  for(i=0;i<N_tot;i++) buf1[i] = P[i].Rho;
  blklen = sizeof(float)*N_tot;
  outshot.write((char *)&blklen, sizeof(int));
  outshot.write((char *)buf1,    N_tot*sizeof(float));
  outshot.write((char *)&blklen, sizeof(int));


  //cout<<"Pressure..."<<endl;
  for(i=0;i<N_tot;i++) buf1[i] = P[i].Pre;
  blklen = sizeof(float)*N_tot;
  outshot.write((char *)&blklen, sizeof(int));
  outshot.write((char *)buf1,    N_tot*sizeof(float));
  outshot.write((char *)&blklen, sizeof(int));

  //cout<<"Smoothing length..."<<endl;
  for(i=0;i<N_tot;i++) buf1[i] = P[i].hsml;
  blklen = sizeof(float)*N_tot;
  outshot.write((char *)&blklen, sizeof(int));
  outshot.write((char *)buf1,    N_tot*sizeof(float));
  outshot.write((char *)&blklen, sizeof(int));


  delete [] buf1;

  outshot.close();

  return 0;
}
