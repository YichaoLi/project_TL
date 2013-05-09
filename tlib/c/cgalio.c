#include <stdio.h>
#include <stdlib.h>




void read_particle_positions(char *file_inp, int npart, float *x, float *y, float *z){
  FILE   *fp;
  int   i;
  int   tempi;
  float tempf;

  fp = fopen(file_inp,"r");
  fread(&tempf,sizeof(float),1,fp);
  fread(&tempi,sizeof(int),1,fp);

  for (i=0 ; i<npart ; i++){            //--- DTFE format
    fread(&tempf,sizeof(float),1,fp);
    x[i] = tempf;
    fread(&tempf,sizeof(float),1,fp);
    y[i] = tempf;
    fread(&tempf,sizeof(float),1,fp);
    z[i] = tempf;
    }
  fclose(fp);
  }



void read_particle_velocity(char *file_inp, int npart, double *den){
  FILE   *fp;
  int    i;
  int    tempi;
  double tempd;

  fp = fopen(file_inp,"r");
  fread(&tempi,sizeof(int),1,fp);
  
  for (i=0 ; i<npart ; i++){            //--- DTFE format
    fread(&tempd,sizeof(double),1,fp);
    den[i] = tempd;
    }
  fclose(fp);
  }



/* ---- >>>  writing  <<< ----*/


void write_particle_positions(char *file_inp, int npart, float *x, float *y, float *z){
  FILE   *fp;
  int   i;
  float tempf=0;

  fp=fopen(file_inp,"wb");
  fwrite(&tempf,sizeof(float),1,fp);
  fwrite(&npart,sizeof(int),1,fp);

  for (i=0 ; i<npart; i++){            //--- DTFE format
    fwrite(&x[i],sizeof(float),1,fp);
    fwrite(&y[i],sizeof(float),1,fp);
    fwrite(&z[i],sizeof(float),1,fp);
    }
  fclose(fp);
  }



void write_particle_velocity(char *file_inp, int npart, double *den){
  FILE   *fp;
  int    i;

  fp = fopen(file_inp,"wb");
  fwrite(&npart,sizeof(int),1,fp);
  
  for (i=0; i<npart; i++){            //--- DTFE format
    fwrite(&den[i],sizeof(double),1,fp);
    }
  fclose(fp);
  }



void read_particle_header(char *file_inp, int *npart, float *boxsize){
  FILE   *fp;

  fp = fopen(file_inp,"r");
  fread(boxsize,sizeof(float),1,fp);
  fread(npart,sizeof(int),1,fp);
  fclose(fp);
  }




/* -->> Miguel's Gadget reader <<--- */


