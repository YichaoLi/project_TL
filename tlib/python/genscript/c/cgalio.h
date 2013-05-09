#ifndef _H_CGALIO_
#define _H_CGALIO_


void read_particle_positions(char *file_inp, int npart, float *x, float *y, float *z);
void read_particle_velocity(char *file_inp, int npart, double *den);


void write_particle_positions(char *file_inp, int npart, float *x, float *y, float *z);
void write_particle_velocity(char *file_inp, int npart, double *den);


void read_particle_header(char *file_inp, int *npart, float *boxsize);

#endif
