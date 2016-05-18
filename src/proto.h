/* INCLUDES */
#include <stdlib.h>		// system
#include <stdio.h>
#include <stdint.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include <float.h>

#include <mpi.h>		// Parallel Env, MPI & OpenMP
#include <omp.h>

#include <gsl/gsl_const_cgsm.h> // GNU Scientific Library
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_math.h>	

#include "config.h"		// Our compiletime options

#include "constants.h" // global header files
#include "macro.h"
#include "cosmo.h"
#include "unit.h"
#include "timing.h"

void read_param_file();
void setup();
void select_cosmology(int);
void read_snapshot(char *);
long guess_snapnum(char *);
void find_DM_densities();
void write_output();
void select_effect_module(int choice);
void Remove_Unused_Particles();
void Set_Barycenter();
extern void Find_SPH_Densities();
void domain_decomposition();
void project();

double Metric_Hsml2(int);
double Metric_One(int);


#ifdef O_FITS_COMPRESSED_CUBE
void move_image_to_cube(const int k);
#endif

/* helpers */
void print_compile_time_settings();

extern void *Malloc_info(size_t, const char *, const char *, const int);
void *Realloc_info(void *, size_t, const char *, const char *, const int);
void Assert_Info(const char *, const char *, int, int, const char *, ...);
void Free_info(void *, const char *, const char *, const int);
double *image_alloc(size_t, size_t, size_t);
void Reallocate_P(int, int *, int);
