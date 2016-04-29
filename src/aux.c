/*Here we initialise the common variables*/
#include "proto.h"	
#include "globals.h"

/* Global Variables */

struct SmacProp Param;
struct Snap_Properties Snap;
struct EffectProp Effect;
struct Particle_Data *P = NULL;
struct Gas_Data *Gas = NULL;

void ((*effect_func_ptr) (int, double *)) = NULL;
double (*weight_func_ptr) (int) = NULL;
double (*metric_func_ptr) (int) = NULL;
void (*pre_func_ptr) (size_t, size_t, double) = NULL;
void (*post_func_ptr) () = NULL;
void (*prep_func_ptr) () = NULL;

struct ParallelInfos Task;

struct OpenMP_infos Omp = { 0, 0 };

/* 
 * Memory Management
 */

void *Malloc_info(size_t size, const char *file, const char *func,
		  const int line)
{
	void *result = malloc(size);

	Assert_Info(file, func, line, result != NULL || size == 0,
		    "Allocation failed, %zu Bytes \n", size);

	return result;
}

void *Realloc_info(void *ptr, size_t size, const char *file, const char *func,
		   const int line)
{
	if (size == 0) {
		
		Free(ptr);

		return NULL;
	}

	void *result = realloc(ptr, size);

	Assert_Info(file, func, line, result != NULL || size == 0,
		    "Reallocation failed: %zu bytes \n", size);

	return result;
}

void Free_info(void *ptr, const char *file, const char *func, const int line)
{
	if (ptr != NULL)
		free(ptr);
	else
		fprintf(stderr,
			"\nWARNING ! Task %d tried to free a NULL pointer at "
			"file %s, \n     function %s() : line %d \n", Task.Rank,
			file, func, line);

	return;
}

double *image_alloc(size_t xSize, size_t ySize, size_t nImage)
{
	size_t nBytes = nImage * xSize * ySize * sizeof(double);
	double *result = Malloc(nBytes);

	memset(result, 0, nBytes);

	return result;
}

/* 
 * Reallocates the Particle structures P and Gas. Takes the relative change
 * as argument, not the total number. Add or Remove via sign argument. 
 */

extern void Reallocate_P(int partTotal, int nPart[N_part_types], int sign)
{
	if (partTotal == 0)
		return;

	Task.PartTotal += sign * partTotal;

	for (int type = 0; type < N_part_types; type++) {

		Task.Npart[type] += sign * nPart[type];

		Assert(Task.Npart[type] >= 0, "Can't allocate negative particles");
	}
	
	P = Realloc((void *)P, sizeof(*P) * Task.PartTotal);

	if (nPart[0] != 0) 
		Gas = Realloc((void *)Gas, sizeof(*Gas) * Task.Npart[0]);
	else
		Gas = NULL;

	return;
}

/* Error Handling, we use variable arguments to be able
 * to print more informative error messages 
 */
void Assert_Info(const char *file, const char *func, int line, int expr,
	    const char *errmsg, ...)
{
	if (expr)
		return;

	va_list varArgList;

	va_start(varArgList, errmsg);

	/* we fucked up, tell them */

	fprintf(stderr, "\nERROR ! Task %d, %s: %s(), %d : \n\n	",
		Task.Rank, file, func, line);

	vfprintf(stderr, errmsg, varArgList);

	fprintf(stderr, "\n\n");

	fflush(stderr);

	va_end(varArgList);

	MPI_Abort(MPI_COMM_WORLD, -1);	// finish him ...

	exit(-1);		// ... fatality

	return;
}
