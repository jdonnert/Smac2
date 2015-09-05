/* MACROS & HELPER FUNCTIONS */
#define Assert(...) Assert_Info(__FILE__, __func__,  __LINE__, __VA_ARGS__)
#define Malloc(x) Malloc_info(x,  __FILE__, __func__, __LINE__)
#define Realloc(x,y) Realloc_info(x, y, __FILE__,__func__, __LINE__)
#define Free(x) Free_info(x,  __FILE__, __func__,__LINE__)

/* define a root printf to output only on the first MPI & OMP task*/
#define rprintf if (Task.Rank == 0 && Omp.ThreadID == 0) printf

#define imin(a,b) ((a)>(b)?(b):(a))
#define imax(a,b) ((a)<(b)?(b):(a))

#define length3(a) sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2])
#define length2(a) sqrt(a[0]*a[0] + a[1]*a[1])

#define length3_sq(a) (a[0]*a[0] + a[1]*a[1] + a[2]*a[2])

#define p2(a) ((a)*(a))
#define p3(a) ((a)*(a)*(a))
