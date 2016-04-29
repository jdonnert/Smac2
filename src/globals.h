/* CODE PARAMETERS */
#define N_part_types 6		// Number of particle species

#define MAXLINELENGTH 500	// of any char buffer !
#define MAXTAGS 300		// In parameter file
#define MAXIMAGES 3		// Maximum Number of images processed at once

#define SAMPLING_SCALE 1	// [pix] where SPH sampl. becomes pixel sampl.

#define VERSION "1.4"

/* 
 * SPH kernel 
 */

#ifdef WK_CUBIC
#define DESNNGB 64			// desired number of neighbours
#define NGBDEV 1			// ngb tolerance
#define NGBMAX (3*DESNNGB)	// SPH : length of neighbour list
#endif

#ifdef WK_QUINTIC
#define DESNNGB 192			
#define NGBDEV 1			
#define NGBMAX (3*DESNNGB)	
#endif

#ifdef WK_WC6
#define DESNNGB 295			
#define NGBDEV 1			
#define NGBMAX (3*DESNNGB)	
#endif

#ifdef WK_WC4
#define DESNNGB 200
#define NGBDEV 1			
#define NGBMAX (3*DESNNGB)	
#endif

/* 
 * MAIN VARIABLES 
 */

double *Image, *Weight_Image;

extern struct ParallelInfos {
	int Rank;		// Rank of local Processor 
	int NTask;		// Number of Processors 
	int PartTotal;		// Local Part. statistics 
	int Npart[N_part_types];
} Task;

extern struct OpenMP_infos {
	int NThreads;		// Number of openMP threads 
	int ThreadID;		// Thread ID of this thread 
} Omp;
#pragma omp threadprivate(Omp)

extern struct SmacProp {	/*parameter from par file */
	int Cosmology;
	double Center[3];
	int Flag_Barycenter;	/* Use Barycenter of Particles */
	double XYSize;		/* Image Size in Unit.Length */
	double ZDepth;		/* Projection maximum depth */
	int XYPix;		/* Image Resolution */
	double EulAng[3];	/* 3 Euler angles for projection */
	int NSide;		/* HEALPIX resolution parameter */
	double Rmin;		/* HEALPIX minimum distance */
	double Rmax;		/* HEALPIX maximum distance */
	int N_IOTasks;		/* Number of read tasks */
	char Input_File[MAXLINELENGTH];
	char Output_File[MAXLINELENGTH];
	int NoClobber;		/* Overwrite output file */
	int Effect_Module;	/* Effect selection */
	int Effect_Flag;	/* Effects sub flag */
	double E_min;		/* Energy Range Min */
	double E_max;		/* Energy Range Max */
	double Freq;		/* Frequency of observation */
	double a_crp;		/* CR spectral index */
	double X_crp;		/* CR energy density rel. to thermal */
	double TurbTime;	/* Timescale of turb.injection */
	double Eta_t;		/* Fraction of Energy in Magnetosonic waves */
	double TurbScale;	/* Scale of turbulent velocity */
	int SynchroPAngInt;	/* Toggle pitch angle integration */
	int SynchroIntrRM;	/* Toggle Intrinsic RM */
	int CubeType;
	int NCube;
	double CubeMin;
	double CubeMax;
	double CubeStep;
	void *CubePtr;
} Param;

extern struct Snap_Properties {
	int Have_Arepo;
	double Boxsize;		/* From Snap Header */
	long SnapNum;		/* Guessed snapshot number */
	double Redshift;
	double Time;
	uint64_t Npart[N_part_types];	/* over all processes */
	uint64_t PartTotal;	/* over all processes */
	double Masstab[N_part_types];
} Snap;

extern struct EffectProp {
	char Name[MAXLINELENGTH];
	char Descr[MAXIMAGES + 3][MAXLINELENGTH];
	char Unit[MAXLINELENGTH];
	struct requirements {	// of effect 
		bool PartType[N_part_types];	// Particle Species 
		char Block[1000][5];	// Blocks required, 4 letters  
		bool Tree;
	} Req;
	int Nimage;		// Number images projected 
} Effect;

extern struct Particle_Data {
	float Pos[3];		// Position vector 
	float Vel[3];		// Velocity vector 
	float Rho;			// Density 
	int Type;
	float Mass;
	float Hsml;			// Smoothing Length 
	float VarHsmlFac;
	unsigned long ID;	// Unique ID 
} *P;

extern struct Gas_Data {
	float U;		/* Internal Energy */
	float Bfld[3];		/* Magnetic Field */
#ifdef DIVB
	float DivB;		/* Magnetic Field Divergence */
#endif
#ifdef I_PIERNIK
	float ECR;		/* CR energy density */
#endif
#ifdef VTURB
	float VTurb;		/* Turbulent Velocity */
	float VRms;			/* RMS around mean  */
	float VTan;			/* RMS around mean tangential comp. */
	float VRad;			/* RMS around mean radial comp. */
	float VBulk[3];		/* Mean Velocity */
	float DivVel;		/* Velocity Divergence */
	float CurlVel;		/* Velocity Curl */
	float TNgb;			/* True Number of Neighbours */
	float Dpp;			/* Magnetosonic Reacceleration coeff. */
#endif
#ifdef RADTRANSFER
	float nHII;
	float nHeII;
	float nHeIII;
#endif
#ifdef SHOCKS
	float NMach;		/* Mach Number */
#endif
#ifdef VORTICITY
	float Vorticity[3];
#endif
#ifdef TIME_DEP_ART_VISC
	float AlphaVisc;
#endif
#ifdef VSPH
	float VelSPH[3];	// kernel weighted velocity
	float VelBulk[3];	// mean of kernel weighted velocity in hsml
	float VelTurb;		// dispersion of kernel weighted velocities in hsml
#endif
} *Gas;

#ifdef HEALPIX
struct Healpix_properties {
	size_t Npix;
	double Area;
	double Angle;
	double ZRange;
	double *Coord[3];
} Healpix;
#endif

struct tree_node {
	int down;		// To first daughter node or target part
	int next;		// To next subnode of father or to unkle 
	float pos[3];	// Center of node 
	int npart;		// Number of particles in node 
	float size;		// Spatial extent of node 
} *tree;

/* Parameter File Tags, also used to write fits header */
int id[MAXTAGS];
void *addr[MAXTAGS];
char tag[MAXTAGS][50], comment[MAXTAGS][50];

struct units Comv2phys;

/* function pointers */
extern void (*prep_func_ptr) ();
extern void ((*effect_func_ptr) (int, double *));
extern void (*post_func_ptr) ();
extern double (*weight_func_ptr) (int);
extern double (*metric_func_ptr) (int);
