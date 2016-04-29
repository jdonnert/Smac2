/* Remove uneeded Particles from Memory.
 * Distribute the rest evxenly over all CPUs
 * */

#include "proto.h"
#include "globals.h"
#include "effects/effects.h"
#include <gsl/gsl_sort.h>

static bool check_particle_in_image(int);
static double load_above_mean(double *);
static void xChange_particles(int, double, int, double);
static void sort_particles_pix();
static void sort_particles_metric();

double Metric_Hsml2(int);
double Metric_One(int);

void domain_decomposition()
{
	rprintf("Starting Domain Decomposition\n");

	start_timing(CPU_DOMAINDECOMP);

	/* Remove particles not in image */
	int target = 0;
	int64_t npart[N_part_types] = { 0 };

	for (int src = 0; src < Task.PartTotal; src++) {

		if (check_particle_in_image(src)) {

			memmove(&P[target], &P[src], sizeof(*P));

			if (P[src].Type == 0)
				memmove(&Gas[target], &Gas[src], sizeof(*Gas));

			target++;

			npart[P[src].Type]++;
		}
	}

	/* Count particle change */
	int dNpart[N_part_types] = { 0 }, dntot = 0;

	for (int i = 0; i < N_part_types; i++) {

		dNpart[i] = Task.Npart[i] - npart[i];

		dntot += dNpart[i];
	}

	Reallocate_P(dntot, dNpart, -1);

	int64_t nTotal = Task.PartTotal;
	for (int i = 0; i < N_part_types; i++)
		npart[i] = Task.Npart[i];

	MPI_Allreduce(&nTotal, &Snap.PartTotal, 1,
		      MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(npart, Snap.Npart, N_part_types,
		      MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

	rprintf("Particles in image:  \n"
		"   Sph   <%9lli>  DM     <%9lli>    \n"
		"   Disk  <%9lli>  Bulge  <%9lli>    \n"
		"   Star  <%9lli>  Bndry  <%9lli>    \n"
		"   Total <%9lli> \n\n",
		Snap.Npart[0], Snap.Npart[1], Snap.Npart[2],
		Snap.Npart[3], Snap.Npart[4], Snap.Npart[5], Snap.PartTotal);
	fflush(stdout);

	MPI_Barrier(MPI_COMM_WORLD);

	Assert(Snap.PartTotal > 0, "No Particles in Image");

	if (Task.NTask == 1) 
		goto exit_func;
		
	/* Distribute particles evenly */
	double *excess = (double *)Malloc(Task.NTask * sizeof *excess);
	size_t *sort = (size_t *) Malloc(Task.NTask * sizeof *sort);

	sort_particles_metric();

	int task = 0, cnt = 0, xchangeTask = 0;
	double imbalance = 0;

	for (;;) {

		if (Task.NTask == 1)
			break;
		
		double total_load = 0;
		
		double myExcess = load_above_mean(&total_load);	// returns total_load 

		MPI_Allgather(&myExcess, 1, MPI_DOUBLE,
			      excess, 1, MPI_DOUBLE, MPI_COMM_WORLD);


		for (imbalance = task = 0; task < Task.NTask; task++)
			imbalance += fabs(excess[task]) / total_load;

		if (imbalance < 1e-6 || cnt > 20)
			break;

		gsl_sort_index(sort, excess, 1, Task.NTask);

		for (int i = 0; i < Task.NTask; i++)
			if (sort[i] == Task.Rank)
				xchangeTask = sort[Task.NTask - i - 1];

		if (excess[Task.Rank] * excess[xchangeTask] < 0)
			xChange_particles(xchangeTask, excess[xchangeTask],
					  Task.Rank, excess[Task.Rank]);

		cnt++;
	}

	Free(excess);
	Free(sort);

	double total_imbalance = 0;

	MPI_Reduce(&imbalance, &total_imbalance, 1,
		   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	rprintf("\nFinal memory imbalance of <%1.4g> "
		"after <%d> iterations \n\n", total_imbalance, cnt);

	MPI_Barrier(MPI_COMM_WORLD);

	sort_particles_pix();	// bad for cache, essential for healpix
	
	exit_func: 
	
	stop_timing(CPU_DOMAINDECOMP);

	return;
}

static bool check_particle_in_image(int ipart)
{
	const double halfXYSize = 0.5 * Param.XYSize;

#ifdef PERIODIC			/* Periodic images have to be considered */
	const double boxSize = Snap.Boxsize;	// From snapshot 

	for (int k = 0; k < 8; k++) {	// Mirror Particle 
		double px = !(k & 0x1) ?
		    P[ipart].Pos[0] : (P[ipart].Pos[0] > 0 ?
				       P[ipart].Pos[0] -
				       boxSize : P[ipart].Pos[0] + boxSize);
		double py =
		    !(k & 0x2) ? P[ipart].Pos[1] : (P[ipart].Pos[1] >
						    0 ? P[ipart].Pos[1] -
						    boxSize : P[ipart].Pos[1] +
						    boxSize);
		double pz =
		    !(k & 0x4) ? P[ipart].Pos[2] : (P[ipart].Pos[2] >
						    0 ? P[ipart].Pos[2] -
						    boxSize : P[ipart].Pos[2] +
						    boxSize);
#else
	double px = P[ipart].Pos[0];
	double py = P[ipart].Pos[1];
	double pz = P[ipart].Pos[2];
#endif

#ifdef INSIDEONLY		// Thow away everything not inside img
	if (!(fabs(px) <= halfXYSize)
	    || !(fabs(py) <= halfXYSize)
	    || !(fabs(pz) <= 0.5 * Param.ZDepth)) {
#else				// Take everything touching img
	double phsml = P[ipart].Hsml * Param.XYPix / Param.XYSize;

	if (px - phsml > halfXYSize || px + phsml < -halfXYSize
	    || py - phsml > halfXYSize || py + phsml < -halfXYSize
	    || pz - phsml > halfXYSize || pz + phsml < -halfXYSize) {
#endif

#ifdef PERIODIC
		continue;
	} else {
		return true;
	}
}				// close k=0,8 - loop 

return false;
#else
		return false;
	} else {
		return true;
	}
#endif // PERIODIC
}

/* 
 * Determine excess computational 
 * load present on this CPU 
 */

double load_above_mean(double *return_ptr)
{
	double myLoad = 0;

	for (int type = 0; type < N_part_types; type++)
		for (int ipart = 0; ipart < Task.Npart[type]; ipart++)
			myLoad += (*metric_func_ptr) (ipart);

	double *load = (double *)Malloc(Task.NTask * sizeof *load);

	MPI_Allgather(&myLoad, 1, MPI_DOUBLE, load, 1, MPI_DOUBLE,
		      MPI_COMM_WORLD);

	double total = 0;

	for (int task = 0; task < Task.NTask; task++)
		total += load[task];

	double mean = total / Task.NTask;

	double excess = load[Task.Rank] - mean;

	free(load);

	*return_ptr = total;	// returns total load and excess load 

	return excess;
}

void xChange_particles(int yourRank, double yourExcess, int myRank, 
		double myExcess)
{
	int partTotal = 0, nPart[N_part_types] = { 0 };
	MPI_Status status;

	if (myExcess > 0) {	// Send particles 

		size_t ipart = Task.PartTotal - 1;

		double maxExcess = fmin(myExcess, fabs(yourExcess));

		double sendLoad = 0;

		while (sendLoad < maxExcess) {

			sendLoad += (*metric_func_ptr) (ipart);

			partTotal++;

			nPart[P[ipart].Type]++;

			ipart--;
		}

		int tag = myRank;

		MPI_Ssend(&partTotal, 1, MPI_INT, yourRank, tag,
			  MPI_COMM_WORLD);
		MPI_Ssend(nPart, N_part_types, MPI_INT, yourRank, tag,
			  MPI_COMM_WORLD);

		MPI_Ssend(&(P[Task.PartTotal - partTotal]),
			  partTotal * sizeof(struct Particle_Data), MPI_BYTE,
			  yourRank, tag, MPI_COMM_WORLD);
		MPI_Ssend(&(Gas[Task.Npart[0] - nPart[0]]),
			  nPart[0] * sizeof(struct Gas_Data), MPI_BYTE,
			  yourRank, tag, MPI_COMM_WORLD);

		Reallocate_P(partTotal, nPart, -1);

	} else {		// Receive particles 

		int tag = yourRank;

		MPI_Recv(&partTotal, 1, MPI_INT, yourRank, tag, MPI_COMM_WORLD,
			 &status);
		MPI_Recv(nPart, N_part_types, MPI_INT, yourRank, tag,
			 MPI_COMM_WORLD, &status);

		Reallocate_P(partTotal, nPart, +1);

		MPI_Recv(&(P[Task.PartTotal - partTotal]),
			 partTotal * sizeof(struct Particle_Data), MPI_BYTE,
			 yourRank, tag, MPI_COMM_WORLD, &status);
		MPI_Recv(&(Gas[Task.Npart[0] - nPart[0]]),
			 nPart[0] * sizeof(struct Gas_Data), MPI_BYTE,
			 yourRank, tag, MPI_COMM_WORLD, &status);
	}

	return;
}

static void sort_particles_metric()
{
	const int nPart = Task.PartTotal;
	const int haveGas = Task.Npart[0];

	double *cost = Malloc(nPart * sizeof(double));
	size_t *idx = Malloc(nPart * sizeof(size_t));

	for (int ipart = 0; ipart < nPart; ipart++)
		cost[ipart] = 1.0 / Metric_Hsml2(ipart);

	gsl_sort_index(idx, cost, 1, nPart);

	struct Gas_Data GasNext, GasSrc;

	for (int ipart = 0; ipart < nPart; ipart++) {

		if (idx[ipart] == ipart)
			continue;

		struct Particle_Data Psrc = P[ipart];

		if (haveGas)
			GasSrc = Gas[ipart];

		size_t dest = idx[ipart];

		for (;;) {

			struct Particle_Data Pnext = P[dest];

			if (haveGas)
				GasNext = Gas[dest];

			size_t idxNext = idx[dest];

			memcpy(&P[dest], &Psrc, sizeof(*P));

			if (haveGas)
				memcpy(&Gas[dest], &GasSrc, sizeof(*Gas));

			idx[dest] = dest;

			if (dest == ipart)
				break;

			memcpy(&Psrc, &Pnext, sizeof(*P));

			if (haveGas)
				memcpy(&GasSrc, &GasNext, sizeof(*Gas));

			dest = idxNext;
		}
	}

	free(idx);
	free(cost);

	return;
}

/* sorting local particles by pixel this might increase cache hits */
static void sort_particles_pix()
{
	const size_t nPart = Task.PartTotal;
	const int haveGas = Task.Npart[0];
	const float len2pix = Param.XYPix / Param.XYSize;	// [pix/Unit.Length] 

	double *pix = Malloc(nPart * sizeof(*pix));
	size_t *idx = Malloc(nPart * sizeof(*idx));

	for (int ipart = 0; ipart < nPart; ipart++) {

		int x = floor(P[ipart].Pos[0] * len2pix);	// find pixel number 
		int y = floor(P[ipart].Pos[1] * len2pix);

		x = imax(x, 0);
		y = imax(y, 0);
		x = imin(x, Param.XYPix - 1);
		y = imin(y, Param.XYPix - 1);

		float inc = fmin(1.0 / P[ipart].Hsml, 1.0 - FLT_MIN);	// inc < 1 !! 

		pix[ipart] = x * Param.XYPix + y + inc;
	}

	gsl_sort_index(idx, pix, 1, nPart);

	struct Gas_Data Gastmp;

	for (size_t i = 0; i < nPart; i++) {

		if (idx[i] == i)
			continue;

		size_t dest = i;

		struct Particle_Data Ptmp = P[dest];

		if (haveGas)
			Gastmp = Gas[dest];

		size_t src = idx[i];

		for (;;) {

			memcpy(&P[dest], &P[src], sizeof(*P));

			if (haveGas)
				memcpy(&Gas[dest], &Gas[src], sizeof(*Gas));

			idx[dest] = dest;

			dest = src;
			src = idx[dest];

			if (src == i)
				break;
		}

		memcpy(&P[dest], &Ptmp, sizeof(*P));

		if (haveGas)
			memcpy(&Gas[dest], &Gastmp, sizeof(*Gas));

		idx[dest] = dest;
	}

	free(pix);
	free(idx);

	return;
}

/* Return a measure of the computational
 * load carried by particle ipart
 * easy effects are dominated be 
 * pixel distribution ~h^2.
 * complex effects are dominated by particle 
 * number.
 * */
double Metric_Hsml2(int ipart)
{
	return P[ipart].Hsml * p2(Param.XYPix / Param.XYSize) * P[ipart].Hsml;
}

double Metric_One(int ipart)
{
	return 1;
}
