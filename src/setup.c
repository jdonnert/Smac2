/* Center the image and apply the projection */

#include "proto.h"
#include "globals.h"

#ifdef HEALPIX
#include "healpix.h"
#include <chealpix.h>
#endif

static void set_barycenter();

#ifndef HEALPIX			// Normal image
static void apply_projection();
#endif

/* 
 * We allocate the image, remove comoving, 
 * put the selected Center and apply the 
 * projection. 
 */

void setup()
{
	start_timing(CPU_SETUP);

	if (Param.XYSize * Param.ZDepth == 0) {

		rprintf("Image size/depth = 0, projecting whole box \n\n");

		Param.XYSize = Param.ZDepth = Snap.Boxsize;

		Param.Center[0] = Param.Center[1] = Param.Center[2] =
		    Snap.Boxsize / 2;
	}

	Param.Center[0] *= Comv2phys.Length;
	Param.Center[1] *= Comv2phys.Length;
	Param.Center[2] *= Comv2phys.Length;

	if (Param.Flag_Barycenter)
		set_barycenter();

	#pragma omp parallel for
	for (int ipart = 0; ipart < Task.PartTotal; ipart++) {

		P[ipart].Pos[0] -= Param.Center[0];
		P[ipart].Pos[1] -= Param.Center[1];
		P[ipart].Pos[2] -= Param.Center[2];
	}

#ifndef HEALPIX			// direct projection
	apply_projection();

	Image = image_alloc(Param.XYPix, Param.XYPix, Effect.Nimage);
	Weight_Image = image_alloc(Param.XYPix, Param.XYPix, 1);

	Param.XYSize *= Comv2phys.Length;	// Make Physical 
	Param.ZDepth *= Comv2phys.Length;

	rprintf("Image Centered to <%6.3f,%6.3f,%6.3f> kpc physical \n"
		"Image Projection Psi = %g°, Theta = %g°, Phi = %g° \n"
		"Image Depth set to <%6.0f> kpc\n"
		"Image Size  set to <%6.0f> kpc x <%6.0f> kpc\n"
		"                or <%6i> pix x <%6i> pix\n\n",
		Param.Center[0] * Unit.Length / kpc2cgs,
		Param.Center[1] * Unit.Length / kpc2cgs,
		Param.Center[2] * Unit.Length / kpc2cgs,
		Param.EulAng[0], Param.EulAng[1], Param.EulAng[2],
		Param.ZDepth * Unit.Length / kpc2cgs,
		Param.XYSize * Unit.Length / kpc2cgs,
		Param.XYSize * Unit.Length / kpc2cgs, Param.XYPix, Param.XYPix);
	fflush(stdout);

#else				// HEALPIX Projection

	Healp.Npix = 12 * pow(Param.NSide, 2);
	Healp.Area = 4 * pi / Healp.Npix;	// [rad] 
	Healp.Angle = sqrt(Healp.Area);
	Healp.ZRange = Param.Rmax - Param.Rmin;

	/* Angular Coordinates of every Pixel */
	Healp.Coord[0] = image_alloc(Healp.Npix, 1, 1);
	Healp.Coord[1] = image_alloc(Healp.Npix, 1, 1);
	Healp.Coord[2] = image_alloc(Healp.Npix, 1, 1);

	for (i = 0; i < Healp.Npix; i++)
		pix2vec_nest(Param.NSide, i, &Healp.Coord[0][i]);

	Image = Image_alloc(Healp.Npix, 1, Effect.Nimage);
	Weight_Image = image_alloc(Healp.Npix, 1, 1);

	rprintf("Position centered to <%6.3f,%6.3f,%6.3f> kpc physical \n\n"
		"Using HEALPIX Projection: Nside = <%d>, Npix = <%d> \n"
		"Data Range : <%g> - <%g> = <%g> kpc \n\n",
		Param.Center[0] * Unit.Length / kpc2cgs,
		Param.Center[1] * Unit.Length / kpc2cgs,
		Param.Center[2] * Unit.Length / kpc2cgs,
		Param.NSide, Healp.Npix,
		Param.Rmax * Unit.Length / kpc2cgs,
		Param.Rmin * Unit.Length / kpc2cgs,
		Healp.ZRange * Unit.Length / kpc2cgs);

#endif				// HEALPIX

	stop_timing(CPU_SETUP);

	return;
}

/* all vector quantities have to be rotated ! */
void apply_projection()
{
	const double psi = Param.EulAng[0] * deg2rad;
	const double theta = Param.EulAng[1] * deg2rad;
	const double phi = Param.EulAng[2] * deg2rad;

	if (!phi && !theta && !psi)
		return;		// nothing to do 

	const int partTotal = Task.PartTotal;
	const int *npart = Task.Npart;

	/* Define rotation matrix */

#ifdef YAWPITCHROLL	// Luftfahrtnorm (DIN 9300) (Yaw-Pitch-Roll, Z, Y’, X’’)
	
	const double A[3][3] = {
		{cos(theta) * cos(phi), cos(theta) * sin(phi), -sin(theta)},
		{sin(psi) * sin(theta) * cos(phi) - cos(psi) * sin(phi),
		 sin(psi) * sin(theta) * sin(phi) + cos(psi) * cos(phi),
		 sin(psi) * cos(theta)},
		{cos(psi) * sin(theta) * cos(phi) + sin(psi) * sin(phi),
		 cos(psi) * sin(theta) * sin(phi) - sin(psi) * cos(phi),
		 cos(psi) * cos(phi)}
	};

#else				// Euler Matrix, y-Convention
	
	const double A[3][3] = {
		{-sin(psi) * sin(phi) + cos(psi) * cos(theta) * cos(phi),
		 -sin(psi) * cos(phi) - cos(psi) * cos(theta) * sin(phi),
		 cos(psi) * sin(theta)},
		{cos(psi) * sin(phi) + sin(psi) * cos(theta) * cos(phi),
		 cos(psi) * cos(phi) - sin(psi) * cos(theta) * sin(phi),
		 sin(psi) * sin(theta)},
		{-sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)}
	};

#endif				// YAWPITCHROLL

	/* P */

	float x[3] = { 0 };
	const size_t nBytes = 3 * sizeof(*x);

	#pragma omp parallel for
	for (int ipart = 0; ipart < partTotal; ipart++) {

		memcpy(&x, P[ipart].Pos, nBytes);
		
		P[ipart].Pos[0] =
		    A[0][0] * x[0] + A[0][1] * x[1] + A[0][2] * x[2];
		P[ipart].Pos[1] =
		    A[1][0] * x[0] + A[1][1] * x[1] + A[1][2] * x[2];
		P[ipart].Pos[2] =
		    A[2][0] * x[0] + A[2][1] * x[1] + A[2][2] * x[2];

		memcpy(&x, P[ipart].Vel, nBytes);
		
		P[ipart].Vel[0] =
		    A[0][0] * x[0] + A[0][1] * x[1] + A[0][2] * x[2];
		P[ipart].Vel[1] =
		    A[1][0] * x[0] + A[1][1] * x[1] + A[1][2] * x[2];
		P[ipart].Vel[2] =
		    A[2][0] * x[0] + A[2][1] * x[1] + A[2][2] * x[2];
	}

	/* Gas */
	#pragma omp parallel for
	for (int ipart = 0; ipart < npart[0]; ipart++) {

		memcpy(&x, Gas[ipart].Bfld, nBytes);
		
		Gas[ipart].Bfld[0] =
		    A[0][0] * x[0] + A[0][1] * x[1] + A[0][2] * x[2];
		Gas[ipart].Bfld[1] =
		    A[1][0] * x[0] + A[1][1] * x[1] + A[1][2] * x[2];
		Gas[ipart].Bfld[2] =
		    A[2][0] * x[0] + A[2][1] * x[1] + A[2][2] * x[2];

#ifdef VTURB
		memcpy(&x, Gas[ipart].VBulk, nBytes);

		Gas[ipart].VBulk[0] =
		    A[0][0] * x[0] + A[0][1] * x[1] + A[0][2] * x[2];
		Gas[ipart].VBulk[1] =
		    A[1][0] * x[0] + A[1][1] * x[1] + A[1][2] * x[2];
		Gas[ipart].VBulk[2] =
		    A[2][0] * x[0] + A[2][1] * x[1] + A[2][2] * x[2];
#endif
	}

	return;
}

/*
 * Find Center of Mass of the snapshot 
 */

static double mtot = 0, buf[3] = { 0 };

void set_barycenter()
{
 	mtot = buf[0] = buf[1] = buf[2] = 0;
	
	for (size_t ipart = 0; ipart < Task.PartTotal; ipart++) {

		buf[0] += P[ipart].Mass * P[ipart].Pos[0];
		buf[1] += P[ipart].Mass * P[ipart].Pos[1];
		buf[2] += P[ipart].Mass * P[ipart].Pos[2];

		mtot += P[ipart].Mass;
	}

	MPI_Allreduce(MPI_IN_PLACE, buf, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	MPI_Allreduce(MPI_IN_PLACE, &mtot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	buf[0] /= mtot;
	buf[1] /= mtot;
	buf[2] /= mtot;

	if (Param.Flag_Barycenter == 2) {

		rprintf("Setting Center relative to Barycenter \n");

		Param.Center[0] += buf[0];
		Param.Center[1] += buf[1];
		Param.Center[2] += buf[2];

	} else {

		rprintf("Setting center to Barycenter \n");

		Param.Center[0] = buf[0];
		Param.Center[1] = buf[1];
		Param.Center[2] = buf[2];
	}

	return;
}



/* This removes unused particle _types_ only. */
void Remove_Unused_Particles()
{
	int nLeft = Task.PartTotal;

	int src = Task.Npart[0], dest = 0;

	rprintf("Removing superfluous particle types : ");

	for (int type = 0; type < N_part_types-1; type++) {

		nLeft -= Task.Npart[type];

		if (Effect.Req.PartType[type] || !Task.Npart[type]) {	// next

			src += Task.Npart[type + 1];

			dest += Task.Npart[type];

			continue;
		}

		rprintf("%d ", type);

		Task.PartTotal -= Task.Npart[type];
		Task.Npart[type] = 0;

		if (nLeft == 0)
			break;

		size_t nBytes = nLeft * sizeof(*P);

		memmove(&P[dest], &P[src], nBytes);
	}

	int type = N_part_types - 1;

	if (Effect.Req.PartType[type] || !Task.Npart[type])
		Task.PartTotal -= Task.Npart[type];

	P = Realloc((void *)P, sizeof(*P) * Task.PartTotal);

	if (Task.Npart[0] == 0) 
		Free(Gas);

	rprintf(" done \n");

	return;
}

