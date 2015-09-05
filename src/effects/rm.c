#include "../globals.h"
#include "effects.h"

static const double faraday_prefac = // Dennison 1980 prefac = 2.6e-17 [cgs]
    e * e * e / (2 * pi * m_e * m_e * c * c * c * c);
static int cmp_part(const void *, const void *);

static struct faraday_data {
	int lastpart;
	double angle;
} rmData = { 0 };

/* Compute Faraday Rotation of particle */
void faraday_rotation(int ipart, double *result)
{
	result[0] = faraday_prefac * NumRho(ipart) * Gas[ipart].Bfld[2];

#ifdef RELATIVISTIC_RM_CORRECTION // Mirnov et al. 07, (50)
	result[0] *= (1 - 2 * k_B * Temperature(ipart) / (m_e * c * c));	 
#endif
	return;
}

 
void set_faraday_factors()
{
	return;
}

/* Faraday rotate pixel ipix of Q and U images
 * considering particle ipart.
 * Angle and sin/cos are computed on every new
 * particle only */
void faraday_rotate_pix(size_t ipart, size_t ipix, double weight)
{
	const int npix2 = p2(Param.XYPix);
	const double Q = Image[npix2 + ipix];
	const double U = Image[2 * npix2 + ipix];

	if (rmData.lastpart != ipart) {	// new particle - recompute        

		double dz = P[ipart].Mass / P[ipart].Rho
		    / (4 * p2(P[ipart].Hsml) );

		faraday_rotation(ipart, &rmData.angle);
		rmData.angle *= p2(c) / p2(Param.Freq) * dz * kpc2cgs;

		rmData.lastpart = ipart;
	}

	const double cos_2ang = cos(2 * weight * rmData.angle);
	const double sin_2ang = sin(2 * weight * rmData.angle);

	/* Mueller Matrix for rotation by $angle */
	Image[npix2 + ipix] = Q * cos_2ang - U * sin_2ang;

	Image[2 * npix2 + ipix] = Q * sin_2ang + U * cos_2ang;

	return;
}

void set_intrinsic_RM()
{
	Assert(Param.SynchroPAngInt == 0,
	       "Intrinsic RM and Pitch Angle Integration can't be combined");

	Assert(Task.NTask == 1, "Intrinsic RM works only on one CPU");

	rprintf("Using Intrinsic RM \n\n");

	set_faraday_factors();

	/* process particles bottom up - sort by z */
	qsort(P, Task.PartTotal, sizeof(*P), cmp_part);

	return;
}

/* comparison function for qsort */
static int cmp_part(const void *p1, const void *p2)
{
	const struct Particle_Data *P1 = (struct Particle_Data *)p1,
	    *P2 = (struct Particle_Data *)p2;

	return P1->Pos[2] > P2->Pos[2] ? 1 : 0;
}
