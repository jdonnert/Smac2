/* Compute reacceleration coefficient */

#include "../globals.h"
#include "effects.h"

#include <gsl/gsl_sf_expint.h>

#define ALPHA -2./3.
#define MAX_EDDY_SIZE 200	/* [kpc] */
#define MIN_EDDY_SIZE 0.1
#define GAMMA adiabatic_index

static double dpp_prefac, tauc_prefac;
double gsl_sf_expint_E1(double x);

void set_Dpp_factors()
{
	dpp_prefac = 4.45 * pi * pi / sqrt(32 * pi * pi * pi) / c;
	tauc_prefac = sqrt(2 * pi * Param.Freq * m_e * c / e);

	if (Param.Effect_Flag == 2)
		rprintf("B = 1muG -> gamma = %g\n", tauc_prefac / sqrt(1e-6));

	return;
}

/* Donnert et al. 2013, Brunetti & Lazarian 2007 */
void Dpp(int ipart, double *result)
{
#ifdef VTURB

	float scale = 2 * P[ipart].Hsml * Unit.Length;

	float vturb = Gas[ipart].VRms * Unit.Vel;

	double csound2 = 0;
	Vsound(ipart, &csound2);
	csound2 *= csound2;
	
	result[0] =
	    9e-8 * Param.Eta_t * p2(vturb) * p2(vturb) / (scale * csound2);
#endif

	return;
}

/* Cassano et al. 2006, eqn 36 */
void Tau_reacc(int ipart, double *result)
{
#ifdef VTURB
	double Dpp_p2;

	Dpp(ipart, &Dpp_p2);

	result[0] = 0.25 / Dpp_p2;
#endif

	return;
}

/* Brunetti & Blasi 2004, eqn 1-3 */
void Tau_cooling(int ipart, double *result)
{
	float B = length3(Gas[ipart].Bfld) / 3.2e-6;

	float gam300 = tauc_prefac / sqrt(B) / 300;	/* Longair 18.51 */

	float nth3 = NumRho(ipart) * 1e3;

	float tau_i = nth3 / gam300 * (1.2 + 1 / 75 * log(gam300 / nth3));

	float tau_rad = 1 / 3 * gam300 * (B * B + pow(1 + Snap.Redshift, 4));

	result[0] = 4 / (tau_i + tau_rad) * Gyr2sec;

	return;
}
