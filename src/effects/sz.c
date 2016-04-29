/* Computes various flavours
 * of the Sunyaev Zeldovich 
 * effect (e.g. Springel,White,Hernquist 2001)
 * */

#include "../proto.h"
#include "../globals.h"
#include "effects.h"

static double x, yPrefac, tSzPrefac, kSzPrefac;

void comptonY(int ipart, double *result)	// = SMAC1 DT/T 
{
	result[0] = yPrefac * NumRho(ipart) * (Temperature(ipart) - Tcmb);

	return;
}

void tSz(int ipart, double *result)
{
	comptonY(ipart, result);

	result[0] *= tSzPrefac;

	return;
}

void kSz(int ipart, double *result)	// We always project along z
{
	result[0] = kSzPrefac * (NumRho(ipart)) * P[ipart].Vel[2] * Unit.Vel;

	return;
}

void set_sz_factors(bool DI_over_I)
{
	x = h_planck * Param.Freq / (k_B * Tcmb);	// Dimensionless Frequency 

	yPrefac = sigma_T * k_B / (m_e * c * c);
	tSzPrefac = (x * (exp(x) + 1) / (exp(x) - 1) - 4);
	kSzPrefac = -1 * sigma_T / c;

	if (DI_over_I) {
		tSzPrefac *= exp(x) - 1 / (x * exp(x));
		kSzPrefac *= exp(x) - 1 / (x * exp(x));
	}

	return;
}
