
#include "../proto.h"
#include "../globals.h"
#include "effects.h"

/* Mean free path from COulomb scattering, Sarazin 1988 */
void coul_mfp(int ipart, double *result)
{
	float T = Temperature(ipart);

	result[0] = 23. * kpc2cgs * p2(T/1e8) / (NumRho(ipart) / 1e3);

	return;
}

void coul_coll(int ipart, double *result)
{
	double mfp;
	coul_mfp(ipart, &mfp);

	result[0] = 2 * P[ipart].Hsml * kpc2cgs / mfp;

	return;
}
