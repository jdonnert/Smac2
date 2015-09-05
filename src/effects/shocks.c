#include "../globals.h"
#include "effects.h"

void Machnumber(int ipart, double *result)
{
#ifdef SHOCKS
	if (Gas[ipart].NMach < 1)	/* Mach < 1 not interesting here */
		result[0] = 0;
	else
		result[0] = Gas[ipart].NMach;
#endif
	return;
}

void AlphaVisc(int ipart, double *result)
{
#ifdef TIME_DEP_ART_VISC
	result[0] = Gas[ipart].AlphaVisc;
#endif
	return;
}
