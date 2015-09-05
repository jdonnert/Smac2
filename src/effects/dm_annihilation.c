#include "../globals.h"
#include "effects.h"

void dm_annihilation(int ipart, double *result)
{
	dm_density(ipart, result);

	result[0] *= result[0];

	return;
}

void dm_density(int ipart, double *result)
{
	result[0] = (double)(P[ipart].Rho) * Unit.Mass /
	    (Unit.Length * Unit.Length) / (Param.XYPix / Param.XYSize);

	return;
}
