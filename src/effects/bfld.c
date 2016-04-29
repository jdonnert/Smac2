#include "../proto.h"
#include "../globals.h"
#include "effects.h"

void bfld(int ipart, double *result)
{
	result[0] = length3(Gas[ipart].Bfld);

	return;
}

void bfld_orth(int ipart, double *result)
{				/* We always integrate along z */
	result[0] = (sqrt(Gas[ipart].Bfld[0] * Gas[ipart].Bfld[0]
			  + Gas[ipart].Bfld[1] * Gas[ipart].Bfld[1]));

	return;
}

void alven_vel(int ipart, double *result)
{
	result[0] = length3(Gas[ipart].Bfld) / sqrt(4 * pi * Rho(ipart));

	return;
}

void divB(int ipart, double *result)
{
#ifdef DIVB
	result[0] =
	    fabs(Gas[ipart].DivB) * P[ipart].Hsml / length3(Gas[ipart].Bfld);
#endif
	return;
}

/* P_therm / P_magnetic */
void plasma_beta(int ipart, double *result)
{
	result[0] =
	    Gas[ipart].U * Rho(ipart) * (adiabatic_index -
					 1) * Unit.Mass / p2(Unit.Time) /
	    (p2(length3(Gas[ipart].Bfld)) / (4 * pi));

	return;
}
