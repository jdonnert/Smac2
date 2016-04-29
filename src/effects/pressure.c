#include "../proto.h"
#include "../globals.h"
#include "effects.h"

void Therm_Pressure(int ipart, double *result)
{
	result[0] = Gas[ipart].U * Rho(ipart) * (adiabatic_index - 1)
	    * Unit.Mass / p2(Unit.Time);

	return;
}

void Magn_Pressure(int ipart, double *result)
{
	result[0] = p2(length3(Gas[ipart].Bfld)) / (4 * pi);

	return;
}

void Turb_Pressure(int ipart, double *result)
{
#ifdef VTURB
	result[0] =
	    0.5 * Rho(ipart) * Gas[ipart].VRms * Gas[ipart].VRms * Unit.Vel *
	    Unit.Vel;
#endif

	return;
}
