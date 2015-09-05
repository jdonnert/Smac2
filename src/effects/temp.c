#include "../globals.h"
#include "effects.h"

void Temp(int ipart, double *result)
{
	result[0] = Temperature(ipart);

	return;
}

void Soundspeed(int ipart, double *result)
{
	result[0] = sqrt(adiabatic_index * (adiabatic_index - 1)
			 * Gas[ipart].U) * Unit.Vel;

	return;
}
