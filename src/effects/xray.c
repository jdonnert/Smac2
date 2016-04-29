/* Compute X-Ray surface brightness
 * of particle between [Emin, Emax]
 * Bartelmann & Steinmetz 1996 eqn 3
 * */

#include "../proto.h"
#include "../globals.h"
#include "effects.h"

#define C_j 2.42E-24

static double prefac = 4.0 * C_j / (1 + H_frac);

void Xray(int ipart, double *result)
{
	const float T_eV = Temperature(ipart)/ eV2cgs;

	result[0] = prefac * sqrt(k_B * T_eV / 1000 )
	    * NumRho(ipart) * NumRho(ipart)
	    * (exp(-Param.E_min / (k_B * T_eV ))
         - exp(-Param.E_max / (k_B * T_eV )));

	return;
}
#undef C_j
