/* Artifiscial CR electron power law spectrum   */

#include "../proto.h"
#include "../globals.h"
#include "effects.h"

#define s Param.a_crp

static double Nprefac;
static double (*cr_energy_density) (int);

/* N(E) dE*/
double cre_spectrum_powerlaw(double Ee, int ipart)
{
	if (Ee < Param.E_min * eV2cgs)
		return 0;
	
	return Nprefac * (*cr_energy_density) (ipart) * pow(Ee, -s);
}

void set_powerlaw_factors(int flag)
{
	Nprefac = (s - 2) * Param.X_crp * pow(Param.E_min * eV2cgs, s - 2);

	switch (flag) {

	case 0:
		rprintf("Using CRe powerlaw\n"
			"Normalisation rel. thermal \n" "s  = %1.4f\n"
			"NPrefac = %g \n\n", s, Nprefac);
		cr_energy_density = &EpsNtherm;
		break;
	case 1:
		rprintf("Using CRe powerlaw\n"
			"Normalisation from simulation \n"
			"Setting X_cr = 1 \n" "s  = %1.4f\n\n", s);
		cr_energy_density = &Epsilon_cre_sim;
		Param.X_crp = 1;
		break;
	}


	return;
}

#undef s
