/* 
 * For a good description see Longair 1981 chapter 8. We implement 8.86
 * Note that for power laws with cut offs this formula breaks down for 
 * large values of B because it is calculated over ]0,infty[.
 * The CRe spectrum is normalised relative to the energy density of the
 * thermal plasma at a normalisation energy. eps_therm = \int E N(E) dE, 
 */

#include "../globals.h"
#include "effects.h"

#include <gsl/gsl_sf_gamma.h>

#define s Param.a_crp
#define nu Param.Freq

static double prefac = 0, cre_spec_norm = 0;

double (*cr_energy_density) (int);

void aSynchrotron(int ipart, double *j_nu)	// [erg/cm^3/Hz] 
{
	const double bx = Gas[ipart].Bfld[0];
	const double by = Gas[ipart].Bfld[1];

	double B = 0;

	if (Param.SynchroPAngInt) 
		B = length3(Gas[ipart].Bfld);
 	else 
		B = sqrt(bx*bx + by*by);

	const double eps_cre = cre_spec_norm * (*cr_energy_density) (ipart);
	
	j_nu[0] = prefac * pow(B, 0.5 * (s + 1)) * eps_cre;

	double sin_2chi, cos_2chi;

	if ((bx != 0) || (by != 0)) {	// Kotarba 2010 

		sin_2chi = -2.0 * bx * by / (bx * bx + by * by);
		
		cos_2chi = (bx * bx - by * by) / (bx * bx + by * by);
	
	} else {

		sin_2chi = cos_2chi = 0;
	}

	double j_pol = j_nu[0] * (s + 1) / (s + 7. / 3.);	// Longair 1981 18.59 

	j_nu[1] = j_pol * cos_2chi;	// Q 
	j_nu[2] = j_pol * sin_2chi;	// U     

	if (Param.SynchroPAngInt == 1)
		j_nu[1] = j_nu[2] = 0;
	
	return;
}

void set_aSynchro_factors()
{
	switch (Param.Effect_Flag) {
	case 0:
		rprintf("Using CR electrons normalisation "
			"relative to thermal \n\n");

		cr_energy_density = &EpsNtherm;
		break;
	case 1:
		rprintf("Using CR electrons normalisation \n"
			"from simulation, setting X_crp = 1 \n");

		cr_energy_density = &Epsilon_cre_sim;	// see unit.c 

		Param.X_crp = 1;

		break;
	default:
		Assert(0, "Effect Flag  not handled");
		break;
	}

	cre_spec_norm = Param.X_crp * (s - 2) * pow(Param.E_min * eV2cgs, s - 2);

	prefac = sqrt(3) * p3(e) / (m_e * c * c * (s + 1)) // Longair eq 8.128
	    * gsl_sf_gamma(s / 4 + 19. / 12.) * gsl_sf_gamma(s / 4. - 1. / 12.)
	    * pow( (3*e)/(p3(m_e) * c * c * c * c * c * 2 * pi * nu), (s - 1)/2);

	if (Param.SynchroPAngInt) // Longair eq 8.87
		prefac *= 0.5 * sqrt(pi) * gsl_sf_gamma(s / 4 + 5. / 4.)
		    / gsl_sf_gamma(s / 4. + 7. / 4.);

	return;
}

#undef s
#undef nu
