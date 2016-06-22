/* 
 * For a good description see Longair 1981 chapter 8. We implement 8.86
 * Note that for power laws with cut offs this formula breaks down for 
 * large values of B because it is calculated over ]0,infty[.
 * The CRe spectrum is normalised relative to the energy density of the
 * thermal plasma at a normalisation energy. eps_therm = \int E N(E) dE, 
 */

#include "../proto.h"
#include "../globals.h"
#include "effects.h"

#include <gsl/gsl_sf_gamma.h>

#define nu Param.Freq

static double cre_spec_norm_parfile(int ipart);
static double cre_spec_idx_parfile(int ipart);
static double cre_spec_norm_particle(int ipart);
static double cre_spec_idx_particle(int ipart);

static double (*cr_energy_density) (int);
static double (*cre_spec_norm) (int);
static double (*cre_spectral_index) (int);

static double nufac = 0;

void aSynchrotron(int ipart, double *j_nu)	// [erg/cm^3/Hz] 
{
	const double bx = Gas[ipart].Bfld[0];
	const double by = Gas[ipart].Bfld[1];

	double B = 0;

	if (Param.SynchroPAngInt) 
		B = length3(Gas[ipart].Bfld);
 	else 
		B = sqrt(bx*bx + by*by);

	if (B == 0)
		B = 3e-6;

	if (Gas[ipart].NMach < 1.3) {
		
		j_nu[0] = j_nu[1] = j_nu[2] = 0;
	
		return ;
	}

	double n0 = (*cre_spec_norm) (ipart) * (*cr_energy_density) (ipart);
	double s = (*cre_spectral_index)(ipart);
		
	double prefac = sqrt(3)*p3(e)/(m_e * p2(c)*(s + 1)) //Longair eq 8.128
	    * gsl_sf_gamma(s / 4 + 19. / 12.) 
		* gsl_sf_gamma(s / 4. - 1. / 12.)
	    * pow(nufac , (s - 1)/2);

	if (Param.SynchroPAngInt) // Longair eq 8.87
		prefac *= 0.5 * sqrt(pi) * gsl_sf_gamma(s / 4 + 5. / 4.)
		    / gsl_sf_gamma(s / 4. + 7. / 4.);

	j_nu[0] = prefac * pow(B, 0.5 * (s + 1)) * n0;
	
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
		
		cre_spec_norm = &cre_spec_norm_parfile;
		
		cre_spectral_index = &cre_spec_idx_parfile;

		break;

	case 1:
		rprintf("Using CR electrons normalisation \n"
			"from simulation, setting X_crp = 1 \n");

		cr_energy_density = &Epsilon_cre_sim;	// see unit.c 

		Param.X_crp = 1;

		cre_spec_norm = &cre_spec_norm_parfile;

		cre_spectral_index = &cre_spec_idx_parfile;

		break;

	case 2:
		rprintf("CR electrons from particles with standard relic model");
	
		cr_energy_density = &EpsNtherm;
	
		cre_spec_norm = &cre_spec_norm_particle;

		cre_spectral_index = &cre_spec_idx_particle;

		break;

	default:
		Assert(0, "Effect Flag  not handled");
		break;
	}

	nufac = (3*e)/(p3(m_e) * c * c * c * c * c * 2 * pi * nu);

	return;
}

static double cre_spec_norm_parfile(int ipart)
{
	return Param.X_crp * (Param.a_crp - 2) * 
			pow(Param.E_min * eV2cgs, Param.a_crp - 2);
};

static double cre_spec_idx_parfile(int ipart)
{
	return Param.a_crp;
}

static double cre_spec_norm_particle(int ipart) // Kang, Ryu, Cen, Ostriker 07
{
	double  M = Gas[ipart].NMach;

	double norm = 0;

	if (M < 2)
		return 0;

	if (M < 2)
		norm = 1.96e-3 * (M*M-1); // eq. A2
	else { // M>2

		norm = 5.46 - 9.78 * (M-1) + 4.17*p2(M-1) - 0.334 * p3(M-1) 
		      + 0.57 * p2(p2(M-1)); //eq. A5+6
		
		norm /= M*M*M*M;
	}

	return norm / 7.6e14; // Donnert et al 2016, eq. 40, p_0 = 0.1 me c
}

static double cre_spec_idx_particle(int ipart)
{
	double  M = Gas[ipart].NMach;

	double s = 0;

	if (M < 1)
		s = 0; // no emission
	else
		s = 2* (M*M + 1) / (M*M-1);

	return s;
}


#undef s
#undef nu
