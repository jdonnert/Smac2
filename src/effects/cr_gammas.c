/* Compute Gamma Ray emission
 * from hadronic interactions
 * of CR protons with thermal
 * gas. (Pfrommer & Ensslin 2004)
 */

#include "../globals.h"
#include "effects.h"
#include <gsl/gsl_sf_gamma.h>

#define xi 2

#define mbarn2cgs (1e-3*barn2cgs)
#define mPi0_GeV (m_pi0*c*c/GeV2cgs)	/* Mass of Pi0 in GeV */
#define mPi0_eV (m_pi0*c*c/eV2cgs)	/* Mass of Pi0 in eV */
#define mp_GeV (m_p*c*c/GeV2cgs)	/* Mass of proton in GeV */

static double crp_number_density(int);

static double q_gamma_prefac, lambda_gamma_prefac, ncrp_prefac;

/* = q_gamma [1/cm^3/GeV] */
void crGamma_dif(int ipart, double *result)
{
	double n_N = NumRhoProton(ipart);
	double n_crp = crp_number_density(ipart);

	result[0] = q_gamma_prefac * n_N * n_crp;	/* [1/cm^3/GeV*GeV] */

	return;
}

/* = lambda_gamma */
void crGamma_int(int ipart, double *result)
{
	double n_N = NumRhoProton(ipart);
	double n_crp = crp_number_density(ipart);

	result[0] = lambda_gamma_prefac * n_N * n_crp;

	return;
}

static double crp_number_density(int ipart)
{
	double X_CRp = Param.X_crp;

#ifdef XCRP_SCALING
	double r = length3(P[ipart].Pos);
	X_CRp = pow(P[ipart].Rho / Cluster_central_density, XCRP_SCALING)
	    * exp(-r / CUTOFFSCALE);
#endif

	return X_CRp * ncrp_prefac * EpsNtherm(ipart);
}

void set_crGamma_factors()
{
	double a_crp = Param.a_crp;

	double a_gam = a_crp;
	double d_gam = 0.14 * pow(a_gam, -1.6) + 0.44;

	double sig_pp = 32 * (0.96 + exp(4.4 - 2.4 * a_gam)) * mbarn2cgs;

	double x1 = 1 / (1 + 0.5 * mPi0_eV / pow(Param.E_min, 2 * d_gam));
	double x2 = 1 / (1 + 0.5 * mPi0_eV / pow(Param.E_max, 2 * d_gam));

	double a = (a_gam + 1) / (2 * d_gam);
	double b = (a_gam - 1) / (2 * d_gam);

	ncrp_prefac = 1.0 / (m_p * c * c) * 2 * (a_crp - 1)
	    * pow(mp_GeV, a_crp - 1)
	    / gsl_sf_beta((0.5 * a_crp - 1), (3 - a_crp) / 2);

	q_gamma_prefac = sig_pp * c * pow(xi, 2. - a_gam)
	    * 4 / (3 * a_gam) * pow(mPi0_GeV, -1 * a_gam)
	    * pow(pow(2 * Param.E_min / (mPi0_eV), d_gam)
		  + pow(2 * Param.E_min / (mPi0_eV), -1 * d_gam),
		  -1 * a_gam / d_gam);

	lambda_gamma_prefac = sig_pp * pow(mPi0_GeV, 1 - a_gam)
	    * c / 3 * gsl_sf_beta(a, b)	/* Denormalise beta_inc */
	    *(gsl_sf_beta_inc(a, b, x2) - gsl_sf_beta_inc(a, b, x1))
	    / (pow(2, a_gam - 2) * a_gam * d_gam);

	return;
}

#undef xi

#undef mbarn2cgs
#undef mPi0_GeV
#undef mPi0_eV
#undef mp_GeV
