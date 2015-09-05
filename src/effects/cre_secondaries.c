/* Computes CR electron spectrum 
 * from Brunetti&Blasi 2005
 * (Crossection High Energy approximation) 
 * */

#include "../globals.h"
#include "effects.h"

#define MIN 0.0			/*Integrals parameter */
#define MAX 1.0
#define NSTEP 500000.0

#define mMu m_mu *c*c*1e-6/eV2cgs	/* [MeV] */
#define mPi m_pi *c*c*1e-6/eV2cgs	/* [MeV] */

#define c1 1.22
#define c2 0.92
#define s Param.a_crp

void compute_P1P2P3();
static double P1, P2, P3;
static double a, A_s, normPrefac, NcrpPrefac;

/* we integrate in electron Energy 
 */
double cre_secondaries_brunetti(double Ee, int ipart)
{
	if (Ee < 1e2 * m_e * c * c)	/* High energy cross section not valid here */
		return (0);

	double Babs2 = length3_sq(Gas[ipart].Bfld);
	double X_CRp = Param.X_crp;

#ifdef XCRP_SCALING
	double r = length3(P[ipart].Pos);
	X_CRp = pow(P[ipart].Rho / Cluster_central_density, XCRP_SCALING)
	    * exp(-r / CUTOFFSCALE);
#endif

	double N_crp = X_CRp * NcrpPrefac * EpsNtherm(ipart);

	double norm = normPrefac / (Babs2 + Bcmb * Bcmb);

	double pwLaw_corr = P1 / (s - 1)
	    + P2 / (s - 1) * (log(a * Ee / (6.4e9 * eV2cgs)) + 1 / (s - 1))
	    - P3 / (1.5 - s) * sqrt(a * Ee / (1e9 * eV2cgs));

	/* e+ & e- => factor 2 */
	double Ne_Brunetti = 2 * N_crp * NumRhoProton(ipart)
	    * norm * pow(Ee, -s - 1) * pwLaw_corr;

	if (Ne_Brunetti < 0)
		Ne_Brunetti = 0;

	return (Ne_Brunetti);
}

void set_brunetti_factors()
{
	compute_P1P2P3();

	a = 2 * mPi * mPi / (mPi * mPi + mMu * mMu);
	A_s = 2 * c * sigma_pp * pow(a, 1.0 - s);

	normPrefac = 3.0 / 4.0 * A_s * (m_e * m_e)
	    * (c * c * c) / sigma_T * 8 * pi;
	NcrpPrefac = (s - 2) / pow(Param.E_min * eV2cgs, 2 - s);

	return;
}

/* The spectrum is rel. 
 * sensitive to these 3 
 * parameters. Brunetti 05 
 * does not define them 
 * accurate enough.
 */
void compute_P1P2P3()
{
	int i;

	double m_fac1, m_fac2;
	double h, aa, bb, cc, P_fac1, P_fac2, P_fac3;
	double x, f1, f2, I0 = 0, I1 = 0, I2 = 0;

	m_fac1 = 2 * mPi * mPi / (mMu * mMu + mPi * mPi);
	m_fac2 = (m_fac1 * m_fac1 - 1) * (1 + (mMu / mPi) * (mMu / mPi))
	    / (1 - (mMu / mPi) * (mMu / mPi));

	aa = 5.0 / 12 * (1 + 0.2 * m_fac2);
	bb = -3.0 / 4 * (1 + m_fac2);
	cc = 1.0 / 3 * (1 + 2.0 * m_fac2);

	P_fac1 = aa / s / s + bb / (s + 2) / (s + 2) + cc / (s + 3) / (s + 3);
	P_fac2 = aa / s + bb / (s + 2) + cc / (s + 3);
	P_fac3 = 1.5 * (aa / (s + 0.5) + bb / (s + 2.5) + cc / (s + 3.5));

	h = fabs(MIN - MAX / NSTEP);

	for (i = 0; i < NSTEP; i++) {
		x = MIN + i * h + 0.5 * h;
		f1 = c1 * pow(1.0 - x, 3.5) + c2 * exp(-18.0 * x);
		f2 = pow(x, s - 2.0);
		I0 += h * f1 * f2;
		I1 += h * f1 * f2 * log(x);
		I2 += h * f1 * pow(x, s - 1.5);
	}
	P1 = P_fac1 * I0 - P_fac2 * I1;
	P2 = P_fac2 * I0;
	P3 = P_fac3 * I2;

	rprintf("Using Brunetti's Secondaries:\n"
		"   s    = %1.4f\n"
		"   P1   = %1.6f\n"
		"   P2   = %1.6f\n" "   P3   = %1.6f\n\n", s, P1, P2, P3);

	return;
}

/* compute mean central density inside max_radius kpc */
void find_central_density()
{
	const double max_radius = 100;

	int nPart = 0;
	double rho = 0;

	for (int ipart = 0; ipart < Task.Npart[0]; ipart++) {

		double r = length3(P[ipart].Pos);

		if (r < max_radius) {	// this assumes the snapshot is centered 

			rho += P[ipart].Rho;

			nPart++;
		}
	}

	double dbl_comm_buf = 0;
	MPI_Allreduce(&rho, &dbl_comm_buf, 1,
		      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	rho = dbl_comm_buf;

	int int_comm_buf = 0;
	MPI_Allreduce(&nPart, &int_comm_buf, 1,
		      MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	nPart = int_comm_buf;

	Cluster_central_density = rho / nPart;
#ifdef XCRP_SCALING
	rprintf("Enabling CRp scaling with density\n"
		"    Scaling coefficient    = %g \n"
		"    Cutoff scale           = %g \n"
		"    central density        = %g \n\n",
		XCRP_SCALING, (double)CUTOFFSCALE, Cluster_central_density);
#endif

	return;
}

#undef MIN
#undef MAX
#undef NSTEP
#undef mMu
#undef mPi
#undef c1
#undef c2
#undef s
