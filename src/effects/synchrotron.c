/* Computes synchroton emission from an arbitrary CR electron spectrum,
 * set in synchro_factors (Donnert & Brunetti 2014)
 * */

#include "../globals.h"
#include "effects.h"

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_synchrotron.h>
#include <gsl/gsl_spline.h>

#define NSTEP 128 // <1% rel error with power-law N(E)
#define TABLE_SIZE NSTEP

#define X_MAX 70.0 // -> F(x) < 2e-6
#define X_MIN 1e-20// -> F(x) < 5e-6

#define THETA_MAX 0.5*pi

double gsl_sf_synchrotron_1(double);	// GSL synchro kernel func 
double gsl_sf_synchrotron_2(double);
static inline double synchro_kernel(double, double *, double *);
static void prepare_kernel();
static double stokes2angle(double Q, double U);
extern void add_Ipol_chi_Pi();

static double (*spectrum_ptr) (double, int);

static double X_table[TABLE_SIZE] = { 0 },Kernel_Table1[TABLE_SIZE] = { 0 }, 
			  Kernel_Table2[TABLE_SIZE] = { 0 },
			  nu_c_prefac, j_nu_prefac, E_cntr_prefac;

static gsl_spline *Synchro_Spline1 = NULL, *Synchro_Spline2 = NULL;
static gsl_interp_accel *Acc[2] = { NULL };
#pragma omp threadprivate(Synchro_Spline1, Synchro_Spline2, Acc)

/* 
 * Return Synchrotron brightness Integrate over energy Ee 
 * (and pitch angle theta) using Trapezoidal rule.
 * see Longair, High Energy Astophysics 1994
 */

void synchrotron(int ipart, double *j_nu)
{
	const double nu = Param.Freq;

	double B = 0;		// calc B strength

	if (Param.SynchroPAngInt)
		B = length3(Gas[ipart].Bfld);
	else
		B = length2(Gas[ipart].Bfld);

	if (B == 0)		// catch boring ones
		return;

	for (int i = 0; i < MAXIMAGES; i++)
		j_nu[i] = 0;

	double E[NSTEP] = { 0 }, dE[NSTEP] = { 0 }, x[NSTEP] = { 0 }, 
		   F[NSTEP] =  { 0 }, F_mid[NSTEP] = { 0 };

#ifdef POLARISATION
	double F_para[NSTEP] = { 0 }, F_mid_para[NSTEP] = { 0 },
		   F_orth[NSTEP] = { 0 }, F_mid_orth[NSTEP] = { 0 };
#endif

	double E_min = sqrt(nu/(nu_c_prefac * B * X_MAX)); // [me*c]
	double E_max = sqrt(nu/(nu_c_prefac * B * X_MIN));

	double di = log(E_max/E_min)/(NSTEP-1);

	x[0] = X_MAX;
	E[0] = E_min; 
	dE[0] = E[0] - E_min * exp(-1*di);

	for (int i = 1; i < NSTEP; i++) {
	
		E[i] = E_min * exp(di * i);
		dE[i] = E[i] - E[i-1];
		x[i] = nu / (nu_c_prefac * p2(E[i]) * B);
	}

	double j_para = 0, j_orth = 0;

	for (int i = 1; i < NSTEP; i++) { // Simpson rule
	
		double K_orth = 0, K_para = 0;

		double N_E = (*spectrum_ptr)(E[i], ipart);
		double K = synchro_kernel(x[i], &K_orth, &K_para);

		F[i] = N_E * K;

#ifdef POLARISATION
		F_orth[i] = N_E * K_orth;
		F_para[i] = N_E * K_para;
#endif

		double x_mid = 0.5 * (x[i-1] + x[i]);
		double e_mid = 0.5 * (E[i-1] + E[i]);

		double N_E_mid = (*spectrum_ptr)(e_mid, ipart);
		double K_mid = synchro_kernel(x_mid, &K_orth, &K_para);

		F_mid[i] = N_E_mid * K_mid;

#ifdef POLARISATION
		F_mid_orth[i] = N_E_mid * K_orth;
		F_mid_para[i] = N_E_mid * K_para;
#endif

		j_nu[0] += dE[i] / 6.0 * (F[i] + F[i-1] + 4*F_mid[i]);

#ifdef POLARISATION
		j_orth += dE[i] / 6.0 * (F_orth[i] + F_orth[i-1] + 4*F_mid_orth[i]);
		j_para += dE[i] / 6.0 * (F_para[i] + F_para[i-1] + 4*F_mid_para[i]);
#endif
	}

	j_nu[0] *= j_nu_prefac * B;

#ifdef POLARISATION
	j_para *= 0.5 * j_nu_prefac * B;
	j_orth *= 0.5 * j_nu_prefac * B;

	const double bx = Gas[ipart].Bfld[0];
	const double by = Gas[ipart].Bfld[1];

	double sin_2chi = 0, cos_2chi = 0;

	if (bx || by) {

		sin_2chi = -2.0 * bx * by / (bx * bx + by * by);
		cos_2chi = (bx * bx - by * by) / (bx * bx + by * by);
	} 

	double j_pol = j_orth - j_para;

	j_nu[1] = j_pol * cos_2chi;	// Q 
	j_nu[2] = j_pol * sin_2chi;	// U 
#endif

	return;
}

/* choose CRe spectrum according to Spec */
void set_synchro_factors()
{
	/* if you change this, include prepare_synchrotron()  */
	prep_func_ptr = NULL;

	switch (Param.Effect_Flag) {

	case 0:		// High Energy Approx. Secondary (Brunetti 05) 
		set_brunetti_factors();
		spectrum_ptr = &cre_secondaries_brunetti;
		break;
	
	case 1:		// Read spectrum from file N(r) 
		prep_func_ptr = &set_tabulated_factors;
		spectrum_ptr = &cre_spectrum_tabulated;
		break;
	
	case 2:		// analytic power law spectrum
		set_powerlaw_factors(0);
		spectrum_ptr = &cre_spectrum_powerlaw;
		break;
	
	case 3:
		Assert(0, "Reading from Simulation not implemented yet");
		break;
	
	case 4:		// power law with norm from sim
		set_powerlaw_factors(1);
		spectrum_ptr = &cre_spectrum_powerlaw;
		break;
	
	case 10:		// Read compressed Spectrum
		prep_func_ptr = &setup_decompression;
		spectrum_ptr = &cre_spectrum_compressed;
		break;

	default:
		Assert(0, "Selected Effect_Flag %d not handled",
		       Param.Effect_Flag);
		break;

	}

	prepare_kernel();

	E_cntr_prefac = m_e*c*c * sqrt(2.0/3.0 * (2*pi * m_e*c)/e /0.29);
	nu_c_prefac = 3.0 * e / (4.0 * pi * pow(m_e, 3) * pow(c, 5));
	j_nu_prefac = e * e * e * sqrt(3) / (m_e * c * c);
	
	return;
}

/* Returns Synchrotron Kernel & Polarisations for x using 
 * limit formulas @ small and large energies. In between the Integral is
 * tabulated.
 */

static inline double synchro_kernel(double x, double *orthPol, double *paraPol)
{
	double Fx = gsl_spline_eval(Synchro_Spline1, x, Acc[0]);

#ifdef POLARISATION

	double Gx = gsl_sf_synchrotron_2(x); //gsl_spline_eval(Synchro_Spline2, x, Acc[1]);

	*orthPol = (Fx + Gx);	// Rybicki & Lightman p179 
	*paraPol = (Fx - Gx);
	
#endif

	return Fx;
}

void prepare_kernel()
{
	double dx = log(X_MAX*1.1/X_MIN/0.9)/(TABLE_SIZE-1);

	#pragma omp parallel
	{
		
	#pragma omp for
	for (int i = 0; i < TABLE_SIZE; i++) {
		
		X_table[i] = X_MIN*0.9 * exp(dx*i);
	
		Kernel_Table1[i] = gsl_sf_synchrotron_1(X_table[i]);
		Kernel_Table2[i] = gsl_sf_synchrotron_2(X_table[i]);
	}

	Acc[0] = gsl_interp_accel_alloc();
	Acc[1] = gsl_interp_accel_alloc();

	Synchro_Spline1 = gsl_spline_alloc(gsl_interp_cspline, TABLE_SIZE);
	Synchro_Spline2 = gsl_spline_alloc(gsl_interp_cspline, TABLE_SIZE);

	gsl_spline_init(Synchro_Spline1, X_table, Kernel_Table1, TABLE_SIZE);
	gsl_spline_init(Synchro_Spline2, X_table, Kernel_Table2, TABLE_SIZE);

	} // omp parallel

	return;
}

/* 
 * Convert stokes parameter to polarisation angle.
 * angle is defined [-pi,pi]
 */

static double stokes2angle(double Q, double U)
{
	return atan2(U, Q) * 0.5;
}

extern void add_Ipol_chi_Pi()
{
	const int npix2 = p2(Param.XYPix);
	size_t i, iQ = npix2, iU = 2 * npix2, iIpol = 3 * npix2,
	    iChi = 4 * npix2, iPi = 5 * npix2;

	rprintf("\nAdding <3> images: Ipol, chi, Pi \n\n");

	Effect.Nimage = 6;

	size_t nBytes = Effect.Nimage * npix2 * sizeof(*Image);

	Image = Realloc(Image, nBytes);

	/* Fill new images */
	for (i = 0; i < npix2; i++) {
		/* pol. Intensity */
		Image[iIpol] = sqrt(Image[iQ] * Image[iQ] + Image[iU] * Image[iU]);
		/* pol. Angle */
		Image[iChi++] = stokes2angle(Image[iQ++], Image[iU++]);
		/* pol. Degree */
		Image[iPi++] = Image[iIpol++] / Image[i];
	}

	return;
}

#undef X_MAX
#undef THETA_MAX
#undef NSTEP
