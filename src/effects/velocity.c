/* All kinds of crazy velocities */
#include "../globals.h"
#include "effects.h"

void Vel(int ipart, double *result)
{
	result[0] = length3(P[ipart].Vel) * Unit.Vel;
	return;
}

void Vel_perp(int ipart, double *result)
{
	result[0] = sqrt(p2(P[ipart].Vel[0]) + p2(P[ipart].Vel[1]));

	return;
}

void Vsound(int ipart, double *result)
{
	result[0] = sqrt(adiabatic_index * (adiabatic_index - 1)
			 * Gas[ipart].U) * Unit.Vel;
	return;
}

void Vorticity(int ipart, double *result)
{
#ifdef VORTICITY
	result[0] = length3(Gas[ipart].Vorticity);
#endif
	return;
}

void Kernel_weighted_velocity(int ipart, double *result)
{
#ifdef VSPH
	result[0] = length3(Gas[ipart].VelSPH) * Unit.Vel;
#endif
	return;
}

void Kernel_weighted_bulk_velocity(int ipart, double *result)
{
#ifdef VSPH
	result[0] = length3(Gas[ipart].VelBulk) * Unit.Vel;
#endif
	return;
}

void Kernel_weighted_velocity_dispersion(int ipart, double *result)
{
#ifdef VSPH
	result[0] = Gas[ipart].VelTurb * Unit.Vel;
#endif
	return;
}

void Turbulent_energy(int ipart, double *result)
{
#ifdef VSPH
	result[0] = 0.5* p2(Gas[ipart].VelTurb * Unit.Vel) 
		* P[ipart].Rho * Unit.Density;
#endif
	return;
}

void Velt(int ipart, double *result)
{
#ifdef VTURB
	result[0] = Gas[ipart].VTurb * Unit.Vel;
#endif
	return;
}

void Vrms(int ipart, double *result)
{
#ifdef VTURB
	result[0] = Gas[ipart].VRms * Unit.Vel;
#endif
	return;
}

void Vbulk(int ipart, double *result)
{
#ifdef VTURB
	result[0] = length3(Gas[ipart].VBulk) * Unit.Vel;
#endif
	return;
}

void Vdiv(int ipart, double *result)
{
#ifdef VTURB
	result[0] = Gas[ipart].DivVel * Unit.Vel / Unit.Length;
#endif
	return;
}

void Vrot(int ipart, double *result)
{
#ifdef VTURB
	result[0] = Gas[ipart].CurlVel * Unit.Vel / Unit.Length;
#endif
	return;
}

void Velt_scaled(int ipart, double *result)
{
#ifdef VTURB
	double khsml = 0.5 / (P[ipart].Hsml * Unit.Length);

	double k0 = 1 / (Param.TurbScale * Unit.Length);

	/* Assume Kolmogorov Law -5/3 */
	Velt(ipart, result);

	result[0] *= pow(k0 / khsml, -1 / 3.);
#endif
	return;
}

void Vrms_scaled(int ipart, double *result)
{
#ifdef VTURB
	double khsml = 0.5 / (P[ipart].Hsml * Unit.Length);

	double k0 = 1 / (Param.TurbScale * Unit.Length);

	/* Assume Kolmogorov Law -5/3 */
	Vrms(ipart, result);

	result[0] *= pow(k0 / khsml, -1 / 3.);
#endif
	return;
}

void Vtan(int ipart, double *result)
{
#ifdef VTURB
	result[0] = Gas[ipart].VTan * Unit.Vel;
#endif
	return;
}

void Vrad(int ipart, double *result)
{
#ifdef VTURB
	result[0] = Gas[ipart].VRad * Unit.Vel;
#endif

	return;
}
