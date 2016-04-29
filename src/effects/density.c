#include "../proto.h"
#include "../globals.h"
#include "effects.h"

void Density(int ipart, double *result)	// 2D in kpc 
{
	result[0] = (double)(P[ipart].Rho) * Unit.Mass /
	    (Unit.Length * Unit.Length) / (Param.XYPix / Param.XYSize);

	return;
}

void NumDensity(int ipart, double *result)
{
	result[0] = (double)(P[ipart].Rho) * Unit.Mass /
	    (Unit.Length * Unit.Length) * n2ne / (u_mol * m_p);

	return;
}
