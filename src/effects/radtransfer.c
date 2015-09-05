#include "../globals.h"
#include "effects.h"

void rt_Hydrogen(int ipart, double *result)
{
#ifdef RADTRANSFER
	result[0] = Gas[ipart].nHII;
#endif
	return;
}

void rt_HeliumII(int ipart, double *result)
{
#ifdef RADTRANSFER
	result[0] = Gas[ipart].nHeII;
#endif
	return;
}

void rt_HeliumIII(int ipart, double *result)
{
#ifdef RADTRANSFER
	result[0] = Gas[ipart].nHeIII;
#endif
	return;
}
