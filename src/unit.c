/* Code Units 
 * Here 3 things have to be done :
 * We need to provide a way to make things 
 * not comoving. Use
 * Comv2phys, which is set in the 
 * according readin routine. 
 * We need to convert things to cgs. For 
 * most variables cgs values cant be stored
 * in a float and taking a double takes to
 * much memory.  
 * Every readin module should 
 * provide helper functions for "complicated" 
 * things like density in a way like Foo(ipart).
 * Further the Unit structure contains the
 * primitive values specified in the parameter
 * file. 
 * */

#include "globals.h"

struct units Unit;
struct units Comv2phys;

double Epsilon_cre_sim(int);

static double Temperature_Prefac = 0;

void set_units()
{
	Unit.Time = Unit.Length / Unit.Vel;
	Unit.Energy = Unit.Mass * Unit.Vel * Unit.Vel;
	Unit.Density = (Unit.Mass / Unit.Length) / Unit.Length / Unit.Length;

	rprintf("Setting System of Units \n"
		"   Unit Length     = %g cm \n"
		"   Unit Time       = %g sec\n"
		"   Unit Mass       = %g g  \n"
		"   Unit Vel        = %g cm/s\n"
		"   Unit Energy     = %g erg\n"
		"   Unit Density    = %g g/cm^3\n"
		"   adiabatic index = %g\n", Unit.Length, Unit.Time, Unit.Mass,
		Unit.Vel, Unit.Energy,Unit.Density, adiabatic_index);

	if (adiabatic_index != 5. / 3.)
		rprintf("\n !!  USING NON STANDARD GAMMA !! \n\n");

	Temperature_Prefac = (adiabatic_index-1) * p2(Unit.Vel) * m_p * 
		mean_mol_weight / k_B;

	return;
}

#if defined I_GADGET2 || defined I_USER || defined I_AREPO || defined I_PIERNIK
/* Gadget Default Units:  cm/kpc, g/(1e10 M_sol), cm/sec/(km/sec), */

extern double NumRho(int ipart) // particle number density
{				
	return (Rho(ipart) * n2ne / (u_mol * m_p));
}

extern double NumRhoProton(int ipart) // proton number density 
{				
	return (Rho(ipart) * H_frac / m_p);
}

extern double Temperature(const int ipart)
{
	return Gas[ipart].U * Temperature_Prefac;
}

extern double T_keV(int ipart)
{
	return ((adiabatic_index - 1) * Gas[ipart].U * Unit.Vel * Unit.Vel * m_p *
		mean_mol_weight / eV2cgs / 1000);
}

extern double EpsNtherm(int ipart) // Thermal energy density in cgs 
{				
	return (Rho(ipart) / (m_p * u_mol) * k_B * Temperature(ipart));
}

extern double Epsilon_cre_sim(int ipart) // CR energy density 
{				
#ifdef I_PIERNIK
	return (Gas[ipart].ECR * Unit.Energy
		/ Unit.Length / Unit.Length / Unit.Length);
#else
	Assert(0, "CR energy Density from Sim not implemented");

	return (-1);
#endif
}

#endif
