/* Code Units */

void set_units();

struct units {			/* For unit Conversion. */
	double Length;
	double Mass;
	double Vel;
	double Time;
	double Energy;
	double Density;
} Unit, Comv2phys;

extern double Temperature(const int ipart);
extern double T_keV(int);
extern double NumRho(int);	// don't consider 25% Helium 
extern double NumRhoProton(int);	// = n_th 
extern double EpsNtherm(int);
extern double Epsilon_cre_sim(int);

#define Rho(i) ((double)(P[i].Rho) * Unit.Density)
