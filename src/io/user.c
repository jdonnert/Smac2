/* 
 * User defined Pseudo Snapshot for testing purposes 
 */

#include "../proto.h"
#include "../globals.h"

#ifdef I_USER

#define NGAS 1
#define BOX 100

extern void read_snapshot(char *filename)
{
	int ngas[6] = { NGAS, 0, 0, 0, 0, 0 };

	start_timing(CPU_READIN);

	Assert(Task.NTask == 1, "USER input only on one CPU \n");

	rprintf("\nUsing USER input: <%d> gas particles \n\n", ngas[0]);

	/* Header Stuff */
	Task.PartTotal = 0;
	Task.Npart[0] = Task.Npart[1] = Task.Npart[2] = 0;
	Task.Npart[3] = Task.Npart[4] = Task.Npart[5] = 0;

	Reallocate_P(NGAS, ngas, +1);

	Snap.SnapNum = 0;
	Snap.Boxsize = BOX;
	Snap.Redshift = 0;
	Snap.Time = 1;
	Snap.PartTotal = Task.PartTotal;
	Snap.Masstab[0] = 1;
	for (int type = 0; type < N_part_types; type++)
		Snap.Npart[type] = Task.Npart[type];

	/* Set Comoving units for Gadget */
	Comv2phys.Length = 1 ; /// (1 + Snap.Redshift) / Cosmo.h;
	Comv2phys.Mass = 1 ; /// Cosmo.h;
	Comv2phys.Vel = 1; //sqrt(1. / (1 + Snap.Redshift));

	/* Particle 0 */
	int ipart = 0;

	P[ipart].Pos[0] = 0;
	P[ipart].Pos[1] = 0;
	P[ipart].Pos[2] = 0;
	P[ipart].Vel[0] = 0;
	P[ipart].Vel[1] = 0;
	P[ipart].Vel[2] = 0;
	P[ipart].ID = ipart + 1;
	P[ipart].Type = 0;
	P[ipart].Mass = Snap.Masstab[0];
	P[ipart].Hsml = 1.2;
	P[ipart].Rho = 2.8075881e-06;

	Gas[ipart].U = 2104884.2;
	Gas[ipart].Bfld[0] = 5e-6;
	Gas[ipart].Bfld[1] = 0;
	Gas[ipart].Bfld[2] = 0;

#ifdef VTURB
	Gas[ipart].VTurb = 0;
	Gas[ipart].VRms = 0;
	Gas[ipart].VBulk[0] = 0;
	Gas[ipart].VBulk[1] = 0;
	Gas[ipart].VBulk[2] = 0;
	Gas[ipart].TNgb = 64;
	Gas[ipart].Dpp = 1e-22;
#endif

	stop_timing(CPU_READIN);

	return;
}

#endif
