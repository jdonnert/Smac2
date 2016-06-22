#include "../proto.h"
#include "../globals.h"
#include "effects.h"

/* Select effect function
 * Select particle weighting function
 * Select metric function
 * Select pix modification
 * Select required particles
 * Select required blocks
 * Put Name and Unit strings
 * Add an optional message
 * Do other special things before main loop */
void select_effect_module(int choice)
{
	int i;
	char message[MAXLINELENGTH] = " ",	// Printed after unit 
	    message2[MAXLINELENGTH] = " ";	// Helper 

	start_timing(CPU_EFFECT);

	/* standard values */
	for (i = 0; i < 1000; i++)	// required blocks from snap
		strcpy(Effect.Req.Block[i], "    ");

	int iEffect = 0;	// These blocks are always needed
	strcpy(Effect.Req.Block[iEffect++], "POS "); // label has to be 4 chars !!
	strcpy(Effect.Req.Block[iEffect++], "MASS");
	strcpy(Effect.Req.Block[iEffect++], "ID  ");

	for (i = 0; i < N_part_types; i++)	// required particle types
		Effect.Req.PartType[i] = false;

	Effect.Req.Tree = false;	// dont build the tree

	Effect.Nimage = 1;	// only one image

	/* pointers */
	effect_func_ptr = NULL;
	metric_func_ptr = &Metric_Hsml2;	// measure cost, h^2 is usually good 

	prep_func_ptr = NULL;	// do stuff before the loop 
	post_func_ptr = NULL;	// do stuff after the loop 

	switch (choice) {	// Choose effect

	case 0:
		
		switch (Param.Effect_Flag) {
		
		case 0:
			
			strcpy(Effect.Name, "2D Gas Density");	/* Name */
			strcpy(Effect.Descr[0], Effect.Name);
			strcpy(Effect.Unit, "[g/cm^2]");	/* Unit of image */
			
			Effect.Req.PartType[0] = true;	/* Particle types required */
			
			effect_func_ptr = &Density;	/* Point to effect function */
			weight_func_ptr = &Part_Weight_One;	/* No weighting */
		
			break;
		
		case 1:
			
			strcpy(Effect.Name, "3D Electron Number Density");
			strcpy(Effect.Descr[0], Effect.Name);
			strcpy(Effect.Unit, "[1/cm^2]");
			
			Effect.Req.PartType[0] = true;
			
			effect_func_ptr = &NumDensity;
			weight_func_ptr = &Part_Weight_One;
			
			break;
		}

		/* Implicitly: No Pixel pre & post modification. 1 image
		 * POS, RHO, HSML, ID needed*/
		break;

	case 1:

#ifndef VTURB
		Assert(0, "Code Compiled without VTURB");
#endif
		strcpy(Effect.Unit, "[cm/s]");

		Effect.Req.PartType[0] = true;

		switch (Param.Effect_Flag) {
		
		case 0:
		
			strcpy(Effect.Name, "Velocity");
			strcpy(Effect.Req.Block[iEffect++], "VEL ");
			
			effect_func_ptr = &Vel;
			
			break;
		
		case 1:
		
			strcpy(Effect.Name, "RMS velocity around central");
			strcpy(Effect.Req.Block[iEffect++], "VELT");
			
			effect_func_ptr = &Velt;
			
			break;

		case 2:
		
			strcpy(Effect.Name, "RMS velocity around mean");
			strcpy(Effect.Req.Block[iEffect++], "VRMS");
			
			effect_func_ptr = &Vrms;
			
			break;

		case 3:
		
			strcpy(Effect.Name, "Mean velocity in hsml");
			strcpy(Effect.Req.Block[iEffect++], "VBULK");
			strcpy(Effect.Req.Block[iEffect++], "TNGB");
			
			effect_func_ptr = &Vbulk;
			
			break;

		case 4:
		
			strcpy(Effect.Name, "Velocity Divergence in hsml");
			strcpy(Effect.Req.Block[iEffect++], "VDIV");
			
			effect_func_ptr = &Vdiv;
			
			break;

		case 5:
		
			strcpy(Effect.Name, "Velocity Rotation in hsml");
			strcpy(Effect.Req.Block[iEffect++], "VROT");
			
			effect_func_ptr = &Vrot;
			
			break;

		case 6:
		
			strcpy(Effect.Name, "Speed of Sound");
			strcpy(Effect.Req.Block[iEffect++], "VEL ");
		
			effect_func_ptr = &Vsound;
			
			break;

		case 7:
		
			strcpy(Effect.Name,
			       "Scaled RMS velocity around central");
			strcpy(Effect.Req.Block[iEffect++], "VELT");
			
			effect_func_ptr = &Velt_scaled;

			sprintf(message, "Scaling Velt to <%g> kpc \n",
				Param.TurbScale * Unit.Length / kpc2cgs);
			
			break;

		case 8:
		
			strcpy(Effect.Name, "Scaled RMS velocity around mean");
			strcpy(Effect.Req.Block[iEffect++], "VRMS");
			
			effect_func_ptr = &Vrms_scaled;

			sprintf(message, "Scaling Vrms to <%g> kpc \n",
				Param.TurbScale * Unit.Length / kpc2cgs);
			
			break;

		case 9:
		
			strcpy(Effect.Name,
			       "RMS tangential velocity around mean");
			strcpy(Effect.Req.Block[iEffect++], "VTAN");
			
			effect_func_ptr = &Vtan;
			
			break;

		case 10:
		
			strcpy(Effect.Name, "RMS radial velocity around mean");
			strcpy(Effect.Req.Block[iEffect++], "VRAD");
			
			effect_func_ptr = &Vrad;
			
			break;

		case 11:
		
			strcpy(Effect.Name, "Velocity perpendicular to LoS");
			strcpy(Effect.Req.Block[iEffect++], "VEL ");
			
			effect_func_ptr = &Vel_perp;
			
			break;

		case 12:
		
			strcpy(Effect.Name, "Vorticity");
			strcpy(Effect.Req.Block[iEffect++], "VORT");
			
			effect_func_ptr = &Vorticity;
			
			break;

		case 13:
		
			strcpy(Effect.Name, "Kernel Weighted Velocity");
			strcpy(Effect.Req.Block[iEffect++], "VEL ");
			
			effect_func_ptr = &Kernel_weighted_velocity;
			
			Effect.Req.Tree = true;
			
			break;

		case 14:
		
			strcpy(Effect.Name, "Kernel Weighted Bulk Velocity");
			strcpy(Effect.Req.Block[iEffect++], "VEL ");
			
			effect_func_ptr = &Kernel_weighted_bulk_velocity;
			
			Effect.Req.Tree = true;
			
			break;

		case 15:
		
			strcpy(Effect.Name,
			       "Kernel Weighted Velocity Dispersion");
			strcpy(Effect.Req.Block[iEffect++], "VEL ");
			
			effect_func_ptr = &Kernel_weighted_velocity_dispersion;
			
			Effect.Req.Tree = true;
			
			break;

		case 16:
		
			strcpy(Effect.Name, "Turbulent Energy");
			strcpy(Effect.Req.Block[iEffect++], "VEL ");

			effect_func_ptr = &Turbulent_energy;
			
			Effect.Req.Tree = true;
			
			break;
		} // switch Effect_Flag

		weight_func_ptr = &Part_Weight_Rho;

		strcpy(Effect.Descr[0], Effect.Name);
		
		break;

	case 2:

		strcpy(Effect.Name, "X-Ray Surface Brightness");
		strcpy(Effect.Descr[0], Effect.Name);
		strcpy(Effect.Unit, "[erg/cm^2/s/Hz]");

		sprintf(message, "E_min=%g keV\nE_max=%g keV\n",
			Param.E_min * 1e-3, Param.E_max * 1e-3);

		Effect.Req.PartType[0] = true;
		strcpy(Effect.Req.Block[iEffect++], "U   ");

		effect_func_ptr = &Xray;
		weight_func_ptr = &Part_Weight_Physical;

		break;

	case 3:
		
		switch (Param.Effect_Flag) {
		
		case 0:
			
			strcpy(Effect.Name, "Reacceleration Coefficient");
			strcpy(Effect.Unit, "[1/s]");

			strcpy(Effect.Descr[0], "Dpp/p^2");

			Effect.Req.PartType[0] = true;
			
			strcpy(Effect.Req.Block[iEffect++], "U   ");
			strcpy(Effect.Req.Block[iEffect++], "VRMS");
			strcpy(Effect.Req.Block[iEffect++], "BFLD");
			strcpy(Effect.Req.Block[iEffect++], "TNGB");

			set_Dpp_factors();

			effect_func_ptr = &Dpp;
			weight_func_ptr = &Part_Weight_Rho;
			metric_func_ptr = &Metric_One;

			break;

		case 1:
		
			strcpy(Effect.Name, "Reacceleration Timescale");
			strcpy(Effect.Unit, "[s]");

			strcpy(Effect.Descr[0], "tau_reacc");

			Effect.Req.PartType[0] = true;
			
			strcpy(Effect.Req.Block[iEffect++], "U   ");
			strcpy(Effect.Req.Block[iEffect++], "VRMS");
			strcpy(Effect.Req.Block[iEffect++], "BFLD");
			strcpy(Effect.Req.Block[iEffect++], "TNGB");

			set_Dpp_factors();

			effect_func_ptr = &Tau_reacc;
			weight_func_ptr = &Part_Weight_Rho;
			metric_func_ptr = &Metric_One;
			
			break;

		case 2:
		
			strcpy(Effect.Name, "Cooling Timescale");
			strcpy(Effect.Unit, "[s]");

			strcpy(Effect.Descr[0], "tau_rad");

			Effect.Req.PartType[0] = true;
			
			strcpy(Effect.Req.Block[iEffect++], "BFLD");

			set_Dpp_factors();

			effect_func_ptr = &Tau_cooling;
			weight_func_ptr = &Part_Weight_Rho;
			metric_func_ptr = &Metric_One;

			sprintf(message, "Frequency = %g GHz\n",
				Param.Freq * 1e-9);

			break;
		} // switch Effect_Flag

		break;

	case 4:

		switch (Param.Effect_Flag) {
		
		case 0:
		
			strcpy(Effect.Name, "Mass Weighted Temperature");
			
			weight_func_ptr = &Part_Weight_Rho;
			
			strcpy(Effect.Unit, "[K]");
			
			effect_func_ptr = &Temp;
			
			break;

		case 1:
			
			strcpy(Effect.Name, "Sound speed");
			
			weight_func_ptr = &Part_Weight_Rho;
			
			effect_func_ptr = &Soundspeed;
			
			strcpy(Effect.Unit, "[cm/s]");
			
			break;

		case 2:	/* Borgani 01 */
		
			strcpy(Effect.Name, "Emission Weighted Temperature");
			
			weight_func_ptr = &Part_Weight_Emission;
			
			strcpy(Effect.Unit, "[K]");
			
			effect_func_ptr = &Temp;
			
			break;

		case 3:	/* Mazzotta 04 */
		
			strcpy(Effect.Name, "Spectroscopic Temperature");
			
			weight_func_ptr = &Part_Weight_Spectroscopic;
			
			strcpy(Effect.Unit, "[K]");
			
			effect_func_ptr = &Temp;
			
			break;

		case 4:	/* like smac1 */
		
			strcpy(Effect.Name, "Xray Band EW Temperature");
			
			weight_func_ptr = &Part_Weight_XrayBand;
			
			strcpy(Effect.Unit, "[K]");
			
			effect_func_ptr = &Temp;
			
			break;
		} // switch Effect_Flag

		strcpy(Effect.Descr[0], Effect.Name);
		
		Effect.Req.PartType[0] = true;
		
		strcpy(Effect.Req.Block[iEffect++], "U   ");

		break;

	case 5:
		
		switch (Param.Effect_Flag) {
		
		case 0:
		
			strcpy(Effect.Name, "Thermal Pressure");
			strcpy(Effect.Descr[0], Effect.Name);
			strcpy(Effect.Unit, "[g/s^2/cm]");

			Effect.Req.PartType[0] = true;

			effect_func_ptr = &Therm_Pressure;
			weight_func_ptr = &Part_Weight_Rho;
			
			break;

		case 1:
		
			strcpy(Effect.Name, "Magnetic Pressure");
			strcpy(Effect.Descr[0], Effect.Name);
			strcpy(Effect.Unit, "[g/s^2/cm]");

			Effect.Req.PartType[0] = true;

			strcpy(Effect.Req.Block[iEffect++], "BFLD");

			effect_func_ptr = &Magn_Pressure;
			weight_func_ptr = &Part_Weight_Rho;
		
			break;

		case 2:
		
			strcpy(Effect.Name, "Turbulent Pressure");
			strcpy(Effect.Descr[0], Effect.Name);
			strcpy(Effect.Unit, "[g/s^2/cm]");

			Effect.Req.PartType[0] = true;

			strcpy(Effect.Req.Block[iEffect++], "VRMS");

			effect_func_ptr = &Turb_Pressure;
			weight_func_ptr = &Part_Weight_Rho;
			
			break;
		} // switch Effect_Flag

		break;

	case 6:
		
		switch (Param.Effect_Flag) {
		
		case 0:
		
			strcpy(Effect.Name, "Total Magnetic Field Strength");
			strcpy(Effect.Unit, "[G]");
			
			effect_func_ptr = &bfld;
			
			weight_func_ptr = &Part_Weight_Rho;
			
			break;

		case 1:
		
			strcpy(Effect.Name,"Orthogonal Magnetic Field Strength");
			strcpy(Effect.Unit, "[G]");
			
			effect_func_ptr = &bfld_orth;
			weight_func_ptr = &Part_Weight_Rho;
			
			break;

		case 2:
		
			strcpy(Effect.Name, "Faraday Rotation Measurement (RM)");
			strcpy(Effect.Unit, "[rad/cm^2]");
			sprintf(message, "\nUNITS CHANGED BY 1E4 FROM SMAC I\n");

			set_faraday_factors();
			
			effect_func_ptr = &faraday_rotation;
			weight_func_ptr = &Part_Weight_Physical;
			
			break;

		case 3:
		
			strcpy(Effect.Name, "Alven Velocity");
			strcpy(Effect.Unit, "[km/s]");
			
			effect_func_ptr = &alven_vel;
			
			weight_func_ptr = &Part_Weight_Rho;
			
			strcpy(Effect.Descr[0], "v_a = sqrt(B^2/rho /4/pi)");
			
			break;

		case 4:
		
			strcpy(Effect.Name, "Div(B)*Hsml/B");
			strcpy(Effect.Unit, "[1]");
			
			effect_func_ptr = &divB;
			
			weight_func_ptr = &Part_Weight_Rho;
			
			strcpy(Effect.Req.Block[iEffect++], "DIVB");
			strcpy(Effect.Descr[0], "");
			
			break;

		case 5:
		
			strcpy(Effect.Name, "Plasma Beta");
			strcpy(Effect.Unit, "[1]");
			
			effect_func_ptr = &plasma_beta;
			weight_func_ptr = &Part_Weight_Rho;
			
			strcpy(Effect.Req.Block[iEffect++], "BFLD");
			strcpy(Effect.Descr[0], "");
			
			break;

		default:
		
			Assert(0, "This Subflag doesn't exist");
			
			break;
		} // switch Effect_Flag


		strcpy(Effect.Descr[0], Effect.Name);
		
		Effect.Req.PartType[0] = true;
		
		strcpy(Effect.Req.Block[iEffect++], "BFLD");

		break;

	case 7:

		Effect.Req.PartType[0] = true;
		
		strcpy(Effect.Req.Block[iEffect++], "U   ");
		
		switch (Param.Effect_Flag) {
		
		case 0:
		
			strcpy(Effect.Name, "Compton-y Parameter");
			
			effect_func_ptr = &comptonY;
			
			set_sz_factors(false);
			
			break;

		case 1:
		
			strcpy(Effect.Name, "Thermal Sunyaev Zeldovich Effect (DT/T)");

			effect_func_ptr = &tSz;
			
			set_sz_factors(false);
			
			break;

		case 2:
		
			strcpy(Effect.Name,"Thermal Sunyaev Zeldovich Effect (DI/I)");

			effect_func_ptr = &tSz;
			
			set_sz_factors(true);
			
			break;

		case 3:
		
			strcpy(Effect.Name, "Kinetic Sunyaev Zeldovich Effect (DT/T)");

			effect_func_ptr = &kSz;
			
			set_sz_factors(false);
			
			break;

		case 4:
			
			strcpy(Effect.Name, "Kinetic Sunyaev Zeldovich Effect (DI/I)");
	
			effect_func_ptr = &kSz;
			
			set_sz_factors(true);
			
			break;
		} // switch Effect_Flag

		strcpy(Effect.Unit, "none");
		strcpy(Effect.Descr[0], Effect.Name);
		
		sprintf(message, "Frequency = %g Ghz => x = %g\n\n"
			"DEFINITIONS DIFFERENT FROM SMAC I\n",
			Param.Freq / 1e9, h_planck * Param.Freq / (k_B * Tcmb));
		
		weight_func_ptr = &Part_Weight_Physical;
		
		break;

	case 8:

		switch (Param.Effect_Flag) {
		
		case 0:
	
			strcpy(Effect.Name, "Differential Gamma Rays from CRs");
			strcpy(Effect.Unit, "[1/cm^2/erg]");
			
			effect_func_ptr = &crGamma_dif;

			sprintf(message, "@ E_min = %g GeV\nX_cr = %g\na_crp = %g\n",
				Param.E_min * 1e-9, Param.X_crp, Param.a_crp);
			
			break;
		
		case 1:
			
			strcpy(Effect.Name, "Integral Gamma Rays from CRs");
			strcpy(Effect.Unit, "[1/cm^2]");
			
			effect_func_ptr = &crGamma_int;
			
			sprintf(message, "E_min = %e GeV\n"
				"E_max = %e GeV\n"
				"X_cr  = %g\n"
				"a_crp = %g\n",
				Param.E_min * 1e-9, Param.E_max * 1e-9,
				Param.X_crp, Param.a_crp);
			
			break;
		} // switch Effect_Flag

		weight_func_ptr = &Part_Weight_Physical;
		metric_func_ptr = &Metric_One;

#ifdef XCRP_SCALING
		prep_func_ptr = &find_central_density;
#endif
		set_crGamma_factors();

		strcpy(Effect.Descr[0], Effect.Name);

		Effect.Req.PartType[0] = true;
		
		break;

	case 9:
		
		strcpy(Effect.Unit, "[erg*s/cm^2]");

		weight_func_ptr = &Part_Weight_Physical;
		effect_func_ptr = &synchrotron;
		metric_func_ptr = &Metric_One;	// particle dominated 

		set_synchro_factors(Param.Effect_Flag);	//include prep_func_ptr here 

		Effect.Req.PartType[0] = true;
		strcpy(Effect.Req.Block[iEffect++], "U   ");
		strcpy(Effect.Req.Block[iEffect++], "BFLD");

		strcpy(Effect.Name, "Synchrotron Total Emission");
		strcpy(Effect.Descr[0], "Total Intensity");

#ifdef POLARISATION
		strcpy(Effect.Name, "Synchrotron Emission Total & "
		       "Polarisation");
		sprintf(Effect.Descr[0], "Stokes I");
		sprintf(Effect.Descr[1], "Stokes Q");
		sprintf(Effect.Descr[2], "Stokes U");
#ifdef ADD_IPOL_CHI_PI
		sprintf(Effect.Descr[3], "I_polarised");
		sprintf(Effect.Descr[4], "Polarisation Angle");
		sprintf(Effect.Descr[5], "Degree of polarisation");
		post_func_ptr = &add_Ipol_chi_Pi;
#endif
		Effect.Nimage = 3;	/* I,Q,U  */

#ifdef INTRINSIC_RM
		set_intrinsic_RM();
#endif
#endif				/* POLARISATION */

		sprintf(message,
			"Freq  = %1.4f GHz\n"
			"E_min = %1.4g GeV\n"
			"X_cr  = %1.4g thermal\n"
			"Pitch Angle Integration : %d \n"
			"%s\n"
			"UNITS CHANGED BY 1E7 (W -> erg) FROM SMAC I\n",
			Param.Freq * 1e-9, Param.E_min * 1e-9,
			Param.X_crp, Param.SynchroPAngInt, message2);

		break;

	case 11:
		
		strcpy(Effect.Name, "DM Annihilation");
		strcpy(Effect.Unit, "[erg*s/cm^2]");
		strcpy(Effect.Descr[0], Effect.Name);

		effect_func_ptr = &dm_annihilation;
		weight_func_ptr = &Part_Weight_Physical;
		metric_func_ptr = &Metric_Hsml2;

		Effect.Req.PartType[1] = true;

		Effect.Req.Tree = true;

		break;

	case 10:

		strcpy(Effect.Name, "Dark Matter Density ");
		strcpy(Effect.Unit, "[cm^-2]");
		strcpy(Effect.Descr[0], Effect.Name);

		effect_func_ptr = &dm_density;
		weight_func_ptr = &Part_Weight_Physical;
		metric_func_ptr = &Metric_Hsml2;

		Effect.Req.PartType[1] = true;

		Effect.Req.Tree = true;

		break;

	case 12:

		strcpy(Effect.Name, "Synchrotron Emission from powerlaw CRe, "
		       "Analytic formula \n");
		strcpy(Effect.Unit, "[erg*s/cm^2]");

		sprintf(Effect.Descr[0], "Stokes I");
		sprintf(Effect.Descr[1], "Stokes Q");
		sprintf(Effect.Descr[2], "Stokes U");
		sprintf(Effect.Descr[3], "I_polarised");
		sprintf(Effect.Descr[4], "Polarisation Angle");
		sprintf(Effect.Descr[5], "Degree of polarisation");

		Effect.Nimage = 3;	/* I,Q,U  */

		set_aSynchro_factors();

		effect_func_ptr = &aSynchrotron;
		weight_func_ptr = &Part_Weight_Physical;
		metric_func_ptr = &Metric_Hsml2;

#ifdef ADD_IPOL_CHI_PI
		post_func_ptr = &add_Ipol_chi_Pi;
#endif
#ifdef INTRINSIC_RM
		set_intrinsic_RM();
#endif

		Effect.Req.PartType[0] = true;
	
		strcpy(Effect.Req.Block[iEffect++], "U   ");
		strcpy(Effect.Req.Block[iEffect++], "BFLD");
		strcpy(Effect.Req.Block[iEffect++], "MACH");
	
		sprintf(message, "CR Spectral Index s = %1.4f \nX_CRp = %g\n", 
				Param.a_crp, Param.X_crp);

		break;

	case 13:

#ifndef RADTRANSFER
		Assert(0, "This effect requires RADTRANSFER");
#endif
		
		strcpy(Effect.Name, "Ionized Hydrogen Fraction");
		strcpy(Effect.Unit, "[1/cm^3]");
		strcpy(Effect.Descr[0], Effect.Name);

		effect_func_ptr = &rt_Hydrogen;
		weight_func_ptr = &Part_Weight_Rho;
		metric_func_ptr = &Metric_Hsml2;

		Effect.Req.PartType[0] = true;
		
		strcpy(Effect.Req.Block[iEffect++], "nHII");
	
		break;

	case 14:

		switch (Param.Effect_Flag) {
		
		case 0:
		
			strcpy(Effect.Name,
			       "Mean Free Path from Coulomb Scattering");
			strcpy(Effect.Unit, "[cm]");
			strcpy(Effect.Descr[0], Effect.Name);

			effect_func_ptr = &coul_mfp;
			weight_func_ptr = &Part_Weight_Rho;
			metric_func_ptr = &Metric_Hsml2;

			Effect.Req.PartType[0] = true;
		
			break;

		case 1:
			
			strcpy(Effect.Name, "Collisionality from Coulomb Scattering");
			strcpy(Effect.Unit, "n/a");
			strcpy(Effect.Descr[0], Effect.Name);

			effect_func_ptr = &coul_coll;
			weight_func_ptr = &Part_Weight_Rho;
			metric_func_ptr = &Metric_Hsml2;

			Effect.Req.PartType[0] = true;

			break;
		} // switch Effect_Flag

		break;

	case 15:
		
		switch (Param.Effect_Flag) {
		
		case 0:
		
			strcpy(Effect.Name, "Mach Number (Beck+ 2011)");
			strcpy(Effect.Unit, " n/a ");
			strcpy(Effect.Descr[0], Effect.Name);

			effect_func_ptr = &Machnumber;
			weight_func_ptr = &Part_Weight_Rho;
			metric_func_ptr = &Metric_Hsml2;

			Effect.Req.PartType[0] = true;
			strcpy(Effect.Req.Block[iEffect++], "MACH");

#ifndef SHOCKS
			Assert(0, "This effect requires SHOCKS");
#endif
			break;

		case 1:

			strcpy(Effect.Name,
			       "Time Dependent Artificial Viscosity Parameter");
			strcpy(Effect.Unit, " n/a ");
			strcpy(Effect.Descr[0], Effect.Name);

			effect_func_ptr = &AlphaVisc;
			weight_func_ptr = &Part_Weight_Rho;
			metric_func_ptr = &Metric_Hsml2;

			Effect.Req.PartType[0] = true;
			strcpy(Effect.Req.Block[iEffect++], "ABVC");

#ifndef TIME_DEP_ART_VISC
			Assert(0, "This effect requires TIME_DEP_ART_VISC");
#endif
			break;

		} // switch Effect_Flag

		break;

	default:
			Assert(0, "Effect Flags %d %d not handled", 
					choice ,Param.Effect_Flag);

			break;

	} // switch choice
	
	if (! Effect.Req.Tree) {

		strcpy(Effect.Req.Block[iEffect++], "RHO ");
		strcpy(Effect.Req.Block[iEffect++], "HSML");
	
	}

	MPI_Barrier(MPI_COMM_WORLD);

	Assert(effect_func_ptr != NULL 
			&& weight_func_ptr != NULL
	       	&& metric_func_ptr != NULL,  "Effect or Weight or Metric not set");

	rprintf("Projecting %s \nUnit %s\n\n%s\nUsing <%d> image(s)\n",
		Effect.Name, Effect.Unit, message, Effect.Nimage);

	fflush(stdout);

	stop_timing(CPU_EFFECT);

	return;
}

/* Following Functions are used as
 * particle weights, where the first
 * corresponds to no weighting
 * */
extern double Part_Weight_One(int ipart)
{
	return 1;
}

extern double Part_Weight_Physical(int ipart)
{
	return pix2cm;		/* We integrate LoS in Pixels, make cm */
}

extern double Part_Weight_Rho(int ipart)
{
	return P[ipart].Rho;
}

extern double Part_Weight_Emission(int ipart)
{
	return p2(P[ipart].Rho) * sqrt(Temperature(ipart));
}

extern double Part_Weight_XrayBand(int ipart)
{
	float T_eV = Temperature(ipart) / eV2cgs;

	return (exp(-Param.E_min / (k_B * T_eV))
		- exp(-Param.E_max / (k_B * T_eV)));
}

/*  Mazotta+ 04*/
extern double Part_Weight_Spectroscopic(int ipart)
{
	return p2(P[ipart].Rho) * pow(Temperature(ipart), 0.75 - 1.5);
}
