/* Common input functions */
#include "../globals.h"

#define LASTPARAMETERID -1111

#define REAL 1
#define STRING 2
#define INT 3

long guess_snapnum(char *fname);

void select_cosmology(int cosmology)
{
	switch (cosmology) {
	case 0:
		Cosmo = null;
		break;
	case 1:
		Cosmo = WMAP3;
		break;
	case 2:
		Cosmo = SMAC1;
		break;
	default:
		Assert(0, "Invalid Cosmology in Parameter File");
		break;
	}

	rprintf("Setting Cosmology : %s \n\n", Cosmo.name);

	return;
}

/* Reads a number of tags from an ascii file
 * the comment sign is  %
 * */
void read_param_file(char *filename)
{

	FILE *fd;
	char buf[MAXLINELENGTH]; 
	int tagDone[MAXTAGS] = { 0 };
	int i, nt = 0;

	Assert(sizeof(long long) == 8,
	       "Type `long long' is not 64 bit on this platform.");
	Assert(sizeof(int) == 4, "Type `int' is not 32 bit on this platform");
	Assert(sizeof(float) == 4,
	       "Type `float' is not 32 bit on this platform");
	Assert(sizeof(double) == 8,
	       "Type `double' is not 64 bit on this platform");

	if (Task.Rank == 0) {	// Task 0 only

		strcpy(tag[nt], "Center_X");
		strcpy(comment[nt], "Center of image - x component");
		addr[nt] = &Param.Center[0];
		id[nt++] = REAL;

		strcpy(tag[nt], "Center_Y");
		strcpy(comment[nt], "Center of image - y component");
		addr[nt] = &Param.Center[1];
		id[nt++] = REAL;

		strcpy(tag[nt], "Center_Z");
		strcpy(comment[nt], "Center of image - z component");
		addr[nt] = &Param.Center[2];
		id[nt++] = REAL;

		strcpy(tag[nt], "Use_Barycenter");
		strcpy(comment[nt], "Center Image on Center of Mass");
		addr[nt] = &Param.Flag_Barycenter;
		id[nt++] = INT;

		strcpy(tag[nt], "XYSize");
		strcpy(comment[nt], "[kpc] Comoving FoV");
		addr[nt] = &Param.XYSize;
		id[nt++] = REAL;

		strcpy(tag[nt], "ZDepth");
		strcpy(comment[nt], "[kpc] Comoving Column Depth");
		addr[nt] = &Param.ZDepth;
		id[nt++] = REAL;

		strcpy(tag[nt], "XYPix");
		strcpy(comment[nt], "Number of pixels a side");
		addr[nt] = &Param.XYPix;
		id[nt++] = INT;

		strcpy(tag[nt], "Euler_Angle_0");
		strcpy(comment[nt], "Euler Angle phi");
		addr[nt] = &Param.EulAng[0];
		id[nt++] = REAL;

		strcpy(tag[nt], "Euler_Angle_1");
		strcpy(comment[nt], "Euler Angle theta");
		addr[nt] = &Param.EulAng[1];
		id[nt++] = REAL;

		strcpy(tag[nt], "Euler_Angle_2");
		strcpy(comment[nt], "Euler Angle psi");
		addr[nt] = &Param.EulAng[2];
		id[nt++] = REAL;

#ifdef HEALPIX
		strcpy(tag[nt], "NSide");
		strcpy(comment[nt],
		       "HEALPIX Resolution Parameter Npix=12*Nside^2");
		addr[nt] = &Param.NSide;
		id[nt++] = INT;

		strcpy(tag[nt], "Rmin");
		strcpy(comment[nt], "HEALPIX minimal distance from observer");
		addr[nt] = &Param.Rmin;
		id[nt++] = REAL;

		strcpy(tag[nt], "Rmax");
		strcpy(comment[nt], "HEALPIX maximal distance from observer");
		addr[nt] = &Param.Rmax;
		id[nt++] = REAL;
#endif

		strcpy(tag[nt], "Output_File");
		strcpy(comment[nt], "Output File Name ");
		addr[nt] = &Param.Output_File;
		id[nt++] = STRING;

		strcpy(tag[nt], "Input_File");
		strcpy(comment[nt], "Input File Name");
		addr[nt] = &Param.Input_File;
		id[nt++] = STRING;

		strcpy(tag[nt], "N_IOTasks");
		strcpy(comment[nt], "Number of files read in parallel");
		addr[nt] = &Param.N_IOTasks;
		id[nt++] = INT;

		strcpy(tag[nt], "Effect_Module");
		strcpy(comment[nt], "Effect Module Number");
		addr[nt] = &Param.Effect_Module;
		id[nt++] = INT;

		strcpy(tag[nt], "Effect_Flag");
		strcpy(comment[nt], "Effect SubFlag");
		addr[nt] = &Param.Effect_Flag;
		id[nt++] = INT;

		strcpy(tag[nt], "Cosmology");
		strcpy(comment[nt], "Cosmology used");
		addr[nt] = &Param.Cosmology;
		id[nt++] = INT;

		strcpy(tag[nt], "NoClobber");
		strcpy(comment[nt], "Overwrite Image");
		addr[nt] = &Param.NoClobber;
		id[nt++] = INT;

		strcpy(tag[nt], "E_min");
		strcpy(comment[nt], "[eV] Energy Range Minimum");
		addr[nt] = &Param.E_min;
		id[nt++] = REAL;

		strcpy(tag[nt], "E_max");
		strcpy(comment[nt], "[eV] Energy Range Maximum");
		addr[nt] = &Param.E_max;
		id[nt++] = REAL;

		strcpy(tag[nt], "Freq");
		strcpy(comment[nt], "[Hz] Frequency of Observation");
		addr[nt] = &Param.Freq;
		id[nt++] = REAL;

		strcpy(tag[nt], "a_cr");
		strcpy(comment[nt], "Cosmic Ray Spectral Index");
		addr[nt] = &Param.a_crp;
		id[nt++] = REAL;

		strcpy(tag[nt], "X_cr");
		strcpy(comment[nt], "CRp Normalisation rel. to thermal");
		addr[nt] = &Param.X_crp;
		id[nt++] = REAL;

		strcpy(tag[nt], "PitchAngInt");
		strcpy(comment[nt], "Toggle Pitch Angle Integration");
		addr[nt] = &Param.SynchroPAngInt;
		id[nt++] = INT;

		strcpy(tag[nt], "IntrRM");
		strcpy(comment[nt], "Intrinsic RM");
		addr[nt] = &Param.SynchroIntrRM;
		id[nt++] = INT;

		strcpy(tag[nt], "t_turb");
		strcpy(comment[nt], "Injection timescale of turbulence");
		addr[nt] = &Param.TurbTime;
		id[nt++] = REAL;

		strcpy(tag[nt], "eta_t");
		strcpy(comment[nt], "Magnetosonic fraction");
		addr[nt] = &Param.Eta_t;
		id[nt++] = REAL;

		strcpy(tag[nt], "turb_scale");
		strcpy(comment[nt], "Scale of turbulent velocity");
		addr[nt] = &Param.TurbScale;
		id[nt++] = REAL;

		strcpy(tag[nt], "UnitLength_in_cm");
		strcpy(comment[nt], "[cm] Unit Length");
		addr[nt] = &Unit.Length;
		id[nt++] = REAL;

		strcpy(tag[nt], "UnitMass_in_g");
		strcpy(comment[nt], "[g] Unit Mass");
		addr[nt] = &Unit.Mass;
		id[nt++] = REAL;

		strcpy(tag[nt], "UnitVelocity_in_cm_per_s");
		strcpy(comment[nt], "[cm/s] Unit Vel");
		addr[nt] = &Unit.Vel;
		id[nt++] = REAL;

		id[nt] = LASTPARAMETERID;

		Param.NCube = 1;
			
		if ((fd = fopen(filename, "r"))) {

			sprintf(buf, "%s%s", filename, "-usedvalues");

			printf("\nReading Parameter file : %s \n", filename);

			while (fgets(buf, MAXLINELENGTH, fd)) {

				char buf1[MAXLINELENGTH] = { "" },buf2[MAXLINELENGTH] = { "" },
					 buf3[MAXLINELENGTH] = { "" },buf4[MAXLINELENGTH] = { "" }, 
					 buf5[MAXLINELENGTH] = { "" };

				int n = sscanf(buf, "%s%s%s%s%s", buf1,buf2,buf3,buf4,buf5);

				if ( n < 2)
					continue;

				if (buf1[0] == '%')
					continue;
	
				bool flag_cube = false;

				if ( ! (buf3[0] == '%' || n == 2) ) {  // make a cube

					flag_cube = true;
				
					printf("\nIterating '%s' from <%s> to <%s>,"
						" step <%s>\n",	buf1, buf2, buf3, buf4);
				}

				int j = -1;

				for (i = 0; i < nt; i++) {

					if ((strcmp(buf1, tag[i]) == 0) && (tagDone[i] != 1)) {

						j = i;
						
						tagDone[i] = 1;
						
						break;
					}
				}

				if (j >= 0) {

					switch (id[j]) {

					case REAL:
						
						*((double *)addr[j]) = atof(buf2);

						if (flag_cube) {

							Param.CubeType = 1;

							Param.CubePtr = addr[j];

							Param.CubeMin = atof(buf2);
							Param.CubeMax = atof(buf3);
							Param.CubeStep = atof(buf4);

							Param.NCube = (Param.CubeMax - Param.CubeMin) 
								/ Param.CubeStep + 1;
						}

						break;
					
					case STRING:

						strcpy((char *)addr[j], buf2);
	
						if (flag_cube) {

							Param.CubeType = 2;
							
							Param.CubePtr = addr[j];

							Param.CubeMin = guess_snapnum(buf2);
							Param.CubeMax = guess_snapnum(buf3);
							Param.CubeStep = atof(buf4);

							Param.NCube = (Param.CubeMax - Param.CubeMin) 
								/ Param.CubeStep + 1;
						}

						break;
					
					case INT:
						
						*((int *)addr[j]) = atoi(buf2);
						
						if (flag_cube) {

							Param.CubeType = 3;
							
							Param.CubePtr = addr[j];

							Param.CubeMin = atof(buf2);
							Param.CubeMax = atof(buf3);
							Param.CubeStep = atof(buf4);

							Param.NCube = (Param.CubeMax - Param.CubeMin) 
								/ Param.CubeStep  + 1;
						}

						break;
					}
				}
			}
		
			fclose(fd);

			printf("\n");
		
		} else
			Assert(0, "Parameter file not found %s", filename);

		bool good = true;

		for (i = 0; i < nt; i++) {

			if (!tagDone[i]) {

				fprintf(stderr,
					"Value for tag '%s' missing in parameter file '%s'.\n",
					tag[i], filename);

				good = false;
			}
		}

		Assert(good, "Parameter file incomplete");
	}

	MPI_Bcast(&Param, sizeof(Param), MPI_BYTE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&Unit, sizeof(Unit), MPI_BYTE, 0, MPI_COMM_WORLD);

	Assert(Param.E_min < Param.E_max, "Parameter Energy range: Emin >= Emax");

	//printf("%d %d %g %g %g \n", Param.CubeType, Param.NCube, Param.CubeMax, Param.CubeMin, Param.CubeStep);

	return;
}

#undef REAL
#undef STRING
#undef INT

long guess_snapnum(char *fname)
{
	char *token, *last_token = NULL, file[MAXLINELENGTH];

	strcpy(file, fname);	// protect fname from strtok 

	token = strtok(file, "_");
	
	while (token != NULL) {

		last_token = token;

		token = strtok(NULL, "_");
	}

	long snapnum = atoi(last_token);

	//rprintf("Guessing snapshot number as : %ld \n\n", snapnum);

	return (snapnum);
}
