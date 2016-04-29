/* This is P-Smac II
 * 
 * Line-Of-Sight integration of simulations in parallel using different 
 * complex emission models (Donnert et al. 2012).
 *
 * Copyright (C) 2010 Julius Donnert: donnert AT ira.inaf.it

 * This program is free software; you can redistribute it and/or modify it 
 * under the terms of the GNU General Public License as published by the 
 * Free Software Foundation; either version 3 of the License, or (at your
 * option) any later version.

 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
 * more details.

 * You should have received a copy of the GNU General Public License along 
 * with this program; if not, see <http://www.gnu.org/licenses/>.
 * */

#include "proto.h"
#include "globals.h"
#include "tree.h"
#include "effects/effects.h"

static void preamble(int, char **);

#ifdef O_FITS_COMPRESSED_CUBE
void advance_parameter(const int i);
#endif

int main(int argc, char *argv[])
{
	preamble(argc, argv);

	init_timing();

	read_param_file(argv[1]);

#ifdef O_FITS_COMPRESSED_CUBE
	for (int i = 0; i < Param.NCube; i++) {

		advance_parameter(i);
#endif
	
	set_units();

	select_cosmology(Param.Cosmology);

	select_effect_module(Param.Effect_Module);

	read_snapshot(Param.Input_File);
	
	Remove_Unused_Particles();

	if (Effect.Req.Tree)
  		Find_SPH_Densities();

	setup();

	domain_decomposition();

	project();

#ifdef O_FITS_COMPRESSED_CUBE
		move_image_to_cube(i);
	}
#endif

	write_output();

	finish_timing();

	MPI_Finalize();

	return EXIT_SUCCESS;	// Finish him of - Fatality 
}

/* technicalities */
static void preamble(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &Task.Rank);
	MPI_Comm_size(MPI_COMM_WORLD, &Task.NTask);

#pragma omp parallel
	{
		Omp.ThreadID = omp_get_thread_num();
		Omp.NThreads = omp_get_num_threads();
	}

	if (Task.Rank == 0) {
		printf("*** This is P-Smac II Version %s ***\n\n", VERSION);

		print_compile_time_settings();

		Assert(argc == 2,
		       "Parameter file missing or too many parameters!");

		printf("\nRunning on <%d> nodes", Task.NTask);

#pragma omp parallel
		if (Omp.ThreadID == 0)
			printf(" and <%d> threads per node", Omp.NThreads);

		printf("\n");

#ifdef PERIODIC
		printf("\n WARNING: You have compiled with PERIODIC.\n"
		       "This is nearly always useless and eats CPU time !\n\n");
#endif
	}

	return;
}

#ifdef O_FITS_COMPRESSED_CUBE
void advance_parameter(const int i) 
{

	char string[MAXLINELENGTH];
	char *token, *last_token, base[MAXLINELENGTH];
	if (i == 0)
		return;
	
	rprintf("New Iteration: %d\n", i);

	switch (Param.CubeType) {
	
	case 1: // REAL

		*((double *)Param.CubePtr) += Param.CubeStep;

		break;
	
	case 2: // STRING

		strncpy(string, Param.CubePtr, MAXLINELENGTH);

		token = strtok(string, "_");

		sprintf(base, "%s", token);
			
		token = strtok(NULL, "_");

		while (token != NULL) {

			last_token = token;
			
			token = strtok(NULL, "_");

			if (token != NULL)
				sprintf(base, "%s_%s", base, last_token);
		}

		sprintf((char *)Param.CubePtr, "%s_%03d",base, 
				atoi(last_token)+(int) Param.CubeStep);

		break;

	case 3: // INT
		
		*((int *)Param.CubePtr) += Param.CubeStep;

		break;
	}
	
	return ;
}
#endif
