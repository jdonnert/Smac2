/* Read Gadget data format 2 */

#include "../proto.h"
#include "../globals.h"
#include "gadget.h"

#if defined I_GADGET2 || defined I_AREPO

#define int4bytes int

#define READF77HEAD  {if(my_fread(&f77block,sizeof(int),1,fd)){swap_Nbyte((char*)&f77block,1,4);}}

#define UINT sizeof(unsigned int)
#define INT sizeof(long)
#define FLOAT sizeof(float)
#define DBL sizeof(double)

int4bytes f77block, Swap = 0;

int block_required(char *);

/* Main Reading :
 * This routine tries first to put one snapshot per task
 * Then N_IOTasks Master distribute to their groups
 * All files left are then read at the same time.
 * */
extern void read_snapshot(char *filename)
{
	long rest_files, groupLast,
	     group, file_number, nFiles;
	char buf[MAXLINELENGTH];

	start_timing(CPU_READIN);

	Snap.Have_Arepo = 0;

	Snap.SnapNum = guess_snapnum(Param.Input_File);

	/* Start parallel reading */
	nFiles = rest_files = find_files(filename);

	if (Param.N_IOTasks > nFiles) 
		Param.N_IOTasks = nFiles;

	int groupSize = Task.NTask / Param.N_IOTasks;
	
	if (Task.NTask % Param.N_IOTasks) 
		groupSize++;
	
	int groupMaster = (Task.Rank / groupSize) * groupSize;

	groupLast = groupMaster + groupSize - 1;

	if (groupLast > Task.NTask - 1)
		groupLast = Task.NTask - 1;

	while (rest_files > 0) {	/* Any left ? */

		if (nFiles == 1) { // read 1 to all 

			sprintf(buf, "%s", filename);

			read_file(buf, 0, Task.NTask - 1);

			MPI_Barrier(MPI_COMM_WORLD);

			rest_files -= 1;
			
		} else if (rest_files >= Task.NTask) { // Read in nIO blocks 

			file_number = Task.Rank + (rest_files - Task.NTask);

			sprintf(buf, "%s.%li", filename, file_number);

			for (group = 0; group < groupSize; group++) {
				if (Task.Rank == (groupMaster + group))
					read_file(buf, Task.Rank, Task.Rank);

				MPI_Barrier(MPI_COMM_WORLD);
			}

			rest_files -= Task.NTask;

			
		} else if (rest_files >= Param.N_IOTasks) { // Read make nIO groups 

			file_number = groupMaster / groupSize 
				+ (rest_files - Param.N_IOTasks);

			sprintf(buf, "%s.%li", filename, file_number);

			read_file(buf, groupMaster, groupLast);

			MPI_Barrier(MPI_COMM_WORLD);

			rest_files -= Param.N_IOTasks;
			
		} else { // reduce nIO to rest files 

			groupSize = Task.NTask / rest_files;
			
			if (Task.NTask % Param.N_IOTasks)
				groupSize++;

			groupMaster = (Task.Rank / groupSize) * groupSize;

			file_number = groupMaster / groupSize;

			sprintf(buf, "%s.%li", filename, file_number);

			read_file(buf, groupMaster, groupLast);

			MPI_Barrier(MPI_COMM_WORLD);

			rest_files -= groupSize;
		}
	}

	rprintf("\nRedshift of Snapshot : %g \n"
		"Boxsize : %g \n", Snap.Redshift, Snap.Boxsize);

	if (Snap.Have_Arepo)
		rprintf("Read an AREPO Snapshot ! \n\n");

	stop_timing(CPU_READIN);

	return;
}

/* Reads and Distributes a file 
 * */
void read_file(char *filename, int ReadTask, int LastTask)
{
	long i = 0, j = 0, task = 0;
	long target, src, nTask;
	int nRead[N_part_types] = { 0 }, nReadTot = 0;
	int nSend[N_part_types] = { 0 };
	unsigned char *comm_buf = NULL;
	size_t nBytes, byteOffset, partOffset, bufOffset;
	FILE *fd = NULL;
	int tag = 0;
	enum iofields blocknr;
	MPI_Status status;

	if (Task.Rank == ReadTask) {

		fd = fopen(filename, "r");

		uint64_t ntot = read_gadget_head(fd);

		printf("\nReading file <%s> on Task <%i-%i> \n"
		       "   Sph   <%9lli>  DM     <%9lli>    \n"
		       "   Disk  <%9lli>  Bulge  <%9lli>    \n"
		       "   Star  <%9lli>  Bndry  <%9lli>    \n"
		       "   Total <%9lli> \n\n",
		       filename, ReadTask, LastTask,
		       Header.Npart[0], Header.Npart[1], Header.Npart[2],
		       Header.Npart[3], Header.Npart[4], Header.Npart[5], ntot);
		fflush(stdout);

		for (task = ReadTask + 1; task <= LastTask; task++) {
			MPI_Ssend(&Header, sizeof(Header), MPI_BYTE, task, tag,
				  MPI_COMM_WORLD);
		}
	} else {
		MPI_Recv(&Header, sizeof(Header), MPI_BYTE, ReadTask, tag,
			 MPI_COMM_WORLD, &status);
	}

	/* Set Snapshot Properties */
	Snap.Boxsize = Header.BoxSize;
	Snap.Redshift = Header.Redshift;
	Snap.Time = Header.Time;

	for (i = 0; i < N_part_types; i++) {
		Snap.Npart[i] += Header.Npart[i];
		Snap.Masstab[i] = Header.Mass[i] * 1 / Cosmo.h;
	}

	MPI_Bcast(&Snap, sizeof(Snap), MPI_BYTE, 0, MPI_COMM_WORLD);

	/* Set Comoving units for Gadget */
	Comv2phys.Length = 1 / (1 + Snap.Redshift) / Cosmo.h;
	Comv2phys.Mass = 1 / Cosmo.h;
	Comv2phys.Vel = 1. / sqrt(1 + Snap.Redshift);

	/* Determine particle distribution over CPUs */
	nTask = LastTask - ReadTask + 1;
	for (i = 0; i < N_part_types; i++) {

		for (j = Task.Rank - ReadTask; j < Header.Npart[i]; j += nTask)
			nRead[i]++;

		nReadTot += nRead[i];
	}

	Reallocate_P(nReadTot, nRead, +1);

	/* Shift collisionless particles if multiple files */
	src = Task.Npart[0] - nRead[0];
	target = src + nReadTot;
	nBytes = (Task.PartTotal - nReadTot - src) * sizeof(*P);

	memmove(&P[target].Pos[0], &P[src].Pos[0], nBytes);

	/* Read blocks  */
	for (blocknr = IO_POS; blocknr < IO_LASTENTRY; blocknr++) {

		set_block_prop(blocknr);
		int blockExist = 0;

		if (!block_required(Block.Label))
			continue;

		if (Task.Rank == ReadTask) {

			blockExist = find_block(fd, Block.Label);
			
			int sendExist = blockExist;

			for (task = ReadTask + 1; task < LastTask + 1; task++)
				MPI_Ssend(&sendExist, 1, MPI_LONG, task, tag,
					  MPI_COMM_WORLD);

			if (blockExist > 0) {

				printf("%28s %8d Mb\n", Block.Name,
				       blockExist / 1024 / 1024);

				comm_buf =
				    Malloc(Block.Ntot *
					   Block.Bytes_per_element);

				blockExist =
				    read_gadget_block(comm_buf, Block.Label, fd,
						      Block.Data_type);
			} else
				for (i = 0; i < 1000; i++)
					Assert(blocknr == IO_MASS
					       || strncmp(Effect.Req.Block[i],
							  Block.Label, 5),
					       "Required Block not found: %s",
					       Block.Label);
		} else
			MPI_Recv(&blockExist, 1, MPI_LONG, ReadTask, tag,
				 MPI_COMM_WORLD, &status);

		byteOffset = 0;

		if (Task.Rank == ReadTask && blockExist) {

			for (i = 0; i < N_part_types; i++) {

				if (!Block.Npart[i])
					continue;

				byteOffset +=
				    nRead[i] * Block.Bytes_per_element;

				for (task = ReadTask + 1; task < LastTask + 1;
				     task++) {

					MPI_Recv(nSend, N_part_types,
						 MPI_LONG_LONG, task, tag,
						 MPI_COMM_WORLD, &status);

					nBytes =
					    nSend[i] * Block.Bytes_per_element;

					MPI_Ssend(comm_buf + byteOffset, nBytes,
						  MPI_BYTE, task, tag,
						  MPI_COMM_WORLD);

					byteOffset += nBytes;
				}
			}
		} else if (blockExist) {

			comm_buf = Malloc(Block.Bytes_per_element * nReadTot);

			for (i = 0; i < N_part_types; i++) {

				if (!Block.Npart[i])
					continue;

				MPI_Ssend(nRead, N_part_types, MPI_LONG_LONG,
					  ReadTask, tag, MPI_COMM_WORLD);

				nBytes = nRead[i] * Block.Bytes_per_element;

				MPI_Recv(comm_buf + byteOffset, nBytes,
					 MPI_BYTE, ReadTask, tag,
					 MPI_COMM_WORLD, &status);

				byteOffset += nBytes;
			}
		}

		partOffset = Task.Npart[0] - nRead[0];
		bufOffset = 0;

		for (i = 0; i < N_part_types; i++) {

			if (Block.Npart[i] && blockExist) {

				for (j = 0; j < nRead[i]; j++)
					empty_comm_buffer((enum iofields)
							  blocknr, comm_buf,
							  partOffset + j,
							  (bufOffset +
							   j) *
							  Block.
							  Val_per_element);

				if (blocknr == IO_ID)	// set type 
					for (j = 0; j < nRead[i]; j++)
						P[partOffset + j].Type = i;

			} else if (blocknr == IO_MASS)
				for (j = 0; j < nRead[i]; j++)	// masses from Header 
					P[partOffset + j].Mass =
					    Snap.Masstab[i];

			if (Block.Npart[i]) {

				if (Task.Rank == ReadTask)
					bufOffset += Block.Npart[i];
				else
					bufOffset += nRead[i];

				partOffset += nRead[i];
			}
		}

		if (blockExist)
			Free(comm_buf);
	}

	if (Task.Rank == ReadTask)
		fclose(fd);

	/* treat AREPO VOL -> HSML */
	if (Snap.Have_Arepo)
		for (j = 0; j < nRead[0]; j++)
			P[j].Hsml = pow(P[j].Hsml / fourpithirds, 1.0 / 3.0);

	return;
}

/*Determine number of files to read
 * */
int find_files(char *filename)
{
	char buf[MAXLINELENGTH];
	int nFiles = 0;
	FILE *fd = NULL;

	if (!(fd = fopen(filename, "r"))) {
		for (;;) {
			sprintf(buf, "%s.%i", filename, nFiles);

			if (!(fd = fopen(buf, "r")))
				break;

			fclose(fd);

			nFiles++;

			Assert(nFiles < 1000, "Found too many files");

		}

		Assert(nFiles > 0, "Can't open input file %s", filename);

		rprintf(" \n Found <%i> file(s) ! \n\n", nFiles);

	} else {
		nFiles = 1;
	}

	return (nFiles);
}

/* Basic routine to read data from a file 
 * */
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
	size_t nRead;

	if ((nRead = fread(ptr, size, nmemb, stream)) != nmemb) {
		if (feof(stream)) {
			nRead = 0;	/*EOF reached */
		} else
			Assert(0, "I/O error (fread) ");
	}
	return nRead;
}

/* Routine to Swap ENDIAN
 * */
void swap_Nbyte(char *data, int n, int m)
{
	int i, j;
	char old_data[16];

	if (Swap > 0) {
		for (j = 0; j < n; j++) {
			memcpy(&old_data[0], &data[j * m], m);
			for (i = 0; i < m; i++) {
				data[j * m + i] = old_data[m - i - 1];
			}
		}
	}
	return;
}

int find_block(FILE * fd, char *label)
{
	int4bytes blocksize = 0;
	char blocklabel[4] = { "    " };

	rewind(fd);

	while (!feof(fd) && blocksize == 0) {

		READF77HEAD;

		if (f77block == 134217728) {

			printf("\nEnabling ENDIAN swapping !\n");

			Swap = 1 - Swap;
			swap_Nbyte((char *)&f77block, 1, 4);
		}

		Assert(f77block == 8, "incorrect FORTRAN77 file format");

		if (my_fread(blocklabel, 4 * sizeof(char), 1, fd)) {

			my_fread(&blocksize, sizeof(int4bytes), 1, fd);
			swap_Nbyte((char *)&blocksize, 1, 4);

			READF77HEAD;

			if (strncmp(label, blocklabel, 4) != 0) {	// go forward 
				fseek(fd, blocksize, 1);
				blocksize = 0;
			}

		} else {
			blocksize = 8;
			break;
		}
	}
	return (blocksize - 8);
}

/* Read the header information, returns total nr of particles in file
 * */
int read_gadget_head(FILE * fd)
{
	int blocksize, dummysize, i;
	unsigned int Npart[N_part_types], Nall[N_part_types],
	    NallHW[N_part_types];
	int ntot = 0;

	blocksize = find_block(fd, "HEAD");

	Assert(blocksize > 0, "Header not found");

	dummysize = blocksize - 2 * N_part_types * sizeof(int) -
	    4 * sizeof(long) - 12 * sizeof(double);
	READF77HEAD;

	my_fread(Npart, N_part_types * sizeof(unsigned int), 1, fd);
	swap_Nbyte((char *)Npart, N_part_types, 4);

	my_fread(Header.Mass, N_part_types * sizeof(double), 1, fd);
	swap_Nbyte((char *)Header.Mass, N_part_types, 8);

	my_fread((void *)&Header.Time, sizeof(double), 1, fd);
	swap_Nbyte((char *)&Header.Time, 1, 8);

	my_fread((void *)&Header.Redshift, sizeof(double), 1, fd);
	swap_Nbyte((char *)&Header.Redshift, 1, 8);

	my_fread((void *)&Header.FlagSfr, sizeof(int), 1, fd);
	swap_Nbyte((char *)&Header.FlagSfr, 1, 4);

	my_fread((void *)&Header.FlagFeedback, sizeof(int), 1, fd);
	swap_Nbyte((char *)&Header.FlagFeedback, 1, 4);

	my_fread(Nall, N_part_types * sizeof(unsigned int), 1, fd);
	swap_Nbyte((char *)Nall, N_part_types, 4);

	my_fread((void *)&Header.FlagCooling, sizeof(int), 1, fd);
	swap_Nbyte((char *)&Header.FlagCooling, 1, 4);

	my_fread((void *)&Header.NumFiles, sizeof(int), 1, fd);
	swap_Nbyte((char *)&Header.NumFiles, 1, 4);

	my_fread((void *)&Header.BoxSize, sizeof(double), 1, fd);
	swap_Nbyte((char *)&Header.BoxSize, 1, 8);

	my_fread((void *)&Header.Omega0, sizeof(double), 1, fd);
	swap_Nbyte((char *)&Header.Omega0, 1, 8);

	my_fread((void *)&Header.OmegaLambda, sizeof(double), 1, fd);
	swap_Nbyte((char *)&Header.OmegaLambda, 1, 8);

	my_fread((void *)&Header.HubbleParam, sizeof(double), 1, fd);
	swap_Nbyte((char *)&Header.HubbleParam, 1, 8);

	my_fread((void *)&Header.FlagAge, sizeof(int), 1, fd);
	swap_Nbyte((char *)&Header.FlagAge, 1, 8);

	my_fread((void *)&Header.FlagMetals, sizeof(int), 1, fd);
	swap_Nbyte((char *)&Header.FlagMetals, 1, 8);

	my_fread((void *)NallHW, N_part_types * sizeof(unsigned int), 1, fd);
	swap_Nbyte((char *)NallHW, N_part_types, 4);

	if (NallHW[0] != 0)
		printf("Nall HW not well tested !! ");

	fseek(fd, dummysize, 1);
	READF77HEAD;

	for (i = 0; i < N_part_types; i++) {	// HighWord 
		Header.Npart[i] = Npart[i];
		Header.Nall[i] =
		    (long long)(Nall[i] + (((long long)NallHW[i]) << 32));
		ntot += Header.Npart[i];
	}

	return ntot;
}

int read_gadget_block(void *data, char *label, FILE * fd, size_t sizeof_type)
{
	unsigned int blocksize = 0;

	blocksize = find_block(fd, label);

	Assert(blocksize, "Block not found");

	READF77HEAD;
	my_fread(data, blocksize, 1, fd);
	swap_Nbyte((char *)data, blocksize / sizeof_type, 4);
	READF77HEAD;

	return (blocksize);
}

/* Set Block characteristics 
 * "You are all different" - 
 * "We are all different" - 
 * "I'm not !" 
 * 				(life of brian) */
void set_block_prop(enum iofields blocknr)
{
	int i;

	for (i = 0; i < N_part_types; i++)
		Block.Npart[i] = 0;

	switch (blocknr) {
	case IO_POS:
		Block.Label = "POS ";	/* Has to be 4 Letters */
		Block.Name = "Coordinates";
		for (i = 0; i < N_part_types; i++)
			Block.Npart[i] = Header.Npart[i];	/* Particle type(s) of block */
		Block.Val_per_element = 3;	/* 3=vector, 1=scalar */
		Block.Data_type = FLOAT;	/* Length in bytes */
		Block.Rmv_comoving = Comv2phys.Length;	/* factor to make not comoving */
		break;
	case IO_VEL:
		Block.Label = "VEL ";
		Block.Name = "Velocities";
		for (i = 0; i < N_part_types; i++)
			Block.Npart[i] = Header.Npart[i];
		Block.Val_per_element = 3;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = Comv2phys.Vel;
		break;
	case IO_ID:
		Block.Label = "ID  ";
		Block.Name = "Particle IDs";
		for (i = 0; i < N_part_types; i++)
			Block.Npart[i] = Header.Npart[i];
		Block.Val_per_element = 1;
		Block.Data_type = UINT;
		Block.Rmv_comoving = 1;
		break;
	case IO_U:
		Block.Label = "U   ";
		Block.Name = "Internal Energy";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 1;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = 1;
		break;
	case IO_RHO:
		Block.Label = "RHO ";
		Block.Name = "Density";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 1;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving =
		    (Comv2phys.Mass / pow(Comv2phys.Length, 3));
		break;
	case IO_HSML:
		Block.Label = "HSML";
		Block.Name = "Smoothing Length";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 1;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = Comv2phys.Length;
		break;
	case IO_VOL:		/* AREPO only */
		Block.Label = "VOL ";
		Block.Name = "AREPO Cell Volume";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 1;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = pow(Comv2phys.Length, 3);
		break;
	case IO_BFLD:
		Block.Label = "BFLD";
		Block.Name = "Magnetic Field";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 3;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = 1;
		break;
	case IO_DIVB:
		Block.Label = "DIVB";
		Block.Name = "Magnetic Field Divergence";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 1;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = 1;
		break;
	case IO_MASS:
		Block.Label = "MASS";
		Block.Name = "Particle Mass";
		for (i = 0; i < N_part_types; i++)
			Block.Npart[i] = Header.Npart[i];
		Block.Val_per_element = 1;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = Comv2phys.Mass;
		break;
	case IO_VELT:
		Block.Label = "VELT";
		Block.Name = "Turbulent Velocity, part";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 1;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = Comv2phys.Vel;
		break;
	case IO_VRMS:
		Block.Label = "VRMS";
		Block.Name = "Turbulent Velocity, mean";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 1;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = Comv2phys.Vel;
		break;
	case IO_VBULK:
		Block.Label = "VBULK";
		Block.Name = "Local Mean Velocity";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 3;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = Comv2phys.Vel;
		break;
	case IO_VTAN:
		Block.Label = "VTAN";
		Block.Name = "Local Mean Velocity Tangential Component";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 1;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = Comv2phys.Vel;
		break;
	case IO_VRAD:
		Block.Label = "VRAD";
		Block.Name = "Local Mean Velocity Radial Component";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 1;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = Comv2phys.Vel;
		break;
	case IO_TNGB:
		Block.Label = "TNGB";
		Block.Name = "True Number Of Neighbours";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 1;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = 1;
		break;
	case IO_VDIV:
		Block.Label = "VDIV";
		Block.Name = "Local Velocity Divergence";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 1;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = Comv2phys.Vel / Comv2phys.Length;
		break;
	case IO_VROT:
		Block.Label = "VROT";
		Block.Name = "Local Velocity Curl";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 1;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = Comv2phys.Vel / Comv2phys.Length;
		break;
	case IO_DPP:
		Block.Label = "DPP ";
		Block.Name = "MagnetosReaccCoefficient";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 1;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = 1;
		break;
	case IO_nHII:
		Block.Label = "nHII";
		Block.Name = "IonisedFractionHydrogen";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 1;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = 1;
		break;
	case IO_nHeII:
		Block.Label = "nHe2";
		Block.Name = "IonisedFractionHelium2";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 1;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = 1;
		break;
	case IO_nHeIII:
		Block.Label = "nHe3";
		Block.Name = "IonisedFactionHelium3";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 1;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = 1;
		break;
	case IO_MACH:
		Block.Label = "MACH";
		Block.Name = "Mach Number";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 1;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = 1;
		break;
	case IO_VORT:
		Block.Label = "VORT";
		Block.Name = "Vorticity";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 3;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = 1;
		break;
	case IO_ABVC:
		Block.Label = "ABVC";
		Block.Name = "AlphaViscosity";
		Block.Npart[0] = Header.Npart[0];
		Block.Val_per_element = 1;
		Block.Data_type = FLOAT;
		Block.Rmv_comoving = 1;
		break;
		/*Add above, not below !! */
	case IO_LASTENTRY:
		Block.Label = "LAST";
		Block.Name = "";
		Block.Val_per_element = 0;
		Block.Data_type = 0;
		break;
	}
	Block.Bytes_per_element = Block.Data_type * Block.Val_per_element;

	for (i = Block.Ntot = 0; i < N_part_types; i++)
		Block.Ntot += Block.Npart[i];

	return;
}

/*Fill P and Gas with data buffer 'fp'. 
 * */
void empty_comm_buffer(enum iofields blocknr, void *fp, size_t iP, size_t ifp)
{
	int i = 0;

	switch (blocknr) {
	case IO_POS:
		for (i = 0; i < 3; i++)
			P[iP].Pos[i] =
			    ((float *)fp)[ifp + i] * Block.Rmv_comoving;
		break;
	case IO_VEL:
		for (i = 0; i < 3; i++)
			P[iP].Vel[i] =
			    ((float *)fp)[ifp + i] * Block.Rmv_comoving;
		break;
	case IO_ID:
		P[iP].ID = ((unsigned int *)fp)[ifp];
		break;
	case IO_U:
		Gas[iP].U = ((float *)fp)[ifp];
		break;
	case IO_RHO:
		P[iP].Rho = ((float *)fp)[ifp] * Block.Rmv_comoving;
		break;
	case IO_VOL:		/* Store VOL in HSML for now */
		Snap.Have_Arepo = 1;
	case IO_HSML:
		P[iP].Hsml = ((float *)fp)[ifp] * Block.Rmv_comoving;
		break;
	case IO_BFLD:
		for (i = 0; i < 3; i++)
			Gas[iP].Bfld[i] =
			    ((float *)fp)[ifp + i] * Block.Rmv_comoving;
		break;
	case IO_DIVB:
#ifdef DIVB
		Gas[iP].DivB = ((float *)fp)[ifp] * Block.Rmv_comoving;
#endif
		break;
	case IO_MASS:
		P[iP].Mass = ((float *)fp)[ifp] * Block.Rmv_comoving;
		break;
	case IO_VELT:
#ifdef VTURB
		Gas[iP].VTurb = ((float *)fp)[ifp] * Block.Rmv_comoving;
#endif
		break;
	case IO_VRMS:
#ifdef VTURB
		Gas[iP].VRms = ((float *)fp)[ifp] * Block.Rmv_comoving;
#endif
		break;
	case IO_VTAN:
#ifdef VTURB
		Gas[iP].VTan = ((float *)fp)[ifp] * Block.Rmv_comoving;
#endif
		break;
	case IO_VRAD:
#ifdef VTURB
		Gas[iP].VRad = ((float *)fp)[ifp] * Block.Rmv_comoving;
#endif
		break;
	case IO_VBULK:
#ifdef VTURB
		for (i = 0; i < 3; i++)
			Gas[iP].VBulk[i] =
			    ((float *)fp)[ifp + i] * Block.Rmv_comoving;
#endif
		break;
	case IO_TNGB:
#ifdef VTURB
		Gas[iP].TNgb = ((float *)fp)[ifp];
#endif
		break;
	case IO_VDIV:
#ifdef VTURB
		Gas[iP].DivVel = ((float *)fp)[ifp] * Block.Rmv_comoving;
#endif
		break;
	case IO_VROT:
#ifdef VTURB
		Gas[iP].CurlVel = ((float *)fp)[ifp] * Block.Rmv_comoving;
#endif
		break;
	case IO_DPP:
#ifdef VTURB
		Gas[iP].Dpp = ((float *)fp)[ifp] * Block.Rmv_comoving;
#endif
		break;
	case IO_nHII:
#ifdef RADTRANSFER
		Gas[iP].nHII = ((float *)fp)[ifp] * Block.Rmv_comoving;
#endif
		break;
	case IO_nHeII:
#ifdef RADTRANSFER
		Gas[iP].nHeII = ((float *)fp)[ifp] * Block.Rmv_comoving;
#endif
		break;
	case IO_nHeIII:
#ifdef RADTRANSFER
		Gas[iP].nHeIII = ((float *)fp)[ifp] * Block.Rmv_comoving;
#endif
		break;
	case IO_MACH:
#ifdef SHOCKS
		Gas[iP].NMach = ((float *)fp)[ifp] * Block.Rmv_comoving;
#endif
		break;
	case IO_VORT:
#ifdef VORTICITY
		for (i = 0; i < 3; i++)
			Gas[iP].Vorticity[i] =
			    ((float *)fp)[ifp] * Block.Rmv_comoving;
#endif
	case IO_ABVC:
#ifdef TIME_DEP_ART_VISC
		Gas[iP].AlphaVisc = ((float *)fp)[ifp] * Block.Rmv_comoving;
#endif
		/*Add above, not below !! */
	case IO_LASTENTRY:
		break;
	}
	return;
}

int block_required(char *label)
{
	const int maxBlocks = 999;

	for (int i = 0; i < maxBlocks; i++) 
		if (!strncmp(label, Effect.Req.Block[i], 4))
			return 1;

	return 0;
}

#endif
