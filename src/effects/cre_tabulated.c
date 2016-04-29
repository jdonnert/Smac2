
#include "../proto.h"
#include "../globals.h"
#include "effects.h"

#define FBASE "spec_uncompressed/spec"
#define MAXFILES 10000

static struct CR_Spectrum_Header {
	long SnapNum;
	long long Nbins;
	unsigned long long Nall;
	unsigned long long Npart;
	long long Nfiles;
	unsigned long long StartID;
	double Pmin;
	double Pmax;
	double Plow;
	double Phigh;
	long FlagCompressed;
	long SpecSizeBytes;
	char Fill[84];
} Head;

static float *cr_spectrum, *comm_buf;
static double *p;

static double Emax, Emin, pstep;

static size_t find_files();

/* return N(E) dE  */
double cre_spectrum_tabulated(double E, int ipart)
{
	if (E < Emin || E > Emax) // out of range 
		return 0;

	size_t ibin = ipart * Head.Nbins + floor(log10(E / Emin) / pstep);
	
	return pow(10, cr_spectrum[ibin]); // n(p): dE/dp = c
}

/* spectrum input */
void set_tabulated_factors()
{
	size_t nBytes, nRead, nFloat, nFiles = 0, ipart;
	FILE *fp = NULL;

	//prepare_synchrotron();

	if (Task.Rank == 0)
		nFiles = find_files();

	rprintf("Reading binary CR spectra from <%zu> files : \n", nFiles);

	MPI_Bcast(&nFiles, sizeof(nFiles), MPI_BYTE, 0, MPI_COMM_WORLD);

	for (int ifile = 0; ifile < nFiles; ifile++) {

		if (Task.Rank == 0) { // Header

			char fname[MAXLINELENGTH] = { "" };

			if (nFiles == 1) 
				sprintf(fname, "%s_%03ld", FBASE, Snap.SnapNum);
			else 
				sprintf(fname, "%s_%03ld.%d", FBASE, Snap.SnapNum, ifile);

			if (!(fp = fopen(fname, "r")))
				Assert(0, "Can't open CRe spectrum file");

			nRead = fread(&Head, sizeof(Head), 1, fp);

			Assert(nRead == 1, "Could not read Header in file %s",
			       fname);

			Assert(Head.FlagCompressed == 0,
			       "File %s contains compressed spectra, rerun with subflag 10",
			       fname);

			printf("<%lld> spectra from <%s> \n", Head.Npart, fname);
		}

		MPI_Barrier(MPI_COMM_WORLD);

		MPI_Bcast(&Head, sizeof(Head), MPI_BYTE, 0, MPI_COMM_WORLD);

		/* Spectrum */
		nBytes = Head.Npart * Head.Nbins * sizeof(*comm_buf);
		comm_buf = Malloc(nBytes);

		if (Task.Rank == 0) {

			nRead = fread(comm_buf, nBytes, 1, fp);

			fclose(fp);
		}

		MPI_Bcast(comm_buf, nBytes, MPI_BYTE, 0, MPI_COMM_WORLD);

		if (ifile == 0) {	// Init spectrum 

			nFloat = Task.Npart[0] * Head.Nbins;
			nBytes = nFloat * sizeof(*cr_spectrum);

			if (cr_spectrum != NULL)
				Free(cr_spectrum);
;
			cr_spectrum = Malloc(nBytes);
		}

		/* copy comm_buffer to spectrum - pick out parts on this CPU */
		nBytes = Head.Nbins * sizeof(*comm_buf);

		for (ipart = 0; ipart < Task.Npart[0]; ipart++) {

			size_t idx_buf = (P[ipart].ID - Head.StartID) * Head.Nbins;

			if (idx_buf < 0)
				continue;

			if (idx_buf > Head.Npart * Head.Nbins)
				continue;

			size_t idx_spec = ipart * Head.Nbins;

			memcpy(&cr_spectrum[idx_spec], &comm_buf[idx_buf], nBytes);
		}

		Free(comm_buf);
	}

	rprintf("Nbins  = %lld \nNall  = %lld \n"
		"Pmin   = %g  \nPmax   = %g \n"
		"Plow   = %g  \nPhigh  = %g \n\n", Head.Nbins, Head.Nall,
		Head.Pmin, Head.Pmax, Head.Plow, Head.Phigh);

	/* momentum grid */
	nBytes = Head.Nbins * sizeof(*p);

	if (p == NULL)
		p = Malloc(nBytes);

	pstep = log10(Head.Pmax / Head.Pmin) / (Head.Nbins - 1);

	for (int i = 0; i < Head.Nbins; i++)
		p[i] = Head.Pmin * pow(10, i * pstep) * m_e * c;

	Emin = Head.Pmin * m_e * c * c;
	Emax = Head.Pmax * m_e * c * c;

	MPI_Barrier(MPI_COMM_WORLD);

	return;
}

static size_t find_files()
{
	char fname[MAXLINELENGTH];
	FILE *fp;

	sprintf(fname, "%s_%03ld", FBASE, Snap.SnapNum);

	if ((fp = fopen(fname, "r"))) {

		fclose(fp);

		return 1;
	}

	for (int i = 0; i < MAXFILES; i++) {

		sprintf(fname, "%s_%03ld.%d", FBASE, Snap.SnapNum, i);

		if (!(fp = fopen(fname, "r"))) {

			Assert(i, "Can't find binary CR spectrum file %s_%03ld "
			       "or %s", FBASE, Snap.SnapNum, fname);

			return i;
		}

		fclose(fp);
	}

	Assert(0, "Too many CR spectrum files: MAXFILES !");

	return -1;		// never reached 
}

#undef MAXFILES
#undef FBASE
