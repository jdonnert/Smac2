#include "../globals.h"
#include "effects.h"
#include "cre_compressed.h"

#define FBASE "spec/spec"
#define MAXFILES 10000

#define SAMPLING 128		// Length of arrays holding curve
#define MAXKNOTS 12
#define SPECSIZE_BYTES 60

int uncompress_knots_binary(const char *, struct Knot *);
void draw_curve(const struct Knot *, const int, double *);
float uncompressFloat_8bit(uint8_t);
float uncompressFloat_16bit(uint16_t);

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

static double log_p[SAMPLING] = { 0 };

static double Emax = 0, Emin = 0, pstep = 0;

static char *compressed_data = NULL;	// compressed data from file

static double current_spectrum[SAMPLING] = { 0 };	// one decompressed spec
static int last_ipart = -1;
#pragma omp threadprivate(current_spectrum, last_ipart)

double cre_spectrum_compressed(double E, int ipart)
{
	if (E < Emin || E > Emax) // out of range
		return 0;	

	if (ipart != last_ipart) {	// decompress new spectrum

		last_ipart = ipart;
	
		const size_t nBytes = (Head.SpecSizeBytes + sizeof(float));

		char *nCRe = &compressed_data[ipart * nBytes];
		char *src_data = &compressed_data[ipart * nBytes + sizeof(float)];

		Uncompress(*((float *)nCRe), src_data, current_spectrum);

		int ibin = floor(log10(E / Emin) / pstep);
	}

	int ibin = floor(log10(E / Emin) / pstep);

	return current_spectrum[ibin];
}

/* Readin and distribution of CRe spectra. Decompression is on the fly 
 * We send the whole file to everyone, then all tasks pick there IDs out 
 * and write them into compressed_data */
void setup_decompression()
{
	const int readTask = 0;

//	prepare_synchrotron();

	int nFiles = 0;
	if (Task.Rank == readTask)
		nFiles = find_files();

	MPI_Bcast(&nFiles, 1, MPI_INT, readTask, MPI_COMM_WORLD);

	rprintf("Reading compressed CR spectra from <%d> files : \n\n", nFiles);

	for (int i = 0; i < nFiles; i++) {

		FILE *fp = NULL;
		char fname[MAXLINELENGTH];

		if (Task.Rank == readTask) {	// open & header

			sprintf(fname, "%s_%03ld", FBASE, Snap.SnapNum);

			if (nFiles > 1)
				sprintf(fname, "%s.%d", fname, i);

			fp = fopen(fname, "r");

			Assert(fp != NULL, "Can't open CRe spectrum file");

			int nRead = fread(&Head, sizeof(Head), 1, fp);

			Assert(nRead == 1, "Could not read Header in file %s",
			       fname);

			Assert(Head.FlagCompressed,
			       "File %s contains uncompressed spectra, "
			       "rerun with subflag 1", fname);
		}

		MPI_Bcast(&Head, sizeof(Head), MPI_BYTE, readTask,
			  MPI_COMM_WORLD);

		const size_t nBytes_per_spec =
		    Head.SpecSizeBytes + sizeof(float);

		if (compressed_data == NULL) {

			size_t nBytes = Task.Npart[0] * nBytes_per_spec;

			compressed_data = Malloc(nBytes);
		}

		size_t nBytes = Head.Npart * nBytes_per_spec;
		char *readBuf = Malloc(nBytes);

		rprintf("  %8lld spectra from file %s (%zu MB)\n",
			Head.Npart, fname, nBytes / 1024 / 1024);

		if (Task.Rank == readTask)
			fread(readBuf, nBytes, 1, fp);

		MPI_Bcast(readBuf, nBytes, MPI_BYTE, readTask, MPI_COMM_WORLD);

		#pragma omp parallel for
		for (size_t ipart = 0; ipart < Task.Npart[0]; ipart++) {	// pick out

			size_t idx =
			    (P[ipart].ID - Head.StartID) * nBytes_per_spec;

			if ((idx >= 0) && (idx < Head.Npart * nBytes_per_spec))
				memcpy(&(compressed_data[ipart * nBytes_per_spec]),
				       &(readBuf[idx]), nBytes_per_spec);
		}

		free(readBuf);
	}

	rprintf("\n  Nbins  = %lld \n  Nall   = %lld \n"
		"  Pmin   = %g \n  Pmax   = %g \n"
		"  Plow   = %g \n  Phigh  = %g \n"
		"  Size   = %zu Bytes\n\n", Head.Nbins, Head.Nall, Head.Pmin,
		Head.Pmax, Head.Plow, Head.Phigh, Head.SpecSizeBytes);

	pstep = log10(Head.Pmax / Head.Pmin) / (SAMPLING - 1);

	for (int i = 0; i < SAMPLING; i++)
		log_p[i] = log10(Head.Pmin) + i * pstep;

	Emin = Head.Pmin * m_e * c * c;
	Emax = Head.Pmax * m_e * c * c;

	return;
}

static size_t find_files()
{
	char fname[MAXLINELENGTH];
	sprintf(fname, "%s_%03ld", FBASE, Snap.SnapNum);

	FILE *fp = NULL;

	if ((fp = fopen(fname, "r"))) {

		fclose(fp);

		return 1;
	}

	for (int i = 0; i < MAXFILES; i++) {

		sprintf(fname, "%s_%03ld.%d", FBASE, Snap.SnapNum, i);

		FILE *fp;

		if (!(fp = fopen(fname, "r"))) {

			Assert(i,
			       "Can't find binary CR spectrum file %s_%03ld or %s",
			       FBASE, Snap.SnapNum, fname);

			return i;
		}

		fclose(fp);
	}

	Assert(0, "Too many CR spectrum files: MAXFILES !");

	return -1;		// never reached 
}

void Uncompress(const float nCRe, const char *spectrum, double *out_spectrum)
{
	Assert(out_spectrum != NULL, "Can't uncompress into NULL");
	Assert(spectrum != NULL, "Can't uncompress from NULL");

	struct Knot K[MAXKNOTS] = { {0} };

	int nKnots = uncompress_knots_binary(spectrum, K);

	double np[SAMPLING] = { 0 };
	
	if (nKnots != 0) {

		for (int i = 0; i < SAMPLING; i++)
			np[i] = -FLT_MAX;

		draw_curve(K, nKnots, np);	// fill raw spectrum
		
		for (int i = 0; i < SAMPLING; i++) 
			np[i] = pow(10, np[i] * nCRe);

	}

	memcpy((void*) out_spectrum, (void *) np, SAMPLING * sizeof(*out_spectrum));

	return;
}

/* find y(x), cubic Hermite spline */
void draw_curve(const struct Knot *K, const int nKnots, double *Curve) 
{
    double c0[2], c1[2], c2[2], c3[2];
    
	for (int i = 0; i < nKnots-1; i++) { // def. Hermitian polynome (cspline)
        
        c0[0] = 2*K[i].P[0] - 2*K[i+1].P[0] + K[i].Mright[0] + K[i+1].Mleft[0];
        c0[1] = 2*K[i].P[1] - 2*K[i+1].P[1] + K[i].Mright[1] + K[i+1].Mleft[1];

        c1[0] = -3*K[i].P[0]+3*K[i+1].P[0]-2*K[i].Mright[0] - K[i+1].Mleft[0];
        c1[1] = -3*K[i].P[1]+3*K[i+1].P[1]-2*K[i].Mright[1] - K[i+1].Mleft[1];

        c2[0] = K[i].Mright[0];
        c2[1] = K[i].Mright[1];

        c3[0] = K[i].P[0];
        c3[1] = K[i].P[1];

		double dPinv = 1.0/(K[i+1].P[0]-K[i].P[0]);

        for (int j = K[i].idx; j < K[i+1].idx+1; j++) {

			double t = (log_p[j]-K[i].P[0])*dPinv;
			int it = 0;

			for (;;) { // Newton-Raphson, damn it's fast 

				double q = c0[0]*t*t*t + c1[0]*t*t + c2[0]*t + c3[0];
				t -= (q - log_p[j]) / (3*c0[0]*t*t + 2 * c1[0]*t + c2[0]);

				if ( fabs(q-log_p[j]) / log_p[j] < 1e-2 || it++ > 50) 
					break;
			}
	
			Curve[j] = c0[1]*t*t*t + c1[1]*t*t + c2[1]*t + c3[1];
        }
    }

	return;
}

int uncompress_knots_binary(const char *spectrum, struct Knot *K)
{
    struct half_point hp = { 0 };
    struct full_point fp = { 0 };

    const char *src = spectrum;
        
    memcpy(&hp, src, sizeof(hp)); // KNOT_START
    src += sizeof(hp);    

    K[0].type = KNOT_START;

    size_t idx = K[0].idx = hp.x;

    K[0].P[0] = log_p[idx];
    K[0].P[1] = uncompressFloat_16bit(hp.y);

	K[0].Mright[0] = uncompressFloat_8bit(hp.M.xyLR[0]);
	K[0].Mright[1] = uncompressFloat_8bit(hp.M.xyLR[1]);

	if (K[0].idx == 0)
		return 0;

	int nKnots = 1;           
    size_t KnotMemSize = sizeof(hp);

	int nMinMaxKnots = 0;

    for (int i = 1; i < MAXKNOTS; i++) {
        
        if (*src & 0x80) { // marker bit set
            
			K[i].type = KNOT_FULL;
            
            memcpy(&fp, src, sizeof(fp));

            src += sizeof(fp);
            KnotMemSize += sizeof(fp);

            idx = K[i].idx = 0x7F & fp.x; // remove marker bit

            K[i].P[0] = log_p[idx];
            K[i].P[1] = uncompressFloat_16bit(fp.y);

            K[i].Mleft[0] = uncompressFloat_8bit(fp.xleft);
            K[i].Mleft[1] = uncompressFloat_16bit(fp.yleft);

            K[i].Mright[0] = uncompressFloat_8bit(fp.xright);
            K[i].Mright[1] = uncompressFloat_16bit(fp.yright);

        } else { // marker bit not set 
            
			K[i].type = KNOT_MINMAX;

            memcpy(&hp, src, sizeof(hp));
            
            src += sizeof(hp);
            KnotMemSize += sizeof(hp);

            idx = K[i].idx = hp.x;

            K[i].P[0] = log_p[idx];
            K[i].P[1] = uncompressFloat_16bit(hp.y);

            K[i].Mleft[0] = uncompressFloat_8bit(hp.M.xLR[0]);
            K[i].Mleft[1] = 0;           

            K[i].Mright[0] = uncompressFloat_8bit(hp.M.xLR[1]);
            K[i].Mright[1] = 0;

			nMinMaxKnots++;
        } 
        
        nKnots++;           

        if (nMinMaxKnots % 2 
			&& KnotMemSize + sizeof(fp)+sizeof(hp) > SPECSIZE_BYTES ) 
            break; // only KNOT_STOP possible, exit loop
    }
    
    memcpy(&hp, src, sizeof(hp)); // KNOT_STOP

	if (hp.x == 0) {// we overshot if START == MINMAX 
		
		src -= sizeof(hp);

		memcpy(&hp, src, sizeof(hp));

		nKnots--;
	}

    K[nKnots].type = KNOT_STOP;
    
    K[nKnots].idx = idx = hp.x;

    K[nKnots].P[0] = log_p[idx];
    K[nKnots].P[1] = uncompressFloat_16bit(hp.y);
     
	K[nKnots].Mleft[0] = uncompressFloat_8bit(hp.M.xyLR[0]);
	K[nKnots].Mleft[1] = uncompressFloat_8bit(hp.M.xyLR[1]);

    nKnots++;           

    return nKnots;
}


float uncompressFloat_8bit(uint8_t input)
{
	return input * 0.0156250 - 1;	// 2^-6
}

float uncompressFloat_16bit(uint16_t input)
{
	return input * 0.000244140625 - 8;	// * 2^-12
}

#undef FBASE
#undef MAXFILES

#undef SAMPLING
#undef MAXKNOTS
#undef SPECSIZE_BYTES
