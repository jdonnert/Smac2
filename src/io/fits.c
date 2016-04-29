#include "../proto.h"
#include "../globals.h"

#include <fitsio.h>

static int fits_write_header(fitsfile *, int);

#if defined(O_FITS) || defined (O_FITS_COMPRESSED)
/* Output image as fits file */
void write_output()
{
	start_timing(CPU_OUTPUT);
	
	const int npix2 = Param.XYPix * Param.XYPix;

	const long naxes[2] = { Param.XYPix, Param.XYPix };
	const long naxis = 2;

	char fname[MAXLINELENGTH];

	if (Task.Rank == 0) {

		size_t nBytes = Effect.Nimage * npix2 * sizeof(double);
		
		double *fitsimage = Malloc(nBytes);

#ifdef O_FITS
		sprintf(fname, "%s", Param.Output_File);
		printf("Writing Fits image to <%s>\n", fname);
#endif // O_FITS

#ifdef O_FITS_COMPRESSED
		sprintf(fname, "%s.fz", Param.Output_File);

		printf("Writing Rice compressed fits image to <%s>\n", fname);
#endif // O_FITS_COMPRESSED

		if (Param.NoClobber == 0)
			remove(fname);

		fitsfile *fptr = NULL;

		int status = 0;

		fits_create_file(&fptr, fname, &status);

		Assert(status == 0, "Can't open file for writing");

#ifdef O_FITS_COMPRESSED
		fits_set_compression_type(fptr, RICE_1, &status);
		fits_set_quantize_level(fptr, 4, &status);	/* accurate enough */
#endif				// O_FITS_COMPRESSED

		/* Write multiple file extensions */
		for (int img = 0; img < Effect.Nimage; img++) {

			for (int i = 0; i < naxes[0]; i++) {

				for (int j = 0; j < naxes[1]; j++) {

					size_t idx = img * npix2 + j * naxes[0] + i;

					fitsimage[i * naxes[1] + j] = Image[idx];
				}
			}

			fits_create_img(fptr, DOUBLE_IMG, naxis, naxes, &status);

			fits_write_header(fptr, img);

			fits_write_img(fptr, TDOUBLE, 1, npix2,
				       fitsimage, &status);

		}

		fits_close_file(fptr, &status);

		free(fitsimage);
	} // if rank == 0

	MPI_Barrier(MPI_COMM_WORLD);

	stop_timing(CPU_OUTPUT);

	return;
}
#endif

#ifdef O_FITS_COMPRESSED_CUBE

static double *cube = NULL;

void move_image_to_cube(const int k)
{
	const int npix2 = Param.XYPix * Param.XYPix;
	const int ncube = Param.NCube;

	rprintf("Move Image into cube \n");

	size_t nBytes = Effect.Nimage * npix2 * ncube * sizeof(*Image);

	if (Task.Rank == 0) {

		if (cube == NULL)
			cube = Malloc(nBytes);

		for (int img = 0; img <  Effect.Nimage; img++) {

			nBytes = npix2 * sizeof(*Image);

			size_t src = img * npix2;
			size_t dest = img * npix2 * ncube + k*npix2;

			memcpy(&cube[dest], &Image[src], nBytes);
		}
	}

	Free(Image);  // prepare for next iteration, prevent memory leak
	Free(Weight_Image);
	
	Image = NULL; Weight_Image = NULL;

	Free(tree); 
	
	tree = NULL;
	
	if (Task.Npart[0] != 0)
		Free(Gas);

	if (Task.PartTotal != 0)
		Free(P);

	P = Gas = NULL;

	Task.PartTotal = 0;
	for (int i = 0; i < N_part_types; i++)
		Task.Npart[i] = 0;

	Param.Center[0] /= Comv2phys.Length; // prepare for next time
	Param.Center[1] /= Comv2phys.Length;
	Param.Center[2] /= Comv2phys.Length;

	Param.XYSize /= Comv2phys.Length;
	Param.ZDepth /= Comv2phys.Length;

	return ;
}

void write_output()
{
	char fname[MAXLINELENGTH];

	start_timing(CPU_OUTPUT);
	
	const int npix2 = Param.XYPix * Param.XYPix;
	const int ncube = Param.NCube;

	long naxes[3] = { Param.XYPix, Param.XYPix, ncube };
	long naxis = 3;

	if (Task.Rank == 0) {

		size_t nBytes = Effect.Nimage * npix2 * ncube * sizeof(double);
		
		double *fitsimage = Malloc(nBytes);

		sprintf(fname, "%s.fz", Param.Output_File);
		printf("Writing Rice compressed fits cubes to <%s>\n", fname);

		if (Param.NoClobber == 0)
			remove(fname);

		fitsfile *fptr = NULL;

		int status = 0;

		fits_create_file(&fptr, fname, &status);

		Assert(status == 0, "Can't open file for writing");

		//fits_set_compression_type(fptr, RICE_1, &status);
		//fits_set_quantize_level(fptr, 4, &status);	// accurate enough 

		/* Write multiple file extensions */
		for (int img = 0; img < Effect.Nimage; img++) {

			for (int k = 0; k < naxes[2]; k++) {

				for (int j = 0; j < naxes[1]; j++) {
					
					for (int i = 0; i < naxes[0]; i++) {

						size_t idx = img * npix2 * ncube 
							+ k * naxes[0] * naxes[1] + j * naxes[0] + i;

						size_t fits_idx = k*naxes[1]*naxes[0] + i*naxes[0] + j;

						fitsimage[fits_idx] = cube[idx];
				
					}
				}
			}

			fits_create_img(fptr, DOUBLE_IMG, naxis, naxes, &status);

			fits_write_header(fptr, img);

			fits_write_img(fptr, TDOUBLE, 1, npix2*ncube, fitsimage, &status);

		}

		fits_close_file(fptr, &status);

		free(fitsimage);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	stop_timing(CPU_OUTPUT);

	return;
}


#endif

#define LASTPARAMETERID -1111
#define REAL 1
#define STRING 2
#define INT 3

int fits_write_header(fitsfile * fptr, int img)
{
	int status = 0, nt;

	fits_write_key(fptr, TSTRING, " ", " ", NULL, &status);
	fits_write_key(fptr, TSTRING, "STRING", "P-Smac II",
		       "(c) 2014 J.Donnert; donnert@ira.inaf.it",
		       &status);
	fits_write_key(fptr, TSTRING, "STRING", "Donnert & Brunetti 2014",
		       "Please cite this paper when publishing results",
		       &status);
	fits_write_key(fptr, TSTRING, "STRING", VERSION, "Version", &status);
	fits_write_key(fptr, TINT, "INT", &Task.NTask,
		       "Number of MPI Tasks", &status);
	fits_write_key(fptr, TSTRING, "STRING", "COMPILETIMEOPTIONS",
		       "COMPILETIMESETTINGS", &status);
	fits_write_key(fptr, TSTRING, "STRING", "Module Name",
		       Effect.Name, &status);
	fits_write_key(fptr, TSTRING, "STRING", "Unit", Effect.Unit, &status);
	fits_write_key(fptr, TSTRING, "STRING", "Description",
		       Effect.Descr[img], &status);
	fits_write_key(fptr, TSTRING, " ", " ", NULL, &status);

	fits_write_key(fptr, TDOUBLE, "DOUBLE", &Snap.Redshift,
		       "Redshift", &status);
	fits_write_key(fptr, TSTRING, " ", " ", NULL, &status);

	fits_write_key(fptr, TSTRING, "      ", "Parameter File Dump :", NULL,
		       &status);
	for (nt = 0; nt < 100; nt++) {
		if (id[nt] == LASTPARAMETERID)
			break;
		switch (id[nt]) {
		case INT:
			fits_write_key
			    (fptr, TINT, "INT", (int *)(addr[nt]),
			     tag[nt], &status);
			break;
		case REAL:
			fits_write_key
			    (fptr, TDOUBLE, "DOUBLE", (double *)(addr[nt]),
			     tag[nt], &status);
			break;
		case STRING:
			fits_write_key
			    (fptr, TSTRING, "STRING", (char *)(addr[nt]),
			     tag[nt], &status);
			break;
		}
	}

	fits_write_chksum(fptr, &status);

	return (status);
}

#undef REAL
#undef STRING
#undef INT

