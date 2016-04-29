#include "proto.h"
#include "globals.h"
#include "effects/effects.h"

static inline float kernel_spline_3d(float r, float h);
static inline float kernel_quintic_3d(float r, float h);
static inline float kernel_WC6_3d(float, float);
static inline float kernel_WC4_3d(float r, float h);

static void find_position_periodic(const int, const int, float *, float *, 
		float *);

static inline void preprocess(int, int, float);
static void reduce_image();
static void apply_weights();

/* 
 * Compute overlap of Kernel with pixels using Gather approximation. As the
 * overlap is only approximated, the total distributed weight will not add up 
 * to one, but we need to distribute all the particle emission. Hence we have 
 * to loop over the weights first to compute that number <1 ! Only then the 
 * particle can be distributed with kernel weights, which is straight forward
 * given there are many pixels in the kernel support radius. Following Dolag,
 * Hansen et al. 2005, particles are oversampled, when the kernel support 
 * radius is larger than a certain scale set at compile time 
 * (hsml > SAMPLING_SCALE). Then we project onto the centers of the pixels.
 * When hsml <= SAMPLING_SCALE, the particle is undersampled 
 * and we project according to the overlapping areas. Undersampled 
 * particles can overlap with only one pixel center or none and with up to 8 
 * cells. This way the effect of more pixel area for smaller hsml is minimized
 * and we a have a consistent way of handling anti aliasing.
 */

void project()
{
	start_timing(CPU_PROJECT);

	const float len2pix = Param.XYPix / Param.XYSize; // convert length to pix

	const int npix = Param.XYPix;
	const int npix2 = Param.XYPix * Param.XYPix;
	const int nImg = Effect.Nimage;

	if (prep_func_ptr != NULL)	// Do something before LoS integration 
		(*prep_func_ptr) ();

	const int counter_const = imax(1, floor(Task.PartTotal/10.0/Omp.NThreads));

	rprintf("Integrating Line of Sight ... \n");

	#pragma omp parallel 
	{

	double * restrict image = image_alloc(npix, npix, nImg);
	double * restrict wimage = image_alloc(npix, npix, 1);
	double * restrict wk = Malloc(npix2 * sizeof(*wk));

	#pragma omp for schedule(dynamic, Task.PartTotal/Omp.NThreads/64) 
	for (int ipart = 0; ipart < Task.PartTotal; ipart++) {	// in [pix]
        
		if ((ipart % counter_const) == 0)
			rprintf(".");

		double effect[MAXIMAGES] = { 0 };
	
		(*effect_func_ptr) (ipart, effect); //  "emission" of the particle
		
		double part_weight = (*weight_func_ptr) (ipart);
		
		float hsml = P[ipart].Hsml;	// now distribute in pix units

		float px = P[ipart].Pos[0];
		float py = P[ipart].Pos[1];
		float pz = P[ipart].Pos[2];

#ifdef PERIODIC
		for (int k = 0; k < 8; k++) {	// Mirror Particle 

		find_position_periodic(ipart, k, &px, &py, &pz);

		if (px - hsml > halfXYSize || px + hsml < -halfXYSize
	  	 || py - hsml > halfXYSize || py + hsml < -halfXYSize
	 	 || pz - hsml > halfXYSize || pz + hsml < -halfXYSize)
			continue; // test for out of image 
#endif // PERIODIC

		px *= len2pix;
		py *= len2pix;
		hsml *= len2pix;

		float area = pi * hsml * hsml;

		float mass = P[ipart].Mass;
		
		float rho = P[ipart].Rho / p3(len2pix);

		float dz = mass / rho / area;

		float x = px + 0.5 * npix;	// [pix] 
		float y = py + 0.5 * npix;

		int iMin = fmax(floor(x - hsml), 0); // particle extend on image [pix]
		int iMax = fmin(floor(x + hsml) + 1, npix);

		int jMin = fmax(floor(y - hsml), 0);
		int jMax = fmin(floor(y + hsml) + 1, npix);

		bool is_undersampled = false; // distribute via area overlap

		if (hsml <= SAMPLING_SCALE)
			is_undersampled = true;

		double distr_weight = 0, distr_area = 0; 
		int n_distr_pix = 0;

		for (int i = iMin; i < iMax; i++) { // find distributed area & weight

			for (int j = jMin; j < jMax; j++) {

				size_t idx = i * npix + j;

				float dx = fmin(x + hsml, i + 1) - fmax(x - hsml, i);
				float dy = fmin(y + hsml, j + 1) - fmax(y - hsml, j);

				distr_area += dx * dy;
	
				if (is_undersampled) { // no kernel 

					wk[idx] = dx * dy;
	
					distr_weight += wk[idx];

					n_distr_pix++;

					continue;
				}

				/* kernel weight via distance to pixel _centre_ */

				float x_dist = x - i - 0.5;
				float y_dist = y - j - 0.5;
	
				float r2_dist = p2(x_dist) + p2(y_dist);

				if (r2_dist >= p2(hsml)) {
	
					wk[idx] = 0;
	
					continue;
				}

				float r = sqrt(r2_dist);

#if defined(WK_CUBIC)
				wk[idx] = kernel_spline_3d(r, hsml);
#elif defined(WK_QUINTIC)
				wk[idx] = kernel_quintic_3d(r,hsml);
#elif defined(WK_WC6)
				wk[idx] = kernel_WC6_3d(r, hsml);
#elif defined(WK_WC4)
				wk[idx] = kernel_WC4_3d(r, hsml);
#endif // WK_

				wk[idx] *= dx * dy;

				distr_weight += wk[idx];

				n_distr_pix++;

			} // j
		} // i

		float weight_per_pix = n_distr_pix / distr_weight;

		float kernel_norm = area / n_distr_pix;

		float area_norm = kernel_norm * weight_per_pix * part_weight * dz;
			
		for (int iImg = 0; iImg < nImg; iImg++) { // distribute
			
			int iOffset = npix2 * iImg;

			for (int i = iMin; i < iMax; i++) {

				for (int j = jMin; j < jMax; j++) {
	
					size_t idx = i * npix + j;

					double pix_weight = area_norm * wk[idx];

					preprocess(ipart, idx, wk[idx]); // rarely needed

					image[iOffset + idx] += effect[iImg] * pix_weight;

					if (iImg == 0)
						wimage[idx] += pix_weight;
				} // j
			} // i
		} // iImg

#ifdef PERIODIC
		} // k
#endif // PERIODIC

	} // for (ipart)

	#pragma omp critical // reduce image
	{ 

	for (int i = 0; i < npix2 * nImg; i++)
	 	Image[i] += image[i];

	for (int i = 0; i < npix2; i++)
		Weight_Image[i] += wimage[i];

	}  			

	Free(wk);
	Free(image);
	Free(wimage);
	
	} // omp parallel

	MPI_Barrier(MPI_COMM_WORLD);

	stop_timing(CPU_PROJECT);

	rprintf(" done\n");

	reduce_image();

	apply_weights();

	if (post_func_ptr != NULL)	// Apply final operations on images
		(*post_func_ptr) ();

	return;
}

/* 
 * Hide ugly PERIODIC k loop. will be optimised out usually 
 */

static void find_position_periodic(const int ipart, const int k, 
		float *x, float *y, float *z)
{
	const float boxsize = Snap.Boxsize;

	*x = !(k & 0x1) ? P[ipart].Pos[0] : (P[ipart].Pos[0] > 0 ? 
			P[ipart].Pos[0] - boxsize : P[ipart].Pos[0] + boxsize);
	
	*y = !(k & 0x2) ? P[ipart].Pos[1] : (P[ipart].Pos[1] > 0 ? 
			P[ipart].Pos[1] - boxsize : P[ipart].Pos[1] + boxsize);

	*z = !(k & 0x4) ? P[ipart].Pos[2] : (P[ipart].Pos[2] > 0 ? 
			P[ipart].Pos[2] - boxsize : P[ipart].Pos[2] + boxsize);

	return;
}


/* 
 * Collect all preprocessing #ifdef 
 */

static inline void preprocess(int ipart, int idx, float kernel_weight)
{

#ifdef INTRINSIC_RM
	faraday_rotate_pix(ipart, idx, kernel_weight);
#endif

	return;
}

static void reduce_image()
{
	if (Task.NTask == 1)
		return;

	start_timing(CPU_REDUCE);

	rprintf("Reducing Image ... "); fflush(stdout);

	int nDouble = Param.XYPix * Param.XYPix * Effect.Nimage; // projected images
	size_t nBytes = nDouble * sizeof(*Image);

	MPI_Allreduce(MPI_IN_PLACE, Image, nDouble, MPI_DOUBLE, MPI_SUM,
		      MPI_COMM_WORLD);

	nDouble = Param.XYPix * Param.XYPix; // weight image 
	nBytes = nDouble * sizeof(*Weight_Image);

	MPI_Allreduce(MPI_IN_PLACE, Weight_Image, nDouble, MPI_DOUBLE, MPI_SUM,
		      MPI_COMM_WORLD);

	rprintf("done\n");

	stop_timing(CPU_REDUCE);

	return;
}

/* 
 * Devide by weight_image if needed 
 */

static void apply_weights()
{
	if (weight_func_ptr == Part_Weight_One
	 || weight_func_ptr == Part_Weight_Physical)
		return;

	const int npix2 = p2(Param.XYPix);

	rprintf("Applying particle weights... "); fflush(stdout);

	for (int iImg = 0; iImg < Effect.Nimage; iImg++) {

		for (int i = 0; i < npix2; i++) {
			
			if (Weight_Image[i] != 0)
				Image[iImg * npix2 + i] /= Weight_Image[i];
		}
	}

	rprintf("done\n");

	return;
}

/* 
 * kernel functions 3D
 */

static inline float kernel_spline_3d(float r, float h)
{
	r /= h;

	const float u = 1 - r;

	if (r > 0.5)
		return 8 / pi * 2 * u * u * u;  
	else
		return 8 / pi * (1 - 6 * r * r * u);
}

static inline float kernel_quintic_3d(float r, float h)
{
	r /= h;

	const double u = 1 - r;
	const double s = 2.0 / 3.0 - r;
	const double t = 1.0 / 3.0 - r;

	if (r < 1 / 3)
		return (2187 / (40 * pi)) * (u * u * u * u * u -
					     6 * s * s * s * s * s +
					     15 * t * t * t * t * t);
	else if (r < 2 / 3)
		return (2187 / (40 * pi)) * (u * u * u * u * u -
					     6 * s * s * s * s * s);
	else if (r < 1)
		return (2187 / (40 * pi)) * u * u * u * u * u;
	else
		return 0;

	return -1;
}

static inline float kernel_WC6_3d(float r, float h)
{
	r /= h;

	const float t = 1 - r;

	return 1365.0 / (64 * pi) / p3(h) * t * t * t * t * t * t * t * t * 
		(1 + 8 * r + 25*r*r + 32*r *r *r);
}

static inline float kernel_WC4_3d(float r, float h)
{
	const float u = r/h;
	const float t = (1-u);

	return 495.0/(32*pi) * t*t*t*t*t*t * (1 + 6*u + 35.0/3.0 * u *u);
}
