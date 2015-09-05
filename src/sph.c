#include "globals.h"
#include "tree.h"

static inline float sph_kernel_WC6(const float r, const float h);
static inline float sph_kernel_derivative_WC6(const float r, const float h);

static bool find_hsml(const int ipart, const int *ngblist, const int ngbcnt,
		      float *dRhodHsml_out, float *hsml_out, float *rho_out);

static void calc_derived_kernel_quantities();

extern void find_SPH_densities()
{
	Build_tree();

	start_timing(CPU_NGBFIND);

	rprintf("Finding Hsml, ");

	const int counter_const = floor(Task.PartTotal/10.0/Omp.NThreads);

	#pragma omp parallel for shared(P, tree) \
		schedule(dynamic, Task.PartTotal/Omp.NThreads/64)
	for (size_t ipart = 0; ipart < Task.PartTotal; ipart++) {
		
		if (ipart % counter_const == 0)
			rprintf(".");

		float hsml = P[ipart].Hsml;

        if (hsml == 0)
            hsml = Guess_hsml(ipart, DESNNGB); // always too large

        float dRhodHsml = 0;
        float rho = 0;

        for (;;) {

            int ngblist[NGBMAX] = { 0 };

            int ngbcnt = Find_ngb_tree(ipart, hsml, ngblist); 

			if (ngbcnt > NGBMAX) {
			
				hsml /= 1.24;

				continue;
			}

            bool part_done = find_hsml(ipart, ngblist, ngbcnt, &dRhodHsml, 
                    &hsml, &rho); 
            
            if (part_done)
                break;
        }

        float varHsmlFac = 1.0 / ( 1 + hsml/(3*rho)* dRhodHsml );

        P[ipart].Hsml = hsml;
        P[ipart].Rho = rho;
        P[ipart].VarHsmlFac = varHsmlFac;

	}

	rprintf("done \n\n");

	calc_derived_kernel_quantities();

	stop_timing(CPU_NGBFIND);

	return;
}

/* 
 * solve SPH continuity eq via bisection and tree search 
 */

static bool find_hsml(const int ipart, const int *ngblist, const int ngbcnt,
	  float *dRhodHsml_out, float *hsml_out, float *rho_out)
{
	const double boxhalf = 0.5 * Snap.Boxsize;
    const double boxsize = Snap.Boxsize;
    const double pos_i[3] = {P[ipart].Pos[0],P[ipart].Pos[1],P[ipart].Pos[2]};
	
	const double mpart = fmax(Snap.Masstab[1], P[ipart].Mass);

    double upper = *hsml_out * sqrt3;
    double lower = 0;
    
    double hsml = upper; //lower + 0.5*(upper-lower); 
    double rho = 0, dRhodHsml = 0;

    int it = 0;

    bool part_done = false;

    for (;;) {  

        double wkNgb = 0; // kernel weight number of neighbours

        rho = dRhodHsml = 0;
        
        it++;

        for (int i = 0; i < ngbcnt; i++) {

	    	int jpart = ngblist[i];	

            double dx = pos_i[0] - P[jpart].Pos[0];
	    	double dy = pos_i[1] - P[jpart].Pos[1];
    		double dz = pos_i[2] - P[jpart].Pos[2];
			
		    if (dx > boxhalf)	// find closest image 
	    		dx -= boxsize;

    		if (dx < -boxhalf)
		    	dx += boxsize;

	    	if (dy > boxhalf)
    			dy -= boxsize;

		    if (dy < -boxhalf)
	    		dy += boxsize;

    		if (dz > boxhalf)
			    dz -= boxsize;

		    if (dz < -boxhalf)
	    		dz += boxsize;

            double r2 = dx*dx + dy*dy + dz*dz;

    		if (r2 > p2(hsml)) 
                continue ;
                    
    		double r = sqrt(r2);

		    double wk = sph_kernel_WC6(r, hsml);
		    double dwk = sph_kernel_derivative_WC6(r, hsml);

            wkNgb += 4*pi/3*wk*p3(hsml);

            rho += mpart * wk;

            dRhodHsml += -mpart * ( 3/hsml*wk + r/hsml * dwk );
        }

        if (it > 60) {

            printf("<%d> Hsml iteration reached 30: ipart=%d up=%g lo=%g h=%g "
                    "Ngbcnt=%d wkNgb=%g DesNngb=%d \n"
                    "x=%g y=%g z=%g upstart=%g, deltas= %g %g %g\n", 
                    Omp.ThreadID, ipart, upper, lower, hsml, ngbcnt, wkNgb, 
                    DESNNGB, pos_i[0], pos_i[1], pos_i[2], *hsml_out, 
                    upper-lower, upper-hsml, lower-hsml);
            fflush(stdout);
            
            break;
        }
        
        if (fabs(wkNgb-DESNNGB) < NGBDEV) {

            part_done = true;

            break;
        }

        if (fabs(upper-lower) < 1e-4) { // find more neighbours ! 
            
            hsml *= 1.26; // double volume
            
            //printf("Repeat %d %d %g %g \n", ipart, ngbcnt, hsml, wkNgb);

            break;
        }

        if (wkNgb > DESNNGB) // bisection
            upper = hsml;

        if (wkNgb < DESNNGB) 
            lower = hsml;
        
        hsml = pow(0.5 * ( p3(lower) + p3(upper) ), 1.0/3.0);
    }
    
    *hsml_out = (float) hsml;

    if (part_done) {
    
        *dRhodHsml_out = (float) dRhodHsml;

		*rho_out = (float) rho;

        double bias_corr = -0.0116 * pow(DESNNGB*0.01, -2.236) 
            * mpart * sph_kernel_WC6(0, hsml); // WC6 (Dehnen+ 12)
    
        *rho_out = (float) (rho + bias_corr );

		if (P[ipart].Type == 0) {
#ifdef VSPH
			Gas[ipart].VelSPH[0] = vSPH[0] / total_weight;
			Gas[ipart].VelSPH[1] = vSPH[1] / total_weight;
			Gas[ipart].VelSPH[2] = vSPH[2] / total_weight;
#endif
		}
	}

	return part_done;
}

static void calc_derived_kernel_quantities()
{

#if defined(VSPH)

	rprintf("Computing derived kernel quantities ...");

#pragma omp parallel for shared(P) \
	schedule(dynamic, Task.PartTotal/Omp.NThreads/64)
	for (int ipart = 0; ipart < Task.PartTotal; ipart++) {

		float hsml = P[ipart].Hsml; // double volume

		int ngblist[NGBMAX] = { 0 };

		int ngbcnt = Find_ngb_tree(ipart, hsml, ngblist);

		float vBulk[3] = { 0 };

		for (int i = 0; i < ngbcnt; i++) {

			int jpart = ngblist[i];

			vBulk[0] += Gas[jpart].VelSPH[0];
			vBulk[1] += Gas[jpart].VelSPH[1];
			vBulk[2] += Gas[jpart].VelSPH[2];
		}

		vBulk[0] /= ngbcnt;
		vBulk[1] /= ngbcnt;
		vBulk[2] /= ngbcnt;

		float vTurb = 0;

		for (int i = 0; i < ngbcnt; i++) {

			int jpart = ngblist[i];

			vTurb += p2(Gas[jpart].VelSPH[0] - vBulk[0])
			    + p2(Gas[jpart].VelSPH[1] - vBulk[1])
			    + p2(Gas[jpart].VelSPH[2] - vBulk[2]);
		}

		vTurb /= ngbcnt;

		Gas[ipart].VelBulk[0] = vBulk[0];
		Gas[ipart].VelBulk[1] = vBulk[1];
		Gas[ipart].VelBulk[2] = vBulk[2];

		Gas[ipart].VelTurb = sqrt(vTurb);
	}

	rprintf(" done\n");
#endif

	return;
}

static inline float sph_kernel_WC6(const float r, const float h)
{   
	const double u= r/h;
    const double t = 1-u;

    return 1365.0/(64*pi)/p3(h) *t*t*t*t*t*t*t*t*(1+8*u + 25*u*u + 32*u*u*u);
}

static inline float sph_kernel_derivative_WC6(const float r, const float h)
{   
	const float u = r/h;
    const double t = 1-u;

    return 1365.0/(64*pi)/(h*h*h*h) * -22.0 *t*t*t*t*t*t*t*u*(16*u*u+7*u+1);
}

