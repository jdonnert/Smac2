/* Cosmologies
 * We set Cosmo to one of the others 
 * */

#include "proto.h"
#include "globals.h"

struct cosmology Cosmo;

/* Spergel et al 2006 */
struct cosmology WMAP3 = { "WMAP3", 0.732, 0.0444, 0.266, 0.732,
	-1.0, 0.079, 0.94E-26, 0.772, 13.73E9
};

/* Old Smac1 settings */
struct cosmology SMAC1 = { "SMAC1", 0.7, 0, 0.3, 0.7, 0, 0, 0, 0, 1.34670e+10 };

/* Non Comoving projections */
struct cosmology null = { "NULL", 1, 0, 0, 0, 0, 0, 1, 0, 0 };
