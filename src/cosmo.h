struct cosmology {		/*cosmologies */
	char name[5];		/*Identifier */
	float h;		/*hubble parameter */
	float Omega_b;		/*baryon overdensity */
	float Omega_m;		/*matter overdensity */
	float Omega_l;		/*dark energy overdensity */
	float w;		/*Equation of state parameter */
	float tau;		/*optical depth of reionization */
	double rho0;		/*critical density */
	float sigma8;		/*galaxy fluctuation amplitude */
	float t0;		/*age of the universe */
} Cosmo, WMAP3, SMAC1, null;	/* Cosmo is set to one of the other */

double NoComov(float);
