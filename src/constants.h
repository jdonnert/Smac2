/*mathematical constants*/
#define pi 			M_PI
#define sqrt2		M_SQRT2
#define sqrt3       1.73205077648163
#define fourpithirds (4.18879032135009765)
#define gamma3      2.678938534707748	/* Gamma(1/3) */

/*physical constants cgs gauss system*/
#define c			GSL_CONST_CGSM_SPEED_OF_LIGHT
#define e			(GSL_CONST_CGSM_ELECTRON_CHARGE*GSL_CONST_CGSM_SPEED_OF_LIGHT)
#define h_planck	GSL_CONST_CGSM_PLANCKS_CONSTANT_H
#define hbar 		GSL_CONST_CGSM_PLANCKS_CONSTANT_HBAR
#define k_B 		GSL_CONST_CGSM_BOLTZMANN
#define m_p 		GSL_CONST_CGSM_MASS_PROTON
#define m_pi0		(134.9766/c/c*1e6*eV2cgs)	/* Eidelman et al. 2004 */
#define m_pi		(139.57018/c/c*1e6*eV2cgs)
#define m_mu		(105.65837/c/c*1e6*eV2cgs)
#define m_e			GSL_CONST_CGSM_MASS_ELECTRON
#define sigma_T	    GSL_CONST_CGSM_THOMSON_CROSS_SECTION
#define sigma_pp	(32*1e-3*barn2cgs)	/*  */

/* unit conversions */
#define barn2cgs	GSL_CONST_CGSM_BARN
#define eV2cgs		GSL_CONST_CGSM_ELECTRON_VOLT
#define GeV2cgs	    (eV2cgs*1e9)
#define Msol2cgs    ((double)(1.98892e33))
#define kpc2cgs 	((double)(3.08568025e21))
#define pc2cgs 	    ((double)(3.08568025e18))
#define pix2cm		(Param.XYSize/Param.XYPix*kpc2cgs)
#define yr2sec		31556926
#define Gyr2sec		(yr2sec*1e9)
#define deg2rad     (pi/180.0)

/* other */
#define Tcmb		(2.728*(1+Snap.Redshift))	/* Temp of the CMB [K] */
#define Bcmb		(3.24516e-6*(1+Snap.Redshift))	/* MF of the CMB [G] */
#define H_frac		0.76	/* Hydrogen fraction */
#define He_frac	    (1.0-H_frac)	/* Helium fraction */
#define u_mol		(4.0/(5.0*H_frac+3.0))	/* Mean mol. weight in hydr. mass */
#define n2ne 		((H_frac+0.5*He_frac)/(2.0*H_frac+0.75*He_frac))
#define yHelium	    (He_frac / (4.0 *H_frac))
#define mean_mol_weight ((1.0+4.0*yHelium)/(1.0+3.0*yHelium+1.0))
#define adiabatic_index (5./3.)
