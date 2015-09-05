/* WEIGHTING FUNCTIONS */
double Part_Weight_One(int);
double Part_Weight_Physical(int);
double Part_Weight_Rho(int);
double Part_Weight_Emission(int);
double Part_Weight_XrayBand(int);
double Part_Weight_Spectroscopic(int);

/* EFFECTS */
void Density(int, double *);	// 2D density  [g/cm^2/kpc] 
void NumDensity(int, double *);

void Dpp(int, double *);	// Reacceleration Coefficient 
void Tau_reacc(int, double *);	// Reacceleration Timescale 
void Tau_cooling(int, double *);	// Cooling Timescale 
void set_Dpp_factors();

void Vel(int, double *);	// Velocity 
void Vsound(int, double *);	// Speed of Sound 
void Velt(int, double *);	//RMS vel around central particle vel 
void Vrms(int, double *);	// RMS vel around mean vel in hsml 
void Velt_scaled(int, double *);	// RMS vel around central particle vel scaled 
void Vrms_scaled(int, double *);	// RMS vel around mean vel in hsml scaled 
void Vtan(int, double *);	//tangential RMS vel around mean vel in hsml 
void Vrad(int, double *);	// radial RMS vel around mean vel in hsml 
void Vbulk(int, double *);	// Mean vel in hsml 
void Vdiv(int, double *);	// vel divergence in hsml 
void Vrot(int, double *);	// Vel curl in hsml 
void Vel_perp(int, double *);	// Velocity perpendicular to LoS 
void Vorticity(int, double *);	// Vorticity 
void Kernel_weighted_velocity(int, double *);
void Kernel_weighted_bulk_velocity(int, double *);
void Kernel_weighted_velocity_dispersion(int, double *);
void Turbulent_energy(int ipart, double *);

void Xray(int, double *);	// continuum X-Ray surface brightness[erg/cm^2]

void Therm_Pressure(int, double *);	// Pressure [g/s^2/cm] 
void Magn_Pressure(int, double *);	// magnetic Pressure 
void Turb_Pressure(int, double *);	// turbulent Pressure 

void Temp(int, double *);	// Temperature [K] 
void Soundspeed(int, double *);	// Speed of sound [cm/s] 

void bfld(int, double *);	// |B| [G] 
void bfld_orth(int, double *);	// B orthogonal to Line of Sight 
void bfld_divergence(int, double *);	// div(B) 
void faraday_rotation(int, double *);	// Faraday Rotation [rad/cm^2]  
void alven_vel(int, double *);	// Alven Velocity [cm/s] 
void divB(int, double *);	// Magnetic field divergence from snap 
void plasma_beta(int, double *);	// Plasma Beta parameter 
void set_faraday_factors();

void comptonY(int, double *);	// Compton y Parameter 
void tSz(int, double *);	// Thermal SZ effect 
void kSz(int, double *);	// Kinetic SZ effect 
void set_sz_factors(bool);	// Also select DT/T or DI/I 

void crGamma_dif(int, double *);	// Gamma Ray emission from hadronic 
void crGamma_int(int, double *);	// interaction of CRp. integral & dif. 
void set_crGamma_factors();

void synchrotron(int, double *);	// Synchrotron emission from CRe 
void set_synchro_factors();
void prepare_synchrotron();
#ifdef ADD_IPOL_CHI_PI
void add_Ipol_chi_Pi();		// Add 3 derived quantities 
#endif

void aSynchrotron(int, double *);	// Analytic synchrotron emission 
void set_aSynchro_factors();

void dm_density(int, double *);	// DM density 
void dm_annihilation(int, double *);	// DM annihilation signal 

void rt_Hydrogen(int, double *);	// Ionised Hydrogen Fractioidon 

void coul_mfp(int, double *);	// Mean Free Path , Coulomb Scattering
void coul_coll(int, double *);	// hsml over mean Free Path 

void Machnumber(int, double *);	// Machnumber 
void AlphaVisc(int, double *);	// Time dependent artif. visc. param 

/* CR SPECTRA */
double cre_secondaries_brunetti(double, int);	// Brunetti  05
void set_brunetti_factors();	// Secondary CRe from Hadronic interactions 

double cre_spectrum_tabulated(double, int);	// Read CRe spectrum from file, radius dep. 
void set_tabulated_factors();

double cre_spectrum_powerlaw(double, int);	// User artifisc. CRe spectrum rel. to thermal 
void set_powerlaw_factors(int);

double cre_spectrum_compressed(double, int);
void setup_decompression();

/* MISC */
void find_central_density();	// for density scaling of CRp norm 
double Cluster_central_density;

void set_intrinsic_RM();
void faraday_rotate_pix();	// intrinsic faraday rot 
