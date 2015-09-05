#ifdef I_GADGET2		/*Additional Stuff for GADGET-2 Input */

int swap;

enum iofields {
	IO_POS,
	IO_VEL,
	IO_ID,
	IO_U,
	IO_RHO,
	IO_HSML,
	IO_VOL,			/* AREPO only */
	IO_BFLD,
	IO_DIVB,
	IO_MASS,
	IO_VELT,
	IO_VRMS,
	IO_VBULK,
	IO_VDIV,
	IO_VROT,
	IO_TNGB,
	IO_DPP,
	IO_VTAN,
	IO_VRAD,
	IO_nHII,
	IO_nHeII,
	IO_nHeIII,
	IO_MACH,
	IO_VORT,
	IO_ABVC,
	IO_LASTENTRY		/* Keep this entry at the end for termination */
};

void swap_Nbyte(char *, int, int);
size_t my_fread(void *, size_t, size_t, FILE *);
int find_block(FILE *, char *);
int find_files(char *);
void read_file(char *, int, int);
int read_gadget_block(void *, char *, FILE *, size_t);
int read_gadget_head(FILE *);
void set_block_prop(enum iofields);
void empty_comm_buffer(enum iofields, void *, size_t, size_t);

struct Gadget_head {		/*Input Header */
	long long Npart[N_part_types];
	double Mass[N_part_types];
	double Time;
	double Redshift;
	long long Nall[N_part_types];
	long NumFiles;
	double BoxSize;
	double Omega0;
	double OmegaLambda;
	double HubbleParam;
	int FlagSfr;
	int FlagFeedback;
	int FlagCooling;
	int FlagAge;
	int FlagMetals;
} Header;

/* Block specifications */
struct blockdef {
	char *Label;
	char *Name;
	void *DataPtr;
	long long Npart[N_part_types];
	long long Ntot;
	int Val_per_element;
	size_t Data_type;
	size_t Bytes_per_element;
	double Rmv_comoving;
} Block;
#endif
