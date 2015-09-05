void read_param_file();
void setup();
void select_cosmology(int);
void read_snapshot(char *);
long guess_snapnum(char *);
void find_DM_densities();
void write_output();
void remove_unsused_particles();

void select_effect_module(int choice);

void domain_decomposition();
double Metric_Hsml2(int);
double Metric_One(int);

void project();

#ifdef O_FITS_COMPRESSED_CUBE
void move_image_to_cube(const int k);
#endif

/* helpers */
void print_compile_time_settings();

void *Malloc_info(size_t, const char *, const char *, const int);
void *Realloc_info(void *, size_t, const char *, const char *, const int);
void Assert_Info(const char *, const char *, int, int, const char *, ...);
void Free_info(void *, const char *, const char *, const int);
double *image_alloc(size_t, size_t, size_t);
void Reallocate_P(int, int *, int);
