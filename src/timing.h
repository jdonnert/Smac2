/* PROFILING
 * We define marks in this enum,
 * which are used to let the start
 * and stop timing functions keep
 * track of what is measured.
 */

enum TimeMarks {
	CPU_READIN,
	CPU_SETUP,
	CPU_DOMAINDECOMP,
	CPU_TREE,
	CPU_NGBFIND,
	CPU_PROJECT,
	CPU_EFFECT,
	CPU_REDUCE,
	CPU_OUTPUT,
	CPU_LASTENTRY		/* Keep this entry at the end */
} mark;

struct Timings {
	double Zero;		/* Absolute Zeropoint */
	double Total[CPU_LASTENTRY];	/* Total Time for mark */
	double Start[CPU_LASTENTRY];	/* Start time for mark */
} Cpu;

void init_timing();
void finish_timing();
void start_timing(enum TimeMarks mark);
void stop_timing(enum TimeMarks mark);
double get_current_time();
