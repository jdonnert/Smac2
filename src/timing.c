/* Timings & Benchmarking */

#include "globals.h"

struct Timings Cpu = { 0 };

const char TimeMarkDescr[CPU_LASTENTRY][20] = {
	"Readin      ",
	"Setup      ",
	"Dom. Decomp.",
	"Tree Constr.",
	"NGB Finding ",
	"Projection  ",
	"Effect      ",
	"Reduction  ",
	"Output      "
};

void init_timing()
{
	for (int i = 0; i < CPU_LASTENTRY; i++) {
		Cpu.Total[i] = 0;
		Cpu.Start[i] = -1;
	}

	Cpu.Zero = get_current_time();

	return;
}

void finish_timing()
{
	double myTime = 0, timeTotal = 0;

	rprintf("\n Timing Statistics [s] on <%d> CPUs : \n"
		"----------------------------\n", Task.NTask);

	for (int mark = 0; mark < CPU_LASTENTRY; mark++) {

		myTime = Cpu.Total[mark];

		timeTotal += Cpu.Total[mark];

		rprintf("%s	%g		\n", TimeMarkDescr[mark],
			myTime);
	}
	rprintf("----------------------------\n"
		"Total		%g		 	        \n\n",
		timeTotal);

	return;
}

void start_timing(enum TimeMarks mark)
{
	if (Cpu.Start[mark] > 0)
		fprintf(stderr, "Timing for Mark %i already running \n",
			(int)mark);
	else
		Cpu.Start[mark] = get_current_time();

	return;
}

void stop_timing(enum TimeMarks mark)
{
	if (Cpu.Start[mark] <= 0)
		fprintf(stderr, "Timing for Mark %i wasn't started \n",
			(int)mark);
	else
		Cpu.Total[mark] += get_current_time() - Cpu.Start[mark];

	Cpu.Start[mark] = 0;

	return;
}

double get_current_time()
{
	return clock();
	//return MPI_Wtime();
}
