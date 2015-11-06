SHELL = /bin/bash	# This should always be present in a Makefile

## Target Computer ##
ifndef SYSTYPE
	SYSTYPE	:= $(shell hostname)
endif

#Std systype
CC		 = mpicc
OPTIMIZE = -Wall -g -O3
MPI_INCL = 
MPI_LIBS = -lmpi 
GSL_INCL =
GSL_LIBS = 
CFITSIO_INCL = 
CFITSIO_LIBS = 
HDF5_INCL =
HDF5_LIBS =

ifeq ($(SYSTYPE),DARWIN)
CC      	 =  mpicc
OPTIMIZE	 = -Ofast -Wall -mtune=native -march=corei7
MPI_LIBS 	 = -lmpich -L/Users/julius/Devel/lib
MPI_INCL 	 = -I/Users/julius/Devel/include
GSL_INCL 	 =  
GSL_LIBS	 = 
CFITSIO_INCL = 
CFITSIO_LIBS = 
HDF5_INCL 	 =
HDF5_LIBS    =
endif

ifeq ($(SYSTYPE),getorin.ira.inaf.it)
CC      	 =  mpicc
OPTIMIZE	 = -O3 -Wall -g -lmpich #-openmp-report=2 -finline -finline-functions -funroll-loops -xhost  -mkl -use-intel-optimized-headers -ipo4 -fast-transcendentals 
MPI_LIBS 	 = -L/homes/donnert/Libs/lib -L/homes/donnert/Libs/opt/intel/lib
MPI_INCL 	 = -I/homes/donnert/Libs/include -I/homes/donnert/Libs/opt/intel/include
GSL_INCL 	 =  
GSL_LIBS	 = 
CFITSIO_INCL = 
CFITSIO_LIBS = 
HDF5_INCL 	 =
HDF5_LIBS    =
endif

ifeq ($(SYSTYPE),mach64.ira.inaf.it)
CC      	 =  mpicc
OPTIMIZE	 =  -g -O2  -march=bdver1 -mtune=native -mprefer-avx128 -mieee-fp -flto 
MPI_LIBS 	 =  -lmpich -L/homes/donnert/Libs/lib 
MPI_INCL 	 = -I/homes/donnert/Libs/include 
GSL_INCL 	 =  
GSL_LIBS	 = 
CFITSIO_INCL = 
CFITSIO_LIBS = 
HDF5_INCL 	 =
HDF5_LIBS    =
endif

ifeq ($(SYSTYPE),geeuw.strw.leidenuniv.nl)
	CC           =  mpicc
	OPTIMIZE     =  -g -O3  -march=native -flto
	MPI_LIBS     =  -lmpich -L/homes/donnert/Libs/lib
	MPI_INCL     = -I/homes/donnert/Libs/include
	GSL_INCL     =
	GSL_LIBS     =
	CFITSIO_INCL =
	CFITSIO_LIBS =
	HDF5_INCL    =
	HDF5_LIBS    =
endif

ifeq ($(SYSTYPE),para33.strw.leidenuniv.nl)
	CC           =  mpicc
	OPTIMIZE     =  -g -O3  -march=native -flto
	MPI_LIBS     =  -lmpich 
	MPI_INCL     = 
	GSL_INCL     = -I/home/jdonnert/data_para/Libs/include
	GSL_LIBS     = -L/home/jdonnert/data_para/Libs/lib
	CFITSIO_INCL =
	CFITSIO_LIBS =
	HDF5_INCL    =
	HDF5_LIBS    =
endif

ifeq ($(SYSTYPE),MSI)
	CC           =  mpicc
	OPTIMIZE     =  -g -O3  -xhost -ipo4
	MPI_LIBS     =  -lmpich 
	MPI_INCL     = 
	GSL_INCL     = 
	GSL_LIBS     = 
	CFITSIO_INCL =
	CFITSIO_LIBS =
	HDF5_INCL    =
	HDF5_LIBS    =
endif

EXEC	= P-Smac2

SRCDIR	= src/

OBJFILES= aux.o setup.o main.o domain_decomp.o \
		  project_sph.o timing.o tree.o  \
		  unit.o cosmo.o healpix.o print_settings.o sph.o \
		  io/input.o \
			io/gadget.o \
			io/user.o \
			io/fits.o \
		  effects/effects.o \
		  	effects/density.o \
			effects/velocity.o \
			effects/fkp.o \
			effects/pressure.o \
			effects/temp.o \
			effects/xray.o \
			effects/bfld.o \
			effects/sz.o \
			effects/cr_gammas.o \
			effects/synchrotron.o \
			effects/analytic_synchrotron.o \
			effects/cre_secondaries.o \
			effects/cre_powerlaw.o \
			effects/cre_tabulated.o \
			effects/cre_compressed.o \
			effects/dm_annihilation.o \
			effects/rm.o \
			effects/radtransfer.o \
			effects/coulomb.o \
			effects/shocks.o


OBJS	= $(addprefix $(SRCDIR),$(OBJFILES))

INCLFILES	= globals.h  tree.h cosmo.h unit.h timing.h  config.h proto.h \
			  macro.h constants.h \
			  io/gadget.h \
			  		io/fits.h \
			  effects/effects.h \
			  ../Makefile

INCL	= $(addprefix $(SRCDIR),$(INCLFILES))

CFLAGS	= -fopenmp -std=c99 $(OPTIMIZE) $(CFITSIO_INCL) $(GSL_INCL) $(MPI_INCL)

LIBS	= -lm -lgsl -lgslcblas -lcfitsio \
		  $(MPI_LIBS) $(GSL_LIBS) $(CFITSIO_LIBS) 

$(EXEC)	: $(OBJS)
	$(CC) $(CFLAGS)  $(OBJS)  $(LIBS) -o $(EXEC)
	cd src && ctags  *.[ch]

$(OBJS)	: $(INCL)

$(SRCDIR)config.h : Config 
	sed '/^#/d; /^$$/d; s/^/#define /g' Config > $(SRCDIR)config.h

$(SRCDIR)print_settings.c : Config
	echo '#include "globals.h"' >  $(SRCDIR)print_settings.c
	echo 'void print_compile_time_settings(){' >> $(SRCDIR)print_settings.c
	echo 'rprintf("Compiled with : \n"' >> $(SRCDIR)print_settings.c
	sed '/^#/d; /^$$/d; s/^/"   /g; s/$$/ \\n"/g;' Config >>  $(SRCDIR)print_settings.c
	echo ');}' >> $(SRCDIR)print_settings.c

clean		: 
	rm -f  $(OBJS) $(EXEC) src/config.h src/print_settings.c

install	: 
	cp -i $(EXEC) ~/bin
