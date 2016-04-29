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

ifeq ($(SYSTYPE),DARWIN) #icc
CC      	 =  mpicc
OPTIMIZE	 = -Ofast -Wall
MPI_LIBS 	 = -lmpich -L/Users/jdonnert/Dev/lib
MPI_INCL 	 = -I/Users/jdonnert/Dev/include
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

ifeq ($(SYSTYPE),mach64.ira.inaf.it) # gcc
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

EXEC	= Smac2

SRCDIR	= src

SRCFILES := ${shell find $(SRCDIR) -name \*.c -print}

ifeq (,$(wildcard $(SRCDIR)/print_settings.c)) # add if missing
SRCFILES += $(SRCDIR)/print_settings.c
endif

OBJFILES = $(SRCFILES:.c=.o)

INCLFILES := ${shell find src -name \*.h -print}
INCLFILES += Config Makefile $(SRCDIR)/config.h

CFLAGS = -fopenmp -std=c99 $(OPTIMIZE) $(CFITSIO_INCL) $(GSL_INCL) $(MPI_INCL)

LIBS = -lm -lgsl -lgslcblas -lcfitsio $(MPI_LIBS) $(GSL_LIBS) $(CFITSIO_LIBS) 

%.o : %.c
	@echo [CC] $@
	@$(CC) $(CFLAGS)  -o $@ -c $<

$(EXEC)	: settings $(OBJFILES)
	$(CC) $(CFLAGS)  $(OBJFILES)  $(LIBS) -o $(EXEC)
	@ctags $(SRCFILES) $(INCLUDEFILES)

$(OBJFILES)	: $(INCLFILES)

$(SRCDIR)/config.h : Config 
	@echo 'Config -> config.h'
	@sed '/^#/d; /^$$/d; s/^/#define /g' Config > $(SRCDIR)/config.h

$(SRCDIR)/print_settings.c : Config
	@echo '-> print_settings.c'
	@echo '#include "proto.h"' >  $(SRCDIR)/print_settings.c
	@echo '#include "globals.h"' >>  $(SRCDIR)/print_settings.c
	@echo 'void print_compile_time_settings(){' >> $(SRCDIR)/print_settings.c
	@echo 'rprintf("Compiled with : \n"' >> $(SRCDIR)/print_settings.c
	@sed '/^#/d; /^$$/d; s/^/"   /g; s/$$/ \\n"/g;' Config >> $(SRCDIR)/print_settings.c
	@echo ');}' >> $(SRCDIR)/print_settings.c

.PHONY : settings
	
settings : 
	@echo " "
	@echo 'CC = ' $(CC)
	@echo 'CFLAGS =' $(CFLAGS)
	@echo 'LDFLAGS =' $(LIBS)
	@echo 'EXEC =' $(EXEC)
	@echo " "

clean : 
	rm -f  $(OBJFILES) $(EXEC) src/config.h src/print_settings.c

