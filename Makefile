#/*
# * PATCHRM
# * a non-Parameteric ApproaCH to constrain tthe transfer function in Reverberation Mapping
# * 
# * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
# *
# * Reference: Li, Y.-R. et al. 2016, ApJ
# *
# */

SHELL=/bin/bash

CC       = gcc -O2 -Wall
OPTIMIZE = 
#OPTIMIZE += -DJAVELINE
#OPTIMIZE += -DTOPHAT

#---------target system
#SYSTEM="Darwin"
SYSTEM="Linux"
#SYSTEM="Cluster"

ifeq ($(SYSTEM), "Darwin")
#GSL_INCL    = -I/opt/local/include
#GSL_LIBS    = -L/opt/local/lib/gsl/lib
GSL_INCL    = $(shell pkg-config --cflags gsl) 
GSL_LIBS    = $(shell pkg-config --libs gsl) 
LAPACK_INCL = -I /usr/local/share/lapack/include -I/opt/local/include
LAPACK_LIBS = -L /usr/local/share/lapack/lib -llapacke
MPFIT_LIBS = -lmpfit
MPFIT_INCS = -I /usr/local/include
OPTIMIZE    = -O2 
#-Wall
endif

ifeq ($(SYSTEM), "Linux")
NCORE      :=$(grep -c ^processor /proc/cpuinfo)
GSL_INCL    = $(shell pkg-config --cflags gsl) 
GSL_LIBS    = $(shell pkg-config --libs gsl) 
#LAPACK_INCL = -I /usr/local/share/lapack/include
#LAPACK_LIBS = /usr/local/share/lapack/lib/liblapacke.a -llapack -L/usr/lib64/atlas -lcblas 
LAPACK_INCL = -I/usr/include/lapacke
LAPACK_LIBS = -L/usr/lib64 -llapacke -llapack -lblas
MPFIT_LIBS = -lmpfit
MPFIT_INCS = -I /usr/local/include
endif

ifeq ($(SYSTEM), "Cluster")
GSL_INCL = -I/mbh/mbhd01/soft/gsl/include
GSL_LIBS = -L/mbh/mbhd01/soft/gsl/lib  -lgsl -lgslcblas -lm
LAPACK_INCL = -I/mbh/mbhd01/user/liyanrong/soft/lapack/include
LAPACK_LIBS = -L/mbh/mbhd01/user/liyanrong/soft/lapack/lib -llapacke -llapack -lblas -lgfortran
MPFIT_LIBS = -L/mbh/mbhd01/user/liyanrong/soft/mpfit -lmpfit
MPFIT_INCS = -I /mbh/mbhd01/user/liyanrong/soft/mpfit
endif

OPTIONS  = $(OPTIMIZE)
CFLAGS   = $(OPTIONS) $(GSL_INCL) $(LAPACK_INCL) $(MPFIT_INCS)
LIBS     = $(GSL_LIBS) $(LAPACK_LIBS) -lm $(MPFIT_LIBS)

EXEC     = rmpatch
SRC      = src/
INCL     = Makefile $(SRC)/allvars.h $(SRC)/proto.h

OBJS  = $(SRC)/main.o $(SRC)/readparam.o $(SRC)/allvars.o $(SRC)/error.o $(SRC)/run.o \
        $(SRC)/init.o $(SRC)/mathfunc.o $(SRC)/mcmc_con.o $(SRC)/mcmc_stats.o         \
        $(SRC)/mcmc_conline.o $(SRC)/mcmc.o $(SRC)/transferfunc.o $(SRC)/lineconv.o   \
        $(SRC)/sim.o $(SRC)/test.o

$(EXEC): $(OBJS)
	cd $(SRC)
	$(CC) $(OPTIMIZE) $(OBJS) $(LIBS) -o $@

$(OBJS): $(INCL)


clean:
	rm $(SRC)/*.o $(EXEC)
