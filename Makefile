#/*
# * PATCHRM
# * a non-Parameteric ApproaCH to constrain tthe transfer function in Reverberation Mapping
# * 
# * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
# *
# * Reference: Li, Y.-R. et al. 2016, ApJ
# *
# */


EXEC     = rmpatch
SRC      = src/
INCL     = Makefile $(SRC)/allvars.h $(SRC)/proto.h

OBJS  = $(SRC)/main.o $(SRC)/readparam.o $(SRC)/allvars.o $(SRC)/error.o $(SRC)/run.o \
        $(SRC)/init.o

$(EXEC): $(OBJS)
	cd $(SRC)
	$(CC) $(OBJS) -o $@

$(OBJS): $(INCL)
