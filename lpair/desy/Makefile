F_FILES = $(wildcard source/*.f)
LIB_FILES = $(patsubst source/%.f,obj/%.o,$(F_FILES))
LIBRARY = -L/afs/cern.ch/sw/lcg/external/cernlib/2005/x86_64-slc5-gcc43-opt/lib \
          -lpdflib804 -lmathlib -lgeant -ljetset74 -lpacklib -lkernlib \
          -lpawlib -lgraflib -lgrafX11 -lphtools
OBJ_FILES = ilpair-cms-pp.o evertex.o $(LIB_FILES)
VPATH = source/

###############################################################################

FC = g77
FFLAGS = -g -fno-automatic -w

RM = /bin/rm
RMFLAGS = -rf

.PHONY: all
all: lpair

#
# Make the executable
#
lpair: $(OBJ_FILES)
	$(FC) $(FFLAGS) -o $@ $^ $(LIBRARY)

nice:
	$(RM) $(RMFLAGS) *.o obj/*.o

clean: nice
	$(RM) $(RMFLAGS) lpair
#
# Make the objects
#
obj/%.o: %.f
	$(FC) -c $(FFLAGS) $< -o $@
