S1      = ilpair-cms-pp.f
INTER1  = lpair
#
LIBRARY = -L/usr/lib/cernlib/2006/lib/ -lkernlib -lpawlib -lgeant321 -lmathlib -lpawlib -lpdflib804 -lg2c
JETSET	= -L/cern/pro/lib/ -ljetset74
LIBINT	= libilpairg.a
#LIBINT	= libilpair.a
OBJ	= evertex.o
#
###############################################################################
#
# Rule for objects:

#FC = g77
FC = gfortran

CFLAGS = -g -fno-automatic -fno-second-underscore -w
CFLAGSLIB = -g -fno-automatic -w

RM = /bin/rm  
RMFLAGS = -rf

# Make the executable

clean: 
	$(RM) $(RMFLAGS) ${INTER1}.exe last.kumac* paw.metafile *.log
Clean:
	$(RM) $(RMFLAGS) ${INTER1}.exe last.kumac* paw.metafile *.log *ps *.out ${INTER1}*.hbook OUT/*.* source/*.o
#
run:
	./lpair.exe < lpair-cms-mumu-pp.card > OUT/ilapir-cms-pp.out &
#
# library
#
lib:
	cd ./source; pwd; $(FC) $(CFLAGS) -c *.f $(LIBRARY) $(JETSET); /usr/bin/ar r libilpairg.a *.o; pwd; mv libilpairg.a ..; cd ..
#
# Make the executable
#
$(INTER1):
	$(FC) $(CFLAGS) -o $@.exe $(S1) $(LIBINT) $(OBJ) $(JETSET) $(LIBRARY)
#
