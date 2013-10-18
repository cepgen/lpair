INTER1  = lpair
S1   = ilpair-cms-pp.f
#
LIBRARY = -L/afs/cern.ch/sw/lcg/external/cernlib/2005/x86_64-slc5-gcc43-opt/lib -lpdflib804 -lmathlib -lgeant -ljetset74  -lpacklib -lkernlib \
          -lpawlib -lgraflib -lgrafX11 -lphtools
LIBINT	= libilpairg.a
OBJ	= evertex.o
#
###############################################################################
#
# Rule for objects:

FC = g77

CFLAGS = -g -fno-automatic -w -o

RM = /bin/rm  
RMFLAGS = -rf

# Make the executable

clean: 
	$(RM) $(RMFLAGS) ${INTER1}.exe last.kumac* paw.metafile *.log
Clean:
	$(RM) $(RMFLAGS) ${INTER1}.exe last.kumac* paw.metafile *.log *ps *.out ${INTER1}*.hbook OUT/*.* source/*.o
#
run:
	./lpair.exe < lpair-cms-mumu-pp.card > lpair-pp.out
#
# library
#
lib:
	cd ./source;pwd;g77  -fno-automatic -c *.f;/usr/bin/ar r libilpairg.a *.o;pwd;mv libilpairg.a ..;cd ..
#
# Make the executable
#
$(INTER1):
	$(FC) $(CFLAGS) $@.exe $(S1) $(LIBINT) $(OBJ) $(LIBRARY)

