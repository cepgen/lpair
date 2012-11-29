####################
# LINUX (Intel PC) #
####################
#
#OPT = -c -C  -Wunused -Wuninitialized -malign-double -m486 -O
#OPT2 = -c -m486
#
OPT = -Wall -Wno-unused -Wno-globals \
      -O -funroll-loops -fomit-frame-pointer \
      -fno-globals \
      -fno-backslash -fno-automatic -finit-local-zero \
      -fno-typeless-boz  -fugly-complex \
      -funderscoring
#
FC = g77
LO = g77
#
# Libraries:
#
#CERNLIB = -L /usr/lib/cernlib/2006/lib \
CERNLIB = -L/nfs/soft/cernlib/2006b_x86_64-slc5-gcc41-opt/lib/ \
          -lpdflib804 -lmathlib -lpacklib -lkernlib \
          -lpawlib -lgraflib -lgrafX11 -lphtools
#-ljetset74
ATLSOFT = ../lib
GENZLIB = $(ATLSOFT)/libgenz.so
#
# Object files:
#
OBJECTS = eemumu.o myfileini.o \
          save.o save2.o vegas.o setgen.o genera.o myaccept.o \
          treat.o ranf.o dangle.o \
          genzini.o genzfil.o genzend.o
#          pawini.o pawfil1.o pawfil2.o pawend.o \

#.f.o:
.o:
	$(FC) $(OPT) -c $<
#genz:
#	$(FC) $(OPT) -shared -o libgenz.so ../lib/*.o $(CERNLIB) -m32

mylpair: $(OBJECTS)
	@echo "   linking program ..."
	$(LO) $(OPT) -o mylpair.exe $(OBJECTS) $(CERNLIB) $(ATLSOFT)/*.o
	@echo "   run program with 'mylpair.job' "

clean:
	rm -f *.exe *.o
