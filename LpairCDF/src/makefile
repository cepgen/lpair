####################
# LINUX (Intel PC) #
####################
#
#OPT = -c -C  -Wunused -Wuninitialized -malign-double -m486 -O
#OPT2 = -c -m486
#
OPT = -Wall -Wno-unused \
      -O \
      -funroll-loops \
      -fomit-frame-pointer \
      -fno-backslash \
      -fno-automatic \
      -funderscoring \
      -c
#      -ff2c-library \
#      -Wno-globals \
#      -O
#      -finit-local-zero \
#      -fno-globals \
#      -fno-typeless-boz \
#      -fugly-complex  \
#      -fsecond-underscore \
      -c
OPT2 = $(OPT)
#
FC = g77
LO = g77
#FC = gfortran
#LO = gfortran
#
# Libraries:
#
CERNLIB_LOCATION = /usr/lib/cernlib/2006/lib
LD_LIBRARY_PATH = $(CERNLIB_LOCATION)
#
CERNLIB = -L$(CERNLIB_LOCATION) -lpawlib -lkernlib -lpdflib804 -lmathlib -lpacklib -lgraflib -lgrafX11 -lphtools 
#	  -L$(CERNLIB_LOCATION)libkernlib.so.2_gfortran \
#          -L$(CERNLIB_LOCATION)libpawlib.so.2_gfortran \
#          -L$(CERNLIB_LOCATION)libpawlib-lesstif \
#          -L$(CERNLIB_LOCATION)libpdflib804 \
#	  -L$(CERNLIB_LOCATION)libmathlib \
#	  -L$(CERNLIB_LOCATION)libpacklib \
#	  -L$(CERNLIB_LOCATION)libgraflib \
#	  -L$(CERNLIB_LOCATION)libgrafX11 \
#	  -L$(CERNLIB_LOCATION)libphtools \
#	  -L$(CERNLIB_LOCATION)libjetset74
ATLSOFT = ../lib
GENZLIB = -L$(ATLSOFT) -lgenz
#
# Object files:
#
OBJECTS = eemumu.o fileini.o \
          pawini.o pawfil1.o pawfil2.o pawend.o \
          genzini.o genzfil.o genzend.o \
          save.o save2.o vegas.o setgen.o genera.o myaccept.o \
          treat.o ranf.o dangle.o

.f.o:
	$(FC) $(OPT) $<

mylpair: $(OBJECTS)
	@echo "   linking program ..."
	$(LO) -o mylpair.exe $(OBJECTS) $(CERNLIB) $(GENZLIB)
	@echo "   run program with 'mylpair.job' "

clean:
	rm -f *.exe *.o
