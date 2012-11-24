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
      -fno-typeless-boz  -fugly-complex  \
      -funderscoring -c
OPT2 = $(OPT)
#
FC = g77
LO = g77
#
# Libraries:
#
CERNLIB = -L /usr/lib/cernlib/2006/lib \
          -lpdflib804 -lmathlib -lpacklib -lkernlib \
          -lpawlib -lgraflib -lgrafX11 -lphtools #-ljetset74
ATLSOFT = ../lib
GENZLIB = $(ATLSOFT)/libgenz.a
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
	$(LO) -o mylpair.exe $(OBJECTS) $(GENZLIB) $(CERNLIB)
	@echo "   run program with 'mylpair.job' "

clean:
	rm -f *.exe *.o
