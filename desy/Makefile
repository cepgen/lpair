CERNLIB_PATH = /usr/lib64/cernlib/2006/lib

F_FILES = $(wildcard source/*.f)
OBJ_DIR = obj
OBJ_FILES = $(patsubst source/%.f,$(OBJ_DIR)/%.o,$(F_FILES))
LIBRARY = -L$(CERNLIB_PATH) \
          -lpdflib804 -lmathlib -lgeant321 -lpacklib -lkernlib \
          -lpawlib -lgraflib -lgrafX11 -lphtools
LIB_FILES = ilpair-cms-pp.o evertex.o $(OBJ_FILES)
VPATH = source/

###############################################################################

FC = gfortran
FFLAGS = -g -fno-automatic -w 

RM = /bin/rm
RMFLAGS = -rf

.PHONY: all
all: lpair

#
# Make the executable
#
lpair: $(LIB_FILES)
	$(FC) $(FFLAGS) -o $@ $^ $(LIBRARY)

nice:
	$(RM) $(RMFLAGS) *.o $(OBJ_DIR)

clean: nice
	$(RM) $(RMFLAGS) lpair
#
# Make the objects
#
$(OBJ_DIR)/%.o: %.f
	$(FC) -c $(FFLAGS) $< -o $@

$(OBJ_FILES): | $(OBJ_DIR)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)
