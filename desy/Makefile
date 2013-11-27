FF_FILES = $(wildcard source/*.f)
EX_FILES = $(wildcard external/*.f)
OBJ_DIR = obj
LIBRARY = -L/usr/lib64/cernlib/2006/lib -lpdflib804
OBJ_FILES = $(patsubst source/%.f,$(OBJ_DIR)/%.o,$(FF_FILES))
EXT_FILES = $(patsubst external/%.f,$(OBJ_DIR)/%.o,$(EX_FILES))
LIB_FILES = ilpair-cms-pp.o $(EXT_FILES) $(OBJ_FILES)
VPATH = source/ external/

###############################################################################

FC = @gfortran
FFLAGS = -g -w 

RM = /bin/rm
RMFLAGS = -rf

.PHONY: all
all: lpair

#
# Make the executable
#
lpair: $(LIB_FILES) $(OBJ_FILES) $(LIB_FILES)
	@echo "Linking..."
	$(FC) $(FFLAGS) -o $@ $^ $(LIBRARY)
	@echo "Done !"

nice:
	$(RM) $(RMFLAGS) *.o $(OBJ_DIR)

clean: nice
	$(RM) $(RMFLAGS) lpair
#
# Make the objects
#
$(OBJ_DIR)/%.o: %.f
	@echo "Building "$<
	$(FC) -c $(FFLAGS) $< -o $@

$(OBJ_FILES): | $(OBJ_DIR)
$(EXT_FILES): | $(OBJ_DIR)

$(OBJ_DIR):
	@mkdir -p $(OBJ_DIR)
