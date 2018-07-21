FF_FILES = $(wildcard src/*.f)
EX_FILES = $(wildcard ../commons/*.f)
OBJ_DIR = obj
LIBRARY = #-L/usr/lib64/cernlib/2006/lib -lpdflib804 -lmathlib -lpacklib -lkernlib
OBJ_FILES = $(patsubst src/%.f,$(OBJ_DIR)/%.o,$(FF_FILES))
EXT_FILES = $(patsubst ../desy/external/%.f,$(OBJ_DIR)/%.o,$(EX_FILES))
LIB_FILES = $(EXT_FILES) $(OBJ_FILES)
VPATH = src/ external/

###############################################################################

FC = @gfortran
FFLAGS = -g -w -static -lc -lm

CC = @g++
CFLAGS = -lgfortran -Wall -lrt

RM = /bin/rm
RMFLAGS = -rf

# ROOT flags
RFLAGS = $(shell root-config --cflags)
RLIBS = $(shell root-config --libs)
RHEAD = $(shell root-config --incdir)

.PHONY: all
all: xsect clpair

#
# Make the executable
#

clpair: $(OBJ_DIR)/lpair.oxx $(LIB_FILES) $(OBJ_FILES) $(LIB_FILES)
	@echo "Linking $@..."
	$(CC) $(CFLAGS) -I$(RHEAD) -o $@ $^ $(LIBRARY) $(RLIBS)

nice:
	$(RM) $(RMFLAGS) *.o $(OBJ_DIR)

xsect: $(OBJ_DIR)/xsect.opp $(LIB_FILES) $(OBJ_FILES)
	@echo "Linking $@..."
	$(CC) $(CFLAGS) -o $@ $^ $(LIBRARY)

clean: nice
	$(RM) $(RMFLAGS) lpair xsect
#
# Make the objects
#
$(OBJ_DIR)/%.o: %.f
	@echo "Building "$<
	$(FC) -c $(FFLAGS) $< -o $@

$(OBJ_DIR)/%.opp: %.cpp | $(OBJ_DIR)
	@echo "Building "$<
	$(CC) -c $(CFLAGS) $< -o $@

$(OBJ_DIR)/%.oxx: %.cxx
	@echo "Building "$<
	$(CC) -c $(CFLAGS) $(RFLAGS) $< -o $@

$(OBJ_FILES): | $(OBJ_DIR)
$(EXT_FILES): | $(OBJ_DIR)

$(OBJ_DIR):
	@mkdir -p $(OBJ_DIR)
