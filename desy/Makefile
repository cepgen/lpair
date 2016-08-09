FF_FILES = $(wildcard source/*.f)
EX_FILES = $(wildcard external/*.f)
OBJ_DIR = obj
LIBRARY = -L/usr/lib64/cernlib/2006/lib #-lpdflib804
OBJ_FILES = $(patsubst source/%.f,$(OBJ_DIR)/%.o,$(FF_FILES))
EXT_FILES = $(patsubst external/%.f,$(OBJ_DIR)/%.o,$(EX_FILES))
LIB_FILES = $(EXT_FILES) $(OBJ_FILES)
VPATH = source/ external/

###############################################################################

FC = gfortran
FFLAGS = -g -w
CC = g++
CFLAGS = -lgfortran -Wall -I$(RHEAD) -lrt
RFLAGS = $(shell root-config --cflags)
RLIBS = $(shell root-config --libs)
RHEAD = $(shell root-config --incdir)

RM = /bin/rm
RMFLAGS = -rf

.PHONY: all
all: lpair

#
# Make the executable
#
lpair: $(OBJ_DIR)/main.opp $(LIB_FILES) $(OBJ_FILES)
	@echo "Linking..."
	@$(CC) $(CFLAGS) -o $@ $^ $(LIBRARY) $(RLIBS)
	@echo "Done!"

xsect: $(OBJ_DIR)/xsect.opp $(LIB_FILES) $(OBJ_FILES)
	@echo "Linking $<..."
	@$(CC) $(CFLAGS) -o $@ $^ $(LIBRARY)
	@echo "Done!"

nice:
	@$(RM) $(RMFLAGS) *.o $(OBJ_DIR)

clean: nice
	@$(RM) $(RMFLAGS) lpair
#
# Make the objects
#
$(OBJ_DIR)/%.o: %.f
	@echo "Building "$<
	@$(FC) -c $(FFLAGS) $< -o $@

$(OBJ_DIR)/%.opp: %.cpp | $(OBJ_DIR)
	@echo "Building "$<
	@$(CC) $(CFLAGS) -c $(RFLAGS) $< -o $@

$(OBJ_FILES): | $(OBJ_DIR)
$(EXT_FILES): | $(OBJ_DIR)

$(OBJ_DIR):
	@mkdir -p $(OBJ_DIR)
