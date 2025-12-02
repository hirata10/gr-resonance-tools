# =========================
# Makefile for tidal resonance pipeline
# =========================

# Compiler and flags
CC     = gcc
CFLAGS = -w
LIBS   = -lm

# Directories
SRC_C   = src/C
SRC_PY  = src/python
OUTDIR  = td_res_exe
TEST_OUTDIR = test_exe

# Perl script and input
PARAMS_SCRIPT = pars.pl
PARAMS_INPUT  = pars.txt
C_HEADER      = globalpars_c.h
PY_FILE       = globalpars.py
TIMESTAMP := $(shell date +"%Y-%m-%d_%H-%M-%S")
DUMP_DIR := dump/dump-$(TIMESTAMP)

# Executables
EXES = $(OUTDIR)/J_evolve_single_start \
       $(OUTDIR)/J_evolve_single_restart \
       $(OUTDIR)/Delta_J_single \
       $(OUTDIR)/Delta_Phi_single

# =========================
# Default target
# =========================
all: td_res

td_res: $(OUTDIR) $(EXES)

# Create output directory
$(OUTDIR):
	mkdir -p $(OUTDIR)

# =========================
# Generate globals
# =========================
# Runs Perl script and moves generated files into proper directories
$(SRC_C)/$(C_HEADER) $(SRC_PY)/$(PY_FILE): $(PARAMS_SCRIPT) $(PARAMS_INPUT)
	perl $(PARAMS_SCRIPT)
	mv -f $(C_HEADER) $(SRC_C)/
	mv -f $(PY_FILE) $(SRC_PY)/

# =========================
# C Executables
# =========================

# J_evolve_single_norestart
$(OUTDIR)/J_evolve_single_start: $(SRC_C)/calling.c $(SRC_C)/resonance_find.c \
                                     $(SRC_C)/kerrtraj.c $(SRC_C)/kerrmode.c \
                                     $(SRC_C)/kerrgwem.c $(SRC_C)/J_dot.c \
                                     $(SRC_C)/J2J_dot.c $(SRC_C)/$(C_HEADER)
	$(CC) $(CFLAGS) -DIS_RK4_J_DOT_start \
        $(SRC_C)/calling.c $(SRC_C)/resonance_find.c \
        $(SRC_C)/kerrtraj.c $(SRC_C)/kerrmode.c \
        $(SRC_C)/kerrgwem.c $(SRC_C)/J_dot.c $(SRC_C)/J2J_dot.c \
        -I$(SRC_C) -o $@ $(LIBS)

# J_evolve_single_restart
$(OUTDIR)/J_evolve_single_restart: $(SRC_C)/calling.c $(SRC_C)/resonance_find.c \
                                   $(SRC_C)/kerrtraj.c $(SRC_C)/kerrmode.c \
                                   $(SRC_C)/kerrgwem.c $(SRC_C)/J_dot.c \
                                   $(SRC_C)/J2J_dot.c $(SRC_C)/$(C_HEADER)
	$(CC) $(CFLAGS) -DIS_RK4_J_DOT_restart \
        $(SRC_C)/calling.c $(SRC_C)/resonance_find.c \
        $(SRC_C)/kerrtraj.c $(SRC_C)/kerrmode.c \
        $(SRC_C)/kerrgwem.c $(SRC_C)/J_dot.c $(SRC_C)/J2J_dot.c \
        -I$(SRC_C) -o $@ $(LIBS)

# Delta_J_single
$(OUTDIR)/Delta_J_single: $(SRC_C)/calling.c $(SRC_C)/kerrtraj.c $(SRC_C)/kerrmode.c \
                          $(SRC_C)/kerrgwem.c $(SRC_C)/resonance_find.c \
                          $(SRC_C)/Gamma.c $(SRC_C)/Delta_J.c $(SRC_C)/J_dot.c \
                          $(SRC_C)/$(C_HEADER)
	$(CC) $(CFLAGS) $(SRC_C)/calling.c -DIS_DELTA_J_SINGLE \
        $(SRC_C)/kerrtraj.c $(SRC_C)/kerrmode.c $(SRC_C)/kerrgwem.c \
        $(SRC_C)/resonance_find.c $(SRC_C)/Gamma.c $(SRC_C)/Delta_J.c \
        $(SRC_C)/J_dot.c -I$(SRC_C) -o $@ $(LIBS)

# Delta_Phi_single
$(OUTDIR)/Delta_Phi_single: $(SRC_C)/kerrphase.c $(SRC_C)/kerrtraj.c \
                            $(SRC_C)/kerrgwem.c $(SRC_C)/kerrmode.c \
                            $(SRC_C)/resonance_find.c $(SRC_C)/J_dot.c \
                            $(SRC_C)/$(C_HEADER)
	$(CC) $(CFLAGS) $(SRC_C)/kerrphase.c $(SRC_C)/kerrtraj.c \
        $(SRC_C)/kerrgwem.c $(SRC_C)/kerrmode.c \
        $(SRC_C)/resonance_find.c $(SRC_C)/J_dot.c \
        -I$(SRC_C) -o $@ $(LIBS)

# ==================================================
# =============  DEBUG SECTION  ===============
# ==================================================

# You can change these compiler options:
TEST_CFLAGS = -Wall    
TEST_LIBS   = -lm

# Test executables
TEST_EXES = $(TEST_OUTDIR)/J2EQL \
            $(TEST_OUTDIR)/J2J_DOT_TIDAL \
            $(TEST_OUTDIR)/J2J_DOT_SF 

# Running "make test" builds all test executables
test: $(TEST_OUTDIR) $(TEST_EXES)

# make directory for test executables
$(TEST_OUTDIR):
	mkdir -p $(TEST_OUTDIR)

# Create executable to compute a set of Js to EQL and orbit data
$(TEST_OUTDIR)/J2EQL: $(SRC_C)/kerrphase.c $(SRC_C)/kerrtraj.c \
                            $(SRC_C)/kerrgwem.c $(SRC_C)/kerrmode.c \
                            $(SRC_C)/resonance_find.c $(SRC_C)/J_dot.c \
                            $(SRC_C)/$(C_HEADER)
	$(CC) $(TEST_CFLAGS) -DIS_J2EQL $(SRC_C)/calling.c $(SRC_C)/kerrtraj.c \
        $(SRC_C)/kerrgwem.c $(SRC_C)/kerrmode.c \
        -I$(SRC_C) -o $@ $(TEST_LIBS)

# Create executable to compute a set of J_dot_tidal of inner body from Js (inner and outer) and resonance mode
$(TEST_OUTDIR)/J2J_DOT_TIDAL: $(SRC_C)/kerrphase.c $(SRC_C)/kerrtraj.c \
                            $(SRC_C)/kerrgwem.c $(SRC_C)/kerrmode.c \
                            $(SRC_C)/resonance_find.c $(SRC_C)/J_dot.c \
                            $(SRC_C)/$(C_HEADER)
	$(CC) $(TEST_CFLAGS) -DIS_J_DOT_TIDAL $(SRC_C)/calling.c $(SRC_C)/kerrtraj.c \
        $(SRC_C)/kerrgwem.c $(SRC_C)/kerrmode.c $(SRC_C)/resonance_find.c $(SRC_C)/J_dot.c \
        -I$(SRC_C) -o $@ $(TEST_LIBS)

# Create executable to compute a set of J_dot_selfforce of a body
$(TEST_OUTDIR)/J2J_DOT_SF: $(SRC_C)/kerrphase.c $(SRC_C)/kerrtraj.c \
                            $(SRC_C)/kerrgwem.c $(SRC_C)/kerrmode.c \
                            $(SRC_C)/resonance_find.c $(SRC_C)/J_dot.c \
                            $(SRC_C)/$(C_HEADER)
	$(CC) $(TEST_CFLAGS) -DIS_J_DOT_SF $(SRC_C)/calling.c $(SRC_C)/kerrtraj.c \
        $(SRC_C)/kerrgwem.c $(SRC_C)/kerrmode.c $(SRC_C)/resonance_find.c $(SRC_C)/J_dot.c \
        -I$(SRC_C) -o $@ $(TEST_LIBS)


# =========================
# Clean Targets
# =========================

# If only editing C code for td_res, remove only td_res build executables
clean-td_res:
	rm -rf $(OUTDIR)

# If only editing C code for test, remove only test build executables
clean-test:
	rm -rf $(TEST_OUTDIR)

# If making changes to python code with existing data directories, move the existing directories to a dated dump directory
clean-data:
	@found=false; \
	# List of all files/directories to archive
	for d in outputs_data Output_Delta_J_* Output_Delta_Phi_* \
	         $(SRC_C)/$(C_HEADER) $(SRC_PY)/$(PY_FILE) \
	         action_angle_pairs.txt precomputed_omegas_* potential_resonances_*; do \
		if [ -e "$$d" ]; then found=true; break; fi; \
	done; \
	if [ "$$found" = true ]; then \
		echo "Archiving Python outputs, global pars files, and additional data into $(DUMP_DIR)"; \
		mkdir -p "$(DUMP_DIR)"; \
		for d in outputs_data Output_Delta_J_* Output_Delta_Phi_* \
		         $(SRC_C)/$(C_HEADER) $(SRC_PY)/$(PY_FILE) \
		         action_angle_pairs.txt precomputed_omegas_* potential_resonances_*; do \
			if [ -e "$$d" ]; then \
				mv -f "$$d" "$(DUMP_DIR)/" || true; \
			fi; \
		done; \
	else \
		echo "No Python outputs, global pars files, or additional data to archive."; \
	fi



# If making changes to the pars.txt, then clean everything and start fresh
clean: clean-td_res clean-test clean-data
	rm -f $(SRC_C)/$(C_HEADER) $(SRC_PY)/$(PY_FILE)

.PHONY: all td_res clean clean-td_res clean-test clean-data clean-py test