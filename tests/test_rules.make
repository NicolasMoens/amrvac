# This is a template Makefile to simplify writing tests
ifndef AMRVAC_DIR
$(error AMRVAC_DIR is not set)
endif

# Can be needed to compile the compare_log utility
ARCH ?= default
export ARCH

# Location of setup script
SETUP_SCRIPT := $(AMRVAC_DIR)/setup.pl

# Disable built-in rules
.SUFFIXES:

# Tool to compare test output numerically
LOG_CMP := $(AMRVAC_DIR)/tools/fortran/compare_logs

# Number of MPI processes to use
NUM_PROCS ?= 4

# force is a dummy to force re-running tests
.PHONY: all force

all: $(TESTS)

# Include architecture and compilation rules for the compare_log utility
include $(AMRVAC_DIR)/arch/$(ARCH).defs
include $(AMRVAC_DIR)/arch/rules.make

$(TESTS): $(LOG_CMP)

%.log: amrvac force
	@$(RM) $@		# Remove log to prevent pass when aborted
	@mpirun -np $(NUM_PROCS) ./amrvac -i $(PARS) > run.log
	@if $(LOG_CMP) 1.0e-5 1.0e-8 $@ correct_output/$@ ; \
	then echo "PASSED $@" ; \
	else echo "** FAILED ** $@" ; \
	fi

amrvac: makefile force		# Always try to build
	@$(MAKE)

makefile: $(AMRVAC_DIR)/arch/amrvac.make
	@$(RM) $@
	@$(SETUP_SCRIPT) $(SETUP_FLAGS) > setup.log

