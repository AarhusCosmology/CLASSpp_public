#Some Makefile for CLASS.
#Julien Lesgourgues, 28.11.2011

MDIR := $(shell pwd)
WRKDIR = $(MDIR)/build

.base:
	if ! [ -e $(WRKDIR) ]; then mkdir $(WRKDIR) ; mkdir $(WRKDIR)/lib; fi;
	touch build/.base

vpath %.h source:tools:main:include:species
vpath %.c source:tools:main
vpath %.cpp source:tools:main:species
vpath %.o build
vpath %.opp build
vpath .base build

########################################################
###### LINES TO ADAPT TO YOUR PLATEFORM ################
########################################################

# your C compiler:
CC        = gcc
CXX       = g++

# Optimization flag. The default uses -O3 which auto-vectorizes using the
# baseline SIMD for your architecture (SSE2 on x86-64, NEON on ARM64).
# On heterogeneous clusters where the build host may differ from the compute
# nodes, avoid -march=native to prevent illegal-instruction crashes.
# To enable AVX2+FMA (available on all x86 nodes since ~2015, 256-bit SIMD):
#   make SIMD_FLAGS="-mavx2 -mfma"
SIMD_FLAGS ?=
OPTFLAG = -O3 $(SIMD_FLAGS) #-ffast-math

# all other compilation flags
CCFLAG = -g -fPIC
CXXFLAG = $(CCFLAG) -std=c++17
LDFLAG = -g -fPIC
LIBRARIES = -lm -lpthread

# leave blank to compile without HyRec, or put path to HyRec directory
# (with no slash at the end: e.g. hyrec or ../hyrec)
HYREC = hyrec

########################################################
###### IN PRINCIPLE THE REST SHOULD BE LEFT UNCHANGED ##
########################################################

# pass current working directory to the code
CCFLAG += -D__CLASSDIR__='"$(MDIR)"'

# where to find include files *.h
INCLUDES = -I../include -I../tools -I../source -I../

# automatically add external programs if needed. First, initialize to blank.
EXTERNAL =

# eventually update flags for including HyRec
ifneq ($(HYREC),)
vpath %.c $(HYREC)
CCFLAG += -DHYREC
#LDFLAGS += -DHYREC
INCLUDES += -I../hyrec
EXTERNAL += hyrectools.o helium.o hydrogen.o history.o
endif
.SUFFIXES: .c .cpp .o .opp .h

# We could let gcc generate dependency information automatically, see this link:
# https://make.mad-scientist.net/papers/advanced-auto-dependency-generation/
# However, a clean build of CLASS is so fast that we just rebuild everything if *any*
# .h-file changed.
H_ALL = $(notdir $(wildcard include/*.h) $(wildcard tools/*.h) $(wildcard source/*.h) $(wildcard species/*.h))

%.o: %.c .base $(H_ALL)
	cd $(WRKDIR);$(CC) $(OPTFLAG) $(CCFLAG) $(INCLUDES) -c ../$< -o $*.o

%.opp: %.cpp .base $(H_ALL)
	cd $(WRKDIR); $(CXX) $(OPTFLAG) $(CXXFLAG) $(INCLUDES) -c ../$< -o $*.opp

TOOLS = growTable.opp dei_rkck.opp sparse.opp evolver_rkck.opp arrays.opp parser.opp quadrature.opp hyperspherical.opp common.opp trigonometric_integrals.opp non_cold_dark_matter.opp dark_radiation.opp exceptions.opp evolver_ndf15.opp

SPECIES_OPP = cdm.opp photons.opp baryons.opp lambda.opp ultra_relativistic.opp fluid.opp dcdm.opp dark_radiation_species.opp ncdm_species.opp scalar_field.opp interacting_species.opp composite_species.opp dcdm_dr_species.opp idm_dr_idr_species.opp idm_drmd_idr_drmd_species.opp dncdm_species.opp dncdm_decay_radiation_species.opp dncdm_dr_species.opp perturb_column_writer.opp

SOURCE = background_column_writer.opp input_module.opp background_module.opp thermodynamics_module.opp perturbations_module.opp primordial_module.opp nonlinear_module.opp transfer_module.opp spectra_module.opp lensing_module.opp cosmology.opp

OUTPUT = output_module.opp

CLASS = class.opp

all: class classy

# Use following line for faster installation if build dependencies (Cython, setuptools, numpy,..) are already installed:
# pip install --no-build-isolation 
# The Python one-liner below renames the compiled library to include the Python version number. This is neccessary
# due to bug introduced in MontePython 3.6, see https://github.com/brinckmann/montepython_public/issues/371

classy:
	rm -rf python/build && mkdir python/build
	pip install .
	cp -r build/lib* python/build
	python -c 'import sys; import os; import glob; \
file=glob.glob(os.path.join("python/build", "lib.*"))[0]; \
new_file=file.replace("lib.", "lib." + sys.version + "."); \
os.rename(file, new_file)'

class: $(TOOLS) $(SPECIES_OPP) $(SOURCE) $(EXTERNAL) $(OUTPUT) $(CLASS)
	$(CXX) $(OPTFLAG) $(LDFLAG) -o class $(addprefix build/,$(notdir $^)) $(LIBRARIES)

clean: .base
	rm -rf $(WRKDIR);
