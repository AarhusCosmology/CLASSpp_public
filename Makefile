#Some Makefile for CLASS.
#Julien Lesgourgues, 28.11.2011

MDIR := $(shell pwd)
WRKDIR = $(MDIR)/build

.base:
	if ! [ -e $(WRKDIR) ]; then mkdir $(WRKDIR) ; mkdir $(WRKDIR)/lib; fi;
	touch build/.base

vpath %.h source:tools:main:include
vpath %.c source:tools:main
vpath %.cpp source:tools:main
vpath %.o build
vpath %.opp build
vpath .base build

########################################################
###### LINES TO ADAPT TO YOUR PLATEFORM ################
########################################################

# your C compiler:
CC        = gcc
CXX       = g++

# your optimization flag
OPTFLAG = -O3 #-ffast-math -march=native

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
INCLUDES = -I../include -I../tools -I../source

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
H_ALL = $(notdir $(wildcard include/*.h) $(wildcard tools/*.h) $(wildcard source/*.h))

%.o: %.c .base $(H_ALL)
	cd $(WRKDIR);$(CC) $(OPTFLAG) $(CCFLAG) $(INCLUDES) -c ../$< -o $*.o

%.opp: %.cpp .base $(H_ALL)
	cd $(WRKDIR); $(CXX) $(OPTFLAG) $(CXXFLAG) $(INCLUDES) -c ../$< -o $*.opp

TOOLS_O = growTable.o dei_rkck.o sparse.o evolver_rkck.o arrays.o parser.opp quadrature.o hyperspherical.o common.o trigonometric_integrals.o

TOOLS_OPP = non_cold_dark_matter.opp dark_radiation.opp exceptions.opp evolver_ndf15.opp

TOOLS = $(TOOLS_O) $(TOOLS_OPP)

SOURCE = input_module.opp background_module.opp thermodynamics_module.opp perturbations_module.opp primordial_module.opp nonlinear_module.opp transfer_module.opp spectra_module.opp lensing_module.opp cosmology.opp

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

class: $(TOOLS) $(SOURCE) $(EXTERNAL) $(OUTPUT) $(CLASS)
	$(CXX) $(OPTFLAG) $(LDFLAG) -o class $(addprefix build/,$(notdir $^)) $(LIBRARIES)

clean: .base
	rm -rf $(WRKDIR);
