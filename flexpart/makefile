SHELL = /bin/bash
################################################################################
#  DESCRIPTION
#    Makefile for FLEXPART.
#
#    minimum requirement:
#     - ECCODES library
#    optional requirements:
#     - NetCDF library
#    environment variables:
#      - include (mod): CPATH
#      - library: LIBRARY_PATH
#    Note: CPATH and LIBRARY path can also be hardcoded below.
#
#  USAGE
#    Compile FLEXPART, add optional V=1 to see include and library paths
#      make [-j] 
#    Compile FLEXPART with eta coordinates
#      make [-j] eta=yes
#    Compile FLEXPART without NetCDF
#      make [-j] eta=<yes/no> ncf=no
#    Compile with shared libraries (LD_LIBRARY_PATH required at runtime)
#     make [-j] shared=yes
#    Compile with debugging flags (no optimization, more checks, debug symbols)
#     make [-j] debug=yes
#    Compile a serial version (no OpenMP parallelism)
#     make [-j] serial=yes
#    Clean object files
#     make clean
#    Clean object files and executables
#     make cleanall
#
################################################################################
MAKEFLAGS += -j2

## PROGRAMS
# Unified executable names
# The same executable is used for both ECMWF and GFS metdata

# meter coordinates + NetCDF executable
FLEXPART = FLEXPART

# meter coordinates + Binary executable
FLEXPART-BIN = FLEXPART_BIN

# Eta coordinates netcdf executable
FLEXPART-ETA = FLEXPART_ETA

# Eta coordinates binary executable
FLEXPART-ETA-BIN = FLEXPART_ETA_BIN

export OMP_NESTED=FALSE

# as default look for GNU Fortran compiler
# override by setting FC in environment
# e.g. FC=ifort make
ifeq ($(FC), f77)
FC := gfortran
endif
FC ?= gfortran
# ensure compiler exists in PATH
ifeq ($(shell command -v $(FC) 2>/dev/null),)
$(error Fortran compiler '$(FC)' not found in PATH. Specify FC=)
endif
################################################################################
#
## Environment Modules set CPATH and LIBRARY_PATH
#
# You can also hardcode the include and library paths here
CPATH = /srv/conda/envs/notebook/include
LIBRARY_PATH = /srv/conda/envs/notebook/lib
#
ifndef CPATH
# default path for system headers and gfortran modules
CPATH = /usr/include:/usr/local/include:/lib64/gfortran/modules/
$(warning CPATH is not set. Set CPATH to include paths for required headers.)
$(warning Using default: $(CPATH))
endif
ifndef LIBRARY_PATH
# default path for system libraries and gfortran libraries
LIBRARY_PATH = /usr/lib64:/usr/local/lib64:/usr/lib:/usr/local/lib
$(warning LIBRARY_PATH is not set. Set LIBRARY_PATH to include paths for required libraries.)
$(warning Using default: $(LIBRARY_PATH))
endif
# use the environmental variable $INCLUDE
INC = $(subst :, ,$(CPATH))
INC := $(INC:%=-I%)  # add a -I/path/to/include
# use the environmental variable $LIBRARY_PATH
LIBPATH := $(subst :, ,$(LIBRARY_PATH))
ifeq ($(filter yes true,$(shared) $(SHARED)),)
# shared=yes ( default is shared=no )
shared := false
LIBRPATH := $(LIBPATH:%=-Wl,-rpath=%)  # we use rpath to link libraries at built time.
endif
LIBPATH := $(LIBPATH:%=-L%)
################################################################################
# Enable netCDF output?
# default is no netCDF output
ifeq ($(filter no false,$(ncf) $(NCF)),)
NCOPT = -DUSE_NCF
NCF := true
LIBS = -leccodes -leccodes_f90 -lm -lnetcdff
else
NCOPT = -UUSE_NCF
NCF := false
LIBS = -leccodes -leccodes_f90 -lm
endif
################################################################################
# Enable eta coordinates?
# default is no eta coordinates
ifeq ($(filter no false,$(eta) $(ETA)),)
ETAOPT = -DETA
ETA := true
else
ETAOPT = -UETA
ETA := false
endif
################################################################################
## OPTIMIZATION LEVEL
O_LEV = 2 # [0,1,2,3,g,s,fast]
# Allow to compile in serial or parallel mode using OpenMP
# e.g. SERIAL=yes make --> OpenMP disabled
ifneq ($(filter yes true,$(SERIAL) $(serial)),)
omp = # -fopenmp_stubs  This flag does not exist with gfortran?
SERIAL := true
else
omp = -fopenmp
SERIAL := false
endif
################################################################################
# Different CPU microarchitectures
# check https://gcc.gnu.org/onlinedocs/gcc/x86-Options.html
# check https://www.intel.com/content/www/us/en/docs/dpcpp-cpp-compiler/developer-guide-reference/2023-1/x-qx.html
# using core-avx2 (AVX2)
# using x86-64 (no vectorization)
# using skylake-avx512 (AVX512)
# using native is the default, which will use whatever it can
# (defined by the machine where you are compiling it)
ifdef arch
FUSER =  -march=$(arch)
else
FUSER =  -march=native -mtune=native
arch = native
endif
################################################################################
#
## Build flags for different compilers
#
ifneq ($(filter %gfortran, $(FC)),)
$(info Using GNU Fortran compiler: $(FC))
# -m64 is needed for large memory models (64 bit addressing)
# -cpp C preprocessor (for ifdefs)
# -convert=little_endian for little endian architectures (little endian is default on x86)
# -record-marker=4 for fixed format files with 4 byte record markers (legacy)
# -mcmodel=large for large memory models (addresses larger than 2GB)
FFLAGS = -O$(O_LEV) -m64 -cpp -mcmodel=large -fconvert=little-endian -frecord-marker=4 $(NCOPT) $(ETAOPT) $(FUSER) $(omp)
# DEBUG FLAGS
# -ffpe-trap=invalid,overflow,denormal,underflow,zero  -Warray-bounds -fcheck=all
# -fbacktrace -fcheck=all -fbounds-check -fdefault-real-8 -freal-8-real-4 -pg
# -O0 no optimization
# -g3 maximum debug information
DBGFLAGS = -O0 -g3 -ggdb3 -m64 -cpp -mcmodel=large -fcheck=all -fconvert=little-endian -frecord-marker=4 -fmessage-length=0 -flto -fbacktrace -Wall -fdump-core $(NCOPT) $(ETAOPT) $(FUSER) $(omp)
else ifneq ($(filter %ifort %ifx, $(FC)),)
$(info Using Intel Fortran compiler: $(FC))
# -m64 is needed for large memory models (64 bit addressing)
# -cpp c preprocessor (for ifdefs)
# -convert little_endian for little endian architectures (little endian is default on x86)
# -assume byterecl for fixed format files with 4 byte record markers (legacy)
# -mcmodel large for large memory models (addresses larger than 2GB)
# -fpmodel fast/precise/strict/source for fast floating point operations (less precise)
## Note: -record-marker=4 is a gfortran option and is not supported by Intel compilers (ifort/ifx).
## For Intel compilers we use -assume byterecl (handles fixed-record files) and omit -record-marker.
FFLAGS = -O$(O_LEV) -m64 -fpp -mcmodel=large -convert little_endian -diag-disable=10448 $(NCOPT) $(ETAOPT) $(FUSER) $(omp)
# DEBUG FLAGS
# -check bounds for debugging array bounds
# -warn all to enable all warnings
# -warn error to treat warnings as errors
# -O0 no optimization
# -g3 maximum debug information
# -traceback to get a backtrace when the program crashes
DBGFLAGS = -O0 -g3 -fpp -m64 -convert little_endian -assume byterecl -mcmodel=large -check bounds -warn all -traceback $(NCOPT) $(ETAOPT) $(FUSER) $(omp)
else
# unknown compiler, use default flags
$(warning WARNING - Unknown Fortran compiler $(FC), using default flags (no optimization))
FFLAGS = $(NCOPT) $(ETAOPT) $(FUSER) $(omp)
DBGFLAGS = -g $(FFLAGS)
endif
# Libraries and flags for linking
# -Wl,-rpath, adds the path to the library to the runtime search path
LDFLAGS = $(LIBS) $(LIBRPATH)
# Debug linking flags
LDDEBUG = $(LIBS)
################################################################################
# Enable debugging?
ifeq ($(filter yes true,$(DEBUG) $(debug)),)
DEBUG := false
else
FFLAGS := $(DBGFLAGS)
LDFLAGS := $(LDDEBUG)
DEBUG := true
endif
ifeq ($(origin V), undefined)
Q := @  # suppress command echo
else
Q :=	# print commands
endif
################################################################################
# Print detected compiler and flags at make startup
ifneq ($(filter help,$(MAKECMDGOALS)), help)
$(info -------------------------------)
$(info Detected Fortran compiler: $(notdir $(FC))  (default: FC=gfortran))
$(info Using architecture flags:  $(arch)  (default: arch=native))
$(info Using debugging flags:     $(DEBUG)  (default: debug=no))
$(info NetCDF output enabled:     $(NCF)  (default: ncf=yes))
$(info Eta coordinates enabled:   $(ETA)  (default: eta=yes))
$(info Serial (OpenMP):           $(SERIAL)  (default: serial=no))
$(info Shared libraries:          $(shared)  (default: shared=no))
$(info FFLAGS: $(FFLAGS))
$(info LDFLAGS: $(LDFLAGS))
ifneq ($(origin V), undefined)
$(info Includes (INC): $(INC))
$(info Library path flags (LIBPATH): $(LIBPATH))
else
$(info To see include and library paths, run make with the V=1 flag)
endif
$(info -------------------------------)
endif

################################################################################
MODOBJS = \
advance_mod.o 		binary_output_mod.o \
cbl_mod.o 		cmapf_mod.o \
com_mod.o 		conv_mod.o \
class_gribfile_mod.o \
date_mod.o 		drydepo_mod.o \
erf_mod.o 		flux_mod.o \
getfields_mod.o 	initialise_mod.o \
interpol_mod.o 		mean_mod.o \
outgrid_mod.o \
output_mod.o 		par_mod.o \
particle_mod.o 		pbl_profile_mod.o \
plume_mod.o 		point_mod.o \
random_mod.o 		readoptions_mod.o \
restart_mod.o 		settling_mod.o \
sort_mod.o \
timemanager_mod.o 	turbulence_mod.o \
txt_output_mod.o 	unc_mod.o \
verttransform_mod.o 	wetdepo_mod.o \
windfields_mod.o 	qvsat_mod.o \
xmass_mod.o		receptor_mod.o 

OBJECTS_NCF = \
netcdf_output_mod.o	chemistry_mod.o \
totals_mod.o		initdomain_mod.o \
receptor_netcdf_mod.o	emissions_mod.o \

OBJECTS_ETA = coord_ecmwf_mod.o 

OBJECTS = \
FLEXPART.o

ifeq ($(ncf), no)
OBJECTS := $(OBJECTS)
ifeq ($(eta), no)
FP := $(FLEXPART-BIN)
else
FP := $(FLEXPART-ETA-BIN)
OBJECTS := $(OBJECTS) $(OBJECTS_ETA)
endif
else
OBJECTS := $(OBJECTS) $(OBJECTS_NCF)
ifeq ($(eta), no)
FP := $(FLEXPART)
else
FP := $(FLEXPART-ETA)
OBJECTS := $(OBJECTS) $(OBJECTS_ETA)
endif
endif
# Print object files (used for this build)
$(info Building: $(FP))

# make runs the first target when non is given!
.PHONY: help clean cleanall recompile gitversion
# First target
.DEFAULT_GOAL := $(FP)
clean:
	\rm -f *.o *.mod
cleanall:
	\rm -f *.o *.mod $(FLEXPART) $(FLEXPART-BIN) $(FLEXPART-ETA) $(FLEXPART-ETA-BIN)

help:
	@echo "Usage: make [options]"
	@echo "Options:"
	@echo "  eta=yes|no        Compile with eta coordinates (default: yes)"
	@echo "  ncf=yes|no        Compile with NetCDF output (default: yes)"
	@echo "  debug=yes|no      Compile with debugging flags (default: no)"
	@echo "  serial=yes|no     Compile a serial version without OpenMP (default: no)"
	@echo "  shared=yes|no     Compile with shared libraries (default: no)"
	@echo "  arch=ARCH         Compile for specific architecture (default: native)"
	@echo "  FC=COMPILER       Use specific Fortran compiler (default: gfortran)"
	@echo "  V=1               Print include and library paths"

#
# Rules to build the program
%.o: %.mod
%.o: %.f90
	+$(Q)$(FC) -c $(FFLAGS) $(INC) $<
# FLEXPART EXE target
$(FP): $(MODOBJS) $(OBJECTS)
	+$(FC) -o $@ $(MODOBJS) $(OBJECTS) $(FFLAGS) $(LIBPATH) $(LDFLAGS)

.SUFFIXES = $(SUFFIXES) .f90
#
##	ETA RECOMPILE
#
etadep := no
ifeq ($(eta), no)
ifneq ("$(wildcard ./coord_ecmwf_mod.o)","")
etadep := "yes"
else
etadep := "no"
endif
else
ifneq ("$(wildcard ./coord_ecmwf_mod.o)","")
etadep := "no"
else
etadep := "yes"
endif
endif
#
# target for forced recompilation if eta dependency changes
#
recompile:
	@echo Why, $(etadep)

#
## GIT VERSIONING
#
# Check if git command is successful (there is a .git directory)
ifeq ($(shell git rev-parse --is-inside-work-tree 2>/dev/null), true)
GIT_VERSION := "$(shell git log --pretty=format:"%h %ad" -n 1)"
else
GIT_VERSION := "undefined"
endif

# Target to update git version in FLEXPART.f90
gitversion:
	@sed -i '0,/gitversion_tmp/s;gitversion_tmp.*;gitversion_tmp=$(GIT_VERSION);' FLEXPART.f90;
	$(info Git version: $(GIT_VERSION))
#
## DEPENDENCIES
#

#1) Independent modules: par_mod.o, com_mod.o, qvsat_mod.o class_gribfile_mod.o
#2) Modules that are only dependent on independent modules
cmapf_mod.o: par_mod.o
com_mod.o: par_mod.o
date_mod.o: par_mod.o
erf_mod.o: par_mod.o
mean_mod.o: par_mod.o
sort_mod.o: par_mod.o
ifeq ($(etadep), "yes")
particle_mod.o: recompile com_mod.o par_mod.o
else
particle_mod.o: com_mod.o par_mod.o
endif
point_mod.o: com_mod.o par_mod.o
unc_mod.o: par_mod.o com_mod.o
cbl_mod.o: com_mod.o par_mod.o random_mod.o
pbl_profile_mod.o: par_mod.o qvsat_mod.o

#3) 3rd level dependencies
turbulence_mod.o: particle_mod.o cbl_mod.o qvsat_mod.o
ifeq ($(etadep), "yes")
windfields_mod.o: recompile cmapf_mod.o point_mod.o class_gribfile_mod.o date_mod.o qvsat_mod.o \
	pbl_profile_mod.o
else
windfields_mod.o: cmapf_mod.o point_mod.o class_gribfile_mod.o date_mod.o qvsat_mod.o \
	pbl_profile_mod.o
endif
#4) 4th level dependencies
verttransform_mod.o: par_mod.o com_mod.o qvsat_mod.o windfields_mod.o
settling_mod.o: windfields_mod.o
interpol_mod.o: windfields_mod.o particle_mod.o turbulence_mod.o

#5)
ifneq ($(eta),no)
coord_ecmwf_mod.o: interpol_mod.o
endif
drydepo_mod.o: unc_mod.o interpol_mod.o erf_mod.o

#6)
ifneq ($(eta),no)
advance_mod.o: coord_ecmwf_mod.o settling_mod.o drydepo_mod.o turbulence_mod.o
plume_mod.o: mean_mod.o coord_ecmwf_mod.o
outgrid_mod.o: unc_mod.o coord_ecmwf_mod.o
wetdepo_mod.o: coord_ecmwf_mod.o unc_mod.o
else
advance_mod.o: interpol_mod.o settling_mod.o drydepo_mod.o turbulence_mod.o
plume_mod.o: mean_mod.o particle_mod.o windfields_mod.o
outgrid_mod.o: unc_mod.o interpol_mod.o
wetdepo_mod.o: interpol_mod.o unc_mod.o
endif
readoptions_mod.o: drydepo_mod.o xmass_mod.o

#7)
getfields_mod.o: wetdepo_mod.o verttransform_mod.o qvsat_mod.o
flux_mod.o: outgrid_mod.o qvsat_mod.o
txt_output_mod.o: outgrid_mod.o
binary_output_mod.o: mean_mod.o outgrid_mod.o
ifneq ($(ncf), no)
receptor_mod.o: par_mod.o com_mod.o point_mod.o particle_mod.o date_mod.o \
        windfields_mod.o receptor_netcdf_mod.o binary_output_mod.o
else
receptor_mod.o: par_mod.o com_mod.o point_mod.o particle_mod.o date_mod.o \
        windfields_mod.o binary_output_mod.o
endif 
ifneq ($(ncf), no)
netcdf_output_mod.o: mean_mod.o outgrid_mod.o readoptions_mod.o drydepo_mod.o
totals_mod.o: par_mod.o com_mod.o netcdf_output_mod.o
chemistry_mod.o: par_mod.o com_mod.o date_mod.o particle_mod.o \
	point_mod.o windfields_mod.o totals_mod.o netcdf_output_mod.o
initdomain_mod.o: par_mod.o com_mod.o point_mod.o random_mod.o outgrid_mod.o \
	initialise_mod.o totals_mod.o date_mod.o windfields_mod.o netcdf_output_mod.o
receptor_netcdf_mod.o: par_mod.o com_mod.o point_mod.o date_mod.o \
	windfields_mod.o netcdf_output_mod.o 
emissions_mod.o: par_mod.o com_mod.o point_mod.o particle_mod.o date_mod.o \
	netcdf_output_mod.o totals_mod.o windfields_mod.o
endif

#8)
conv_mod.o: flux_mod.o class_gribfile_mod.o qvsat_mod.o sort_mod.o
ifeq ($(ncf), no)
output_mod.o: binary_output_mod.o txt_output_mod.o
ifneq ($(eta),no)
restart_mod.o: coord_ecmwf_mod.o unc_mod.o outgrid_mod.o
else
restart_mod.o: interpol_mod.o unc_mod.o outgrid_mod.o
endif
else
output_mod.o: netcdf_output_mod.o binary_output_mod.o txt_output_mod.o
ifneq ($(eta),no)
restart_mod.o: coord_ecmwf_mod.o unc_mod.o netcdf_output_mod.o
else
restart_mod.o: interpol_mod.o unc_mod.o netcdf_output_mod.o
endif
endif

#9)
initialise_mod.o: interpol_mod.o xmass_mod.o output_mod.o turbulence_mod.o

#10)
ifneq ($(ncf), no)
timemanager_mod.o: advance_mod.o conv_mod.o plume_mod.o getfields_mod.o restart_mod.o \
	initialise_mod.o chemistry_mod.o initdomain_mod.o receptor_mod.o \
	emissions_mod.o totals_mod.o
else
timemanager_mod.o: advance_mod.o conv_mod.o plume_mod.o getfields_mod.o restart_mod.o \
        initialise_mod.o receptor_mod.o
endif


#11) Final target
FLEXPART.o: timemanager_mod.o gitversion
